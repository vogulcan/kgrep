#!/usr/bin/env python3
"""
kgrep.py — Fast positional k-mer grep over RocksDict/RocksDB.

This script scans a RocksDict/RocksDB where:
- Keys are fixed-length k-mers (bytes/ASCII).
- Values are counts (big-endian u32, little-endian u64, or ASCII int).

It applies position-specific rules (via --mask or --rule) and performs
bounded, prefix-guided scans to minimize IO, using iterate_lower_bound /
iterate_upper_bound and prefix_same_as_start when available.

Default alphabet : DNA = b"ACGT"

Typical uses
------------
- Mask syntax with count filter: (See below on masking rules.)
    uv run main.py --db /path/dbfolder --mask "A[CT]?AA" \
        --values-are-u32-be --min-count 10

- Rule syntax with output file:
    uv run main.py --db /path/dbfolder --values-are-u32-be --rule "1:A;2:CT" \
               --min-count 1 --output output.txt

- Basic k-mer search without constraints: (Lists all)
    uv run main.py --db /path/dbfolder --k 10

- Complex rule with multiple positions:
    uv run main.py --db /path/dbfolder --rule "3:A;5:CG;10:T" \
               --values-are-u64 --limit 1000

Mask syntax table
=================

Symbol     Meaning
---------- -----------------------------------------------------------------
?          Any character from the alphabet (wildcard)
[AC]       Allowed set at this position (only listed characters are allowed)
G          Exact match for a single required character (here: 'G')

=================

Design notes
------------
- We enumerate only as many prefix combinations as needed (bounded by --max-enum).
- For each enumerated prefix, we seek to a tight (lo, hi) range and iterate.
- Within that range, we still verify full per-position constraints in-memory.
- Value decoding supports multiple encodings; ASCII fallback is tolerant.

Caveats
-------
- If keys are not raw bytes / ASCII kmers, pass --no-raw.
- If your DB uses a prefix extractor, prefix_same_as_start helps limit iteration.
"""

from __future__ import annotations

import argparse
import itertools
import sys
from typing import List, Dict, Iterable, Tuple, Optional

from rocksdict import Rdict, Options, ReadOptions

# Default alphabet for DNA. You can override with --alphabet.
DNA = b"ACGT"


# ---------- Portable helpers over rocksdict API differences ----------
def _try_call(obj, candidates: Iterable[str], *args) -> bool:
    """
    Attempt to call the first available method from `candidates` on `obj`.

    This smooths over version/API diffs between rocksdict bindings:
    some environments expose `set_iterate_lower_bound`, others `iterate_lower_bound`, etc.

    Returns:
        True if a candidate method exists and was invoked, else False.
    """
    for name in candidates:
        fn = getattr(obj, name, None)
        if callable(fn):
            fn(*args)
            return True
    return False


def _set_bound(ro: ReadOptions, kind: str, key: bytes) -> None:
    """
    Set iterate lower/upper bound in a portable way.

    Args:
        ro: ReadOptions
        kind: "lower" or "upper"
        key: boundary key (inclusive lower, exclusive upper in RocksDB semantics)
    """
    if kind == "lower":
        if _try_call(ro, ["set_iterate_lower_bound", "iterate_lower_bound"], key):
            return
    else:
        if _try_call(ro, ["set_iterate_upper_bound", "iterate_upper_bound"], key):
            return
    # If neither bound setter exists, we proceed without bounds.
    # NOTE: This is unusual; performance will degrade but correctness remains.


# ---------- Rule parsing ----------
def parse_rules(args: argparse.Namespace, k: int, alphabet: bytes) -> Dict[int, bytes]:
    """
    Parse rule specifications from --mask or --rule into a dict: pos -> allowed bytes.

    Priority:
        --mask takes precedence if provided; otherwise --rule.

    Mask syntax:
        - '?'    : any char
        - '[AC]' : allowed set at that position
        - 'G'    : literal char constraint
      Positions are 1-based in the mask string.

    Rule syntax:
        "3:A;5:CG;10:T"  -> position: allowed set

    Args:
        args: parsed CLI args
        k: k-mer length (for clipping; out-of-range positions ignored)
        alphabet: default allowed set for unconstrained positions

    Returns:
        Dict mapping 1-based position -> sorted bytes of allowed characters.
    """
    if args.mask:
        pos = 1
        i = 0
        s = args.mask
        rules: Dict[int, bytes] = {}
        while i < len(s):
            c = s[i]
            if c == '?':
                pos += 1
                i += 1
                continue
            if c == '[':
                j = s.find(']', i + 1)
                if j == -1:
                    raise ValueError("Unclosed '[' in --mask")
                allowed = s[i + 1 : j].encode('ascii').upper()
                rules[pos] = bytes(sorted(set(allowed)))
                pos += 1
                i = j + 1
                continue
            # Literal char
            rules[pos] = bytes([ord(c.upper())])
            pos += 1
            i += 1
        return rules

    # --rule pathway
    rules: Dict[int, bytes] = {}
    if args.rule:
        for term in args.rule.split(';'):
            term = term.strip()
            if not term:
                continue
            p_str, vals = term.split(':', 1)
            p = int(p_str)
            allowed = bytes(sorted(set(vals.upper().encode('ascii'))))
            rules[p] = allowed
    return rules


# ---------- Prefix planning ----------
def choose_prefix_len(
    rules: Dict[int, bytes],
    k: int,
    max_enum: int,
    alphabet: bytes,
) -> int:
    """
    Choose how many leading positions (P) to enumerate as explicit prefixes.

    We accumulate the cartesian product size as we move left-to-right. Stop when:
      - enumerations would exceed --max-enum, OR
      - we covered the last constrained position (early exit), OR
      - we reached k.

    This balances enumeration cost vs scan tightness.

    Returns:
        P in [0..k]. P=0 means no explicit prefix enumeration.
    """
    if not rules:
        return 0

    last_constrained = max((p for p in rules.keys() if 1 <= p <= k), default=0)
    P = 0
    combos = 1
    for pos in range(1, k + 1):
        P = pos
        card = len(rules.get(pos, alphabet))  # constrained set or full alphabet
        combos *= card
        if combos > max_enum:
            return P - 1  # previous P fits within the enumeration budget
        if pos >= last_constrained:
            # We already pinned every constrained position in the prefix;
            # no need to go further.
            return P
    return P


def enumerated_prefixes(P: int, rules: Dict[int, bytes], alphabet: bytes) -> Iterable[bytes]:
    """
    Yield all prefixes of length P consistent with the given rules.

    For P=0, yields a single empty prefix (i.e., full keyspace).
    """
    if P <= 0:
        yield b""
        return
    slots = [rules.get(pos, alphabet) for pos in range(1, P + 1)]
    for tup in itertools.product(*slots):
        yield bytes(tup)


# ---------- Fast checking ----------
def build_allowed_arrays(k: int, rules: Dict[int, bytes]) -> List[Optional[bytes]]:
    """
    Build a random-access array of allowed sets per 1-based position.

    Index 0 is unused to allow 1-based indexing: arr[pos] -> allowed bytes or None.
    """
    arr: List[Optional[bytes]] = [None] * (k + 1)
    for pos, allowed in rules.items():
        if 1 <= pos <= k:
            arr[pos] = allowed
    return arr


def match_kmer_bytes(kmer: bytes, allowed_by_pos: List[Optional[bytes]]) -> bool:
    """
    Check whether a k-mer (bytes) satisfies per-position allowed sets.

    Args:
        kmer: key bytes of length k
        allowed_by_pos: list where each entry is None (no constraint)
                        or a bytes object containing the allowed set.

    Returns:
        True if all constrained positions pass, else False.
    """
    for i, allowed in enumerate(allowed_by_pos[1:], start=0):  # i: 0-based index into kmer
        if allowed is None:
            continue
        if kmer[i] not in allowed:
            return False
    return True


# ---------- Bounds / range calculation ----------
def upper_bound_for_prefix(prefix: bytes) -> bytes:
    """
    Compute the exclusive upper bound key for a given prefix.

    RocksDB iterators treat upper bounds as exclusive. We increment the
    last non-0xFF byte to get the smallest key strictly greater than any
    key starting with `prefix`. If all bytes are 0xFF, append a NUL.
    """
    p = bytearray(prefix)
    for i in reversed(range(len(p))):
        if p[i] != 0xFF:
            p[i] += 1
            return bytes(p[: i + 1])
    return prefix + b"\x00"


def seek_range_for_prefix(prefix: bytes, k: int) -> Tuple[bytes, bytes]:
    """
    For a given prefix, derive a (lo, hi) that tightly bounds all keys
    with that prefix and total length k.

    lo: prefix + NUL padding up to k (inclusive lower bound)
    hi: exclusive upper bound computed from prefix
    """
    lo = prefix + b"\x00" * (k - len(prefix))
    hi = upper_bound_for_prefix(prefix)
    return lo, hi


# ---------- Main entrypoint ----------
def main() -> None:
    ap = argparse.ArgumentParser(description="k-mer rule grep over RocksDict")
    ap.add_argument("--db", required=True, help="Path to RocksDict / RocksDB")
    ap.add_argument("--cf", default=None, help="Column family name (default CF if omitted)")
    ap.add_argument("--k", type=int, help="k-mer length in bytes (auto from first key if omitted)")
    ap.add_argument("--rule", help='Rules like "3:A;5:CG;10:T"', default=None)
    ap.add_argument("--mask", help='Compact mask like "A[CT]?G??T" (?=any, [..]=set)', default=None)
    ap.add_argument("--alphabet", default="ACGT", help="Alphabet, default ACGT")
    ap.add_argument("--min-count", type=int, default=0, help="Only output k-mers with count >= this")
    ap.add_argument("--max-enum", type=int, default=200000, help="Max prefix combinations to enumerate")
    ap.add_argument("--limit", type=int, default=0, help="Stop after N matches (0 = unlimited)")
    ap.add_argument(
        "--values-are-u32-be",
        action="store_true",
        help="Interpret value as big-endian u32 count (builder used struct.pack('>I', cnt))",
    )
    ap.add_argument("--values-are-u64", action="store_true", help="Interpret value as little-endian u64 count")
    ap.add_argument("--no-raw", action="store_true", help="Open without raw_mode (use if your DB isn’t raw bytes)")
    ap.add_argument("--readahead", type=int, default=8 << 20, help="Iterator readahead bytes")
    ap.add_argument("--pin", action="store_true", help="Pin data blocks during iteration")
    ap.add_argument("--async-io", action="store_true", help="Use async IO for iterator")
    ap.add_argument(
        "--output",
        default=None,
        help="Write results to this file instead of stdout. Helpful for larger DBs.",
    )
    args = ap.parse_args()

    # --- Open DB / CF ---
    # Prefer raw_mode for performance with binary keys; disable if your DB stores other encodings.
    opts = Options() if args.no_raw else Options(raw_mode=True)
    db = Rdict(args.db, options=opts)
    cf = db if args.cf is None else db.get_column_family(args.cf)

    # --- Infer k if absent ---
    # We peek at the first key to derive length. Empty DB is a hard error.
    if args.k is None:
        it = iter(cf.keys())
        try:
            first_key = next(it)
        except StopIteration:
            print("DB empty; cannot infer k", file=sys.stderr)
            sys.exit(2)
        k = len(first_key)
    else:
        k = args.k

    # --- Build constraints / fast check arrays ---
    alphabet = args.alphabet.encode("ascii").upper()
    rules = parse_rules(args, k, alphabet)
    allowed_by_pos = build_allowed_arrays(k, rules)

    # --- Plan enumeration size vs scan tightness ---
    P = choose_prefix_len(rules, k, args.max_enum, alphabet)
    prefixes = list(enumerated_prefixes(P, rules, alphabet))

    # --- Configure read options for tight, prefix-bounded scans ---
    ro = ReadOptions()
    # Readahead (useful for sequential scans over a bounded key range)
    _try_call(ro, ["set_readahead_size", "readahead_size"], args.readahead)
    # Restrict iteration within the key prefix if a prefix extractor is configured
    _try_call(ro, ["set_prefix_same_as_start", "prefix_same_as_start"], True)
    # Prefer prefix/index-guided seeks over total-order scans
    _try_call(ro, ["set_total_order_seek", "total_order_seek"], False)
    # Optional speed knobs
    if args.async_io:
        _try_call(ro, ["set_async_io", "async_io"], True)
    if args.pin:
        _try_call(ro, ["set_pin_data", "pin_data"], True)
    _try_call(ro, ["set_verify_checksums", "verify_checksums"], False)
    _try_call(ro, ["set_ignore_range_deletions", "ignore_range_deletions"], True)
    _try_call(ro, ["set_tailing", "tailing"], False)
    _try_call(ro, ["set_fill_cache", "fill_cache"], False)

    # --- Output stream (file or stdout) ---
    out = sys.stdout
    _fh = None
    if args.output:
        _fh = open(args.output, "w")
        out = _fh

    # --- Core loop: range-restricted iteration per enumerated prefix ---
    printed = 0
    for pref in prefixes:
        lo, hi = seek_range_for_prefix(pref, k)
        _set_bound(ro, "lower", lo)
        _set_bound(ro, "upper", hi)

        it = cf.iter(ro)
        it.seek(lo)

        while it.valid():
            key = it.key()

            # Exit this prefix range early once we step out of the (lo, hi) window.
            # - len(key) != k: skip malformed entries (defensive check).
            # - not key.startswith(pref): we moved past this prefix.
            if len(key) != k or not key.startswith(pref):
                break

            # Apply per-position constraints cheaply in-memory.
            if match_kmer_bytes(key, allowed_by_pos):
                val = it.value()
                cnt: Optional[int] = None

                # Decode value per chosen format (fall back to ASCII int)
                if args.values_are_u32_be and len(val) == 4:
                    cnt = int.from_bytes(val, "big", signed=False)
                elif args.values_are_u64 and len(val) == 8:
                    cnt = int.from_bytes(val, "little", signed=False)
                else:
                    try:
                        cnt = int(val)
                    except Exception:
                        # Leave as None => treat as "unknown count" and print if not filtering it out.
                        cnt = None

                # If we couldn't decode a count, we still print the k-mer unless min-count excludes it.
                if cnt is None or cnt >= args.min_count:
                    # Keys may be ASCII (preferred) or arbitrary bytes.
                    try:
                        s = key.decode("ascii")
                    except UnicodeDecodeError:
                        s = key.hex()

                    if cnt is None:
                        print(s, file=out)
                    else:
                        print(f"{s}\t{cnt}", file=out)

                    printed += 1
                    if args.limit and printed >= args.limit:
                        db.close()
                        if _fh:
                            _fh.close()
                        return

            it.next()

    # --- Cleanup ---
    db.close()
    if _fh:
        _fh.close()


if __name__ == "__main__":
    main()
