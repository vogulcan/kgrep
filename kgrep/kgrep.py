#!/usr/bin/env python3
"""
kgrep.py — Fast positional k-mer grep over RocksDict/RocksDB.

This script scans a RocksDict/RocksDB where:
- Keys are fixed-length k-mers (bytes/ASCII).
- Values are counts (auto-inferred as ASCII int, little-endian u64, or big-endian u32),
  unless explicitly overridden via flags.

It applies position-specific rules (via --mask or --rule) and performs
bounded, prefix-guided scans to minimize IO, using iterate_lower_bound /
iterate_upper_bound and prefix_same_as_start when available.

Default alphabet : "ACGT". You can override with --alphabet.

Typical uses
------------
- Mask syntax with count filter (auto-infer values):
    kgrep --db /dbpath/dbfolder --mask "A[CT]?AA" --min-count 10

- Rule syntax with output file:
    kgrep --db /dbpath/dbfolder --rule "1:A;2:CT" \
               --min-count 1 --output output.txt

- Basic k-mer search without constraints (lists all):
    kgrep --db /dbpath/dbfolder --k 10

- Complex rule with multiple positions:
    kgrep --db /dbpath/dbfolder --rule "3:A;5:CG;10:T" --limit 1000

- **NEW: OR logic examples**:
    # Match k-mers with (7:A AND 8:A) OR (17:T AND 18:T)
    kgrep --db /dbpath/dbfolder --rule "7:A;8:A|17:T;18:T"

    # Same as above, using multiple --rule flags:
    kgrep --db /dbpath/dbfolder --rule "7:A;8:A" --rule "17:T;18:T"

Mask syntax table
=================
Symbol     Meaning
---------- -----------------------------------------------------------------
?          Any character from the alphabet (wildcard)
[AC]       Allowed set at this position (only listed characters are allowed)
G          Exact match for a single required character (here: 'G')
"""


from __future__ import annotations

import argparse
import itertools
import sys
from typing import List, Dict, Iterable, Tuple, Optional, Set

from rocksdict import Rdict, Options, ReadOptions


# ---------- Portable helpers over rocksdict API differences ----------
def _try_call(obj, candidates: Iterable[str], *args) -> bool:
    for name in candidates:
        fn = getattr(obj, name, None)
        if callable(fn):
            fn(*args)
            return True
    return False


def _set_bound(ro: ReadOptions, kind: str, key: bytes) -> None:
    if kind == "lower":
        if _try_call(ro, ["set_iterate_lower_bound", "iterate_lower_bound"], key):
            return
    else:
        if _try_call(ro, ["set_iterate_upper_bound", "iterate_upper_bound"], key):
            return
    # If neither bound setter exists, we proceed without bounds (correct but slower).


# ---------- Rule parsing (now supports OR groups) ----------
def _parse_rule_group(group: str) -> Dict[int, bytes]:
    """
    Parse a single AND-group like '3:A;5:CG;10:T' into {pos: allowed_bytes}.
    """
    rules: Dict[int, bytes] = {}
    group = group.strip()
    if not group:
        return rules
    for term in group.split(";"):
        term = term.strip()
        if not term:
            continue
        p_str, vals = term.split(":", 1)
        p = int(p_str)
        allowed = bytes(sorted(set(vals.upper().encode("ascii"))))
        rules[p] = allowed
    return rules


def parse_rules_or(
    args: argparse.Namespace, k: int, alphabet: bytes
) -> List[Dict[int, bytes]]:
    """
    Return a list of rule dicts (OR groups). Each dict is an AND of its positions.

    Priority:
        --mask (single group, AND semantics) overrides --rule if provided.
        Otherwise, --rule may be given once with '|' separators or multiple times.
    """
    # Mask path (unchanged semantics: single AND group)
    if args.mask:
        pos = 1
        i = 0
        s = args.mask
        rules: Dict[int, bytes] = {}
        while i < len(s):
            c = s[i]
            if c == "?":
                pos += 1
                i += 1
                continue
            if c == "[":
                j = s.find("]", i + 1)
                if j == -1:
                    raise ValueError("Unclosed '[' in --mask")
                allowed = s[i + 1 : j].encode("ascii").upper()
                rules[pos] = bytes(sorted(set(allowed)))
                pos += 1
                i = j + 1
                continue
            rules[pos] = bytes([ord(c.upper())])
            pos += 1
            i += 1
        return [rules]  # single AND group

    # Rule path with OR support
    groups: List[Dict[int, bytes]] = []
    if args.rule:
        # args.rule may be a list (when repeated) or a single string
        rule_args = args.rule if isinstance(args.rule, list) else [args.rule]
        for ra in rule_args:
            # Split by '|' to allow inline OR
            for grp in ra.split("|"):
                parsed = _parse_rule_group(grp)
                if parsed:
                    groups.append(parsed)
    return groups


# ---------- Prefix planning ----------
def choose_prefix_len(
    rules: Dict[int, bytes],
    k: int,
    max_enum: int,
    alphabet: bytes,
) -> int:
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
            return P - 1
        if pos >= last_constrained:
            return P
    return P


def enumerated_prefixes(
    P: int, rules: Dict[int, bytes], alphabet: bytes
) -> Iterable[bytes]:
    if P <= 0:
        yield b""
        return
    slots = [rules.get(pos, alphabet) for pos in range(1, P + 1)]
    for tup in itertools.product(*slots):
        yield bytes(tup)


# ---------- Fast checking ----------
def build_allowed_arrays(k: int, rules: Dict[int, bytes]) -> List[Optional[bytes]]:
    arr: List[Optional[bytes]] = [None] * (k + 1)
    for pos, allowed in rules.items():
        if 1 <= pos <= k:
            arr[pos] = allowed
    return arr


def match_kmer_bytes_single(kmer: bytes, allowed_by_pos: List[Optional[bytes]]) -> bool:
    for i, allowed in enumerate(allowed_by_pos[1:], start=0):
        if allowed is None:
            continue
        if kmer[i] not in allowed:
            return False
    return True


def match_kmer_bytes_any(kmer: bytes, groups: List[List[Optional[bytes]]]) -> bool:
    """
    OR across groups: match if kmer matches ANY group's per-position constraints.
    If no groups provided, match everything (back-compat with no rules).
    """
    if not groups:
        return True
    for arr in groups:
        if match_kmer_bytes_single(kmer, arr):
            return True
    return False


# ---------- Bounds / range calculation ----------
def upper_bound_for_prefix(prefix: bytes) -> bytes:
    p = bytearray(prefix)
    for i in reversed(range(len(p))):
        if p[i] != 0xFF:
            p[i] += 1
            return bytes(p[: i + 1])
    return prefix + b"\x00"


def seek_range_for_prefix(prefix: bytes, k: int) -> Tuple[bytes, bytes]:
    if len(prefix) == 0:
        lo = b"\x00" * k
        hi = upper_bound_for_prefix(b"\xff" * k)
        return lo, hi
    lo = prefix + b"\x00" * (k - len(prefix))
    hi = upper_bound_for_prefix(prefix)
    return lo, hi


# ---------- Value decoding (auto-infer with optional overrides) ----------
def _looks_like_ascii_int(b: bytes) -> Optional[int]:
    try:
        s = b.decode("ascii", errors="strict").strip()
    except UnicodeDecodeError:
        return None
    if not s:
        return None
    core = s[1:] if s[0] in "+-" else s
    if not core or any(ch not in "0123456789" for ch in core):
        return None
    try:
        return int(s)
    except Exception:
        return None


def decode_count(val: bytes, args: argparse.Namespace) -> Optional[int]:
    if args.values_are_u64:
        if len(val) != 8:
            ascii_int = _looks_like_ascii_int(val)
            return ascii_int
        return int.from_bytes(val, "little", signed=False)

    if args.values_are_u32_be:
        if len(val) != 4:
            ascii_int = _looks_like_ascii_int(val)
            return ascii_int
        return int.from_bytes(val, "big", signed=False)

    if len(val) == 8:
        return int.from_bytes(val, "little", signed=False)
    ascii_int = _looks_like_ascii_int(val)
    if ascii_int is not None:
        return ascii_int
    if len(val) == 4:
        return int.from_bytes(val, "big", signed=False)
    return None


# ---------- Main entrypoint ----------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="k-mer rule grep over RocksDict (auto-infers value encoding)"
    )
    ap.add_argument("--db", required=True, help="Path to RocksDict / RocksDB")
    ap.add_argument(
        "--cf", default=None, help="Column family name (default CF if omitted)"
    )
    ap.add_argument(
        "--k", type=int, help="k-mer length in bytes (auto from first key if omitted)"
    )
    # NOTE: --rule now supports OR via multiple flags or '|' separators.
    ap.add_argument(
        "--rule",
        action="append",
        help='Rules like "3:A;5:CG;10:T". Use multiple --rule flags or "|" to OR groups.',
        default=None,
    )
    ap.add_argument(
        "--mask",
        help='Compact mask like "A[CT]?G??T" (?=any, [..]=set). Overrides --rule (single AND).',
        default=None,
    )
    ap.add_argument("--alphabet", default="ACGT", help="Alphabet, default ACGT")
    ap.add_argument(
        "--min-count", type=int, default=0, help="Only output k-mers with count >= this"
    )
    ap.add_argument(
        "--max-enum",
        type=int,
        default=200000,
        help="Max prefix combinations to enumerate",
    )
    ap.add_argument(
        "--limit", type=int, default=0, help="Stop after N matches (0 = unlimited)"
    )
    ap.add_argument(
        "--values-are-u32-be",
        action="store_true",
        help="Override: interpret value as big-endian u32 count",
    )
    ap.add_argument(
        "--values-are-u64",
        action="store_true",
        help="Override: interpret value as little-endian u64 count",
    )
    ap.add_argument(
        "--no-raw",
        action="store_true",
        help="Open without raw_mode (use if your DB isn’t raw bytes)",
    )
    ap.add_argument(
        "--readahead", type=int, default=8 << 20, help="Iterator readahead bytes"
    )
    ap.add_argument(
        "--pin", action="store_true", help="Pin data blocks during iteration"
    )
    ap.add_argument("--async-io", action="store_true", help="Use async IO for iterator")
    ap.add_argument(
        "--output",
        default=None,
        help="Write results to this file instead of stdout. Helpful for larger DBs.",
    )
    args = ap.parse_args()

    # --- Open DB / CF ---
    opts = Options() if args.no_raw else Options(raw_mode=True)
    db = Rdict(args.db, options=opts)
    cf = db if args.cf is None else db.get_column_family(args.cf)

    # --- Infer k if absent ---
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

    # --- Build constraints (now potentially multiple OR groups) ---
    alphabet = args.alphabet.encode("ascii").upper()
    rule_groups: List[Dict[int, bytes]] = parse_rules_or(args, k, alphabet)
    allowed_groups: List[List[Optional[bytes]]] = [
        build_allowed_arrays(k, g) for g in rule_groups
    ]

    # --- Plan enumeration: union of prefixes across groups (dedup) ---
    prefixes_set: Set[bytes] = set()
    if not rule_groups:
        # No constraints: one empty prefix to cover full keyspace
        prefixes = [b""]
    else:
        for g in rule_groups:
            P_g = choose_prefix_len(g, k, args.max_enum, alphabet)
            for pref in enumerated_prefixes(P_g, g, alphabet):
                prefixes_set.add(pref)
        prefixes = sorted(prefixes_set) if prefixes_set else [b""]

    # --- Configure read options for tight, prefix-bounded scans ---
    ro = ReadOptions()
    _try_call(ro, ["set_readahead_size", "readahead_size"], args.readahead)
    _try_call(ro, ["set_prefix_same_as_start", "prefix_same_as_start"], True)
    _try_call(ro, ["set_total_order_seek", "total_order_seek"], False)
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
            # Exit this prefix range early once we step out of (lo, hi) or the prefix
            if len(key) != k or not key.startswith(pref):
                break

            # Apply OR-of-ANDs positional constraints.
            if match_kmer_bytes_any(key, allowed_groups):
                val = it.value()
                cnt = decode_count(val, args)

                if cnt is None or cnt >= args.min_count:
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
