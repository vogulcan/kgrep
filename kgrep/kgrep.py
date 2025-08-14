#!/usr/bin/env python3

from __future__ import annotations

import argparse
import itertools
import re
import sys
from typing import List, Dict, Iterable, Tuple, Optional, Set

from rocksdict import Rdict, Options, ReadOptions


# ---------- Helpers for rocksdict API differences ----------
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


# ---------- Rule parsing (supports OR groups) ----------
_TERM_RE = re.compile(r"^\s*\d+\s*:\s*[!-~]+\s*$")


def _parse_rule_group(group: str) -> Dict[int, bytes]:
    rules: Dict[int, bytes] = {}
    group = group.strip()
    if not group:
        return rules
    for term in group.split(";"):
        term = term.strip()
        if not term:
            continue
        if not _TERM_RE.match(term):
            continue
        p_str, vals = term.split(":", 1)
        try:
            p = int(p_str)
        except ValueError:
            continue
        allowed = bytes(sorted(set(vals.upper().encode("ascii", "ignore"))))
        if not allowed:
            continue
        rules[p] = allowed
    return rules


def parse_rules_or(args: argparse.Namespace) -> List[Dict[int, bytes]]:
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
                allowed = s[i + 1:j].encode("ascii", "ignore").upper()
                allowed = bytes(sorted(set(allowed)))
                if allowed:
                    rules[pos] = allowed
                pos += 1
                i = j + 1
                continue
            rules[pos] = bytes([ord(c.upper())])
            pos += 1
            i += 1
        return [rules]

    groups: List[Dict[int, bytes]] = []
    if args.rule:
        rule_args = args.rule if isinstance(args.rule, list) else [args.rule]
        for ra in rule_args:
            for grp in str(ra).split("|"):
                parsed = _parse_rule_group(grp)
                if parsed:
                    groups.append(parsed)
    return groups


# ---------- Alphabet inference (optional) ----------
PRINTABLE_UPPER = set(range(65, 91))
PRINTABLE = set(range(32, 127))


def guess_alphabet(cf, k: int, sample: int = 512) -> bytes:
    ro = ReadOptions()
    _try_call(ro, ["set_total_order_seek", "total_order_seek"], True)
    _try_call(ro, ["set_fill_cache", "fill_cache"], False)
    it = cf.iter(ro)
    it.seek_to_first()

    seen: Set[int] = set()
    n = 0
    while it.valid() and n < sample:
        key = it.key()
        if len(key) == k:
            for b in key:
                seen.add(b)
        n += 1
        it.next()

    if not seen:
        return b"ACGT"
    if seen.issubset(PRINTABLE_UPPER) and len(seen) <= 32:
        return bytes(sorted(seen))
    if seen.issubset(PRINTABLE) and len(seen) <= 32:
        return bytes(sorted(seen))
    return b"ACGT"


# ---------- Prefix planning ----------
def choose_prefix_len(rules: Dict[int, bytes], k: int, max_enum: int, alphabet: bytes) -> int:
    if not rules:
        return 0
    last_constrained = max((p for p in rules.keys() if 1 <= p <= k), default=0)
    P = 0
    combos = 1
    for pos in range(1, k + 1):
        P = pos
        card = len(rules.get(pos, alphabet))
        combos *= card
        if combos > max_enum:
            return P - 1
        if pos >= last_constrained:
            return P
    return P


def enumerated_prefixes(P: int, rules: Dict[int, bytes], alphabet: bytes) -> Iterable[bytes]:
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


def match_kmer_bytes_any(kmer: bytes, groups: List[List[Optional[bytes]]], *, user_provided_rules: bool) -> bool:
    if not groups:
        return (not user_provided_rules)
    for arr in groups:
        if match_kmer_bytes_single(kmer, arr):
            return True
    return False


# ---------- Bounds ----------
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


# ---------- Value decoding ----------
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


class ValueDecoder:
    def __init__(self, mode: str) -> None:
        self.mode = mode

    def decode(self, val: bytes) -> Optional[int]:
        m = self.mode
        if m == "none":
            return None
        if m == "ascii":
            return _looks_like_ascii_int(val)
        if m == "u64le":
            if len(val) != 8:
                return _looks_like_ascii_int(val)
            return int.from_bytes(val, "little", signed=False)
        if m == "u64be":
            if len(val) != 8:
                return _looks_like_ascii_int(val)
            return int.from_bytes(val, "big", signed=False)
        if m == "u32le":
            if len(val) != 4:
                return _looks_like_ascii_int(val)
            return int.from_bytes(val, "little", signed=False)
        if m == "u32be":
            if len(val) != 4:
                return _looks_like_ascii_int(val)
            return int.from_bytes(val, "big", signed=False)
        raise RuntimeError(f"Unknown mode {m}")

    @staticmethod
    def infer_from_cf(cf, sample: int = 512) -> "ValueDecoder":
        ro = ReadOptions()
        _try_call(ro, ["set_total_order_seek", "total_order_seek"], True)
        _try_call(ro, ["set_fill_cache", "fill_cache"], False)
        it = cf.iter(ro)
        it.seek_to_first()

        vals: List[bytes] = []
        while it.valid() and len(vals) < sample:
            vals.append(it.value())
            it.next()

        if not vals:
            return ValueDecoder("u64le")

        ascii_hits = sum(1 for v in vals if _looks_like_ascii_int(v) is not None)
        if ascii_hits / len(vals) >= 0.9:
            return ValueDecoder("ascii")

        def score_numeric(length: int, endian: str) -> Tuple[float, str]:
            mode = f"u{length*8}{'le' if endian=='little' else 'be'}"
            subset = [v for v in vals if len(v) == length]
            if not subset:
                return (float("inf"), mode)
            decoded = [int.from_bytes(v, endian, signed=False) for v in subset]
            huge = sum(1 for x in decoded if x > 10**15)
            median = sorted(decoded)[len(decoded) // 2]
            return (huge * 1_000_000 + median, mode)

        scores = []
        for L in (8, 4):
            for e in ("little", "big"):
                scores.append(score_numeric(L, e))

        best_score, best_mode = min(scores, key=lambda t: t[0])
        if best_score == float("inf"):
            return ValueDecoder("u64le")
        return ValueDecoder(best_mode)


# ---------- Main ----------
def main() -> None:
    ap = argparse.ArgumentParser(description="k-mer grep over RocksDict")
    ap.add_argument("--db", required=True)
    ap.add_argument("--cf", default=None)
    ap.add_argument("--k", type=int)
    ap.add_argument("--rule", action="append", default=None)
    ap.add_argument("--mask", default=None)
    ap.add_argument("--alphabet", default="ACGT", help="Use 'auto' to infer from keys")
    ap.add_argument("--min-count", type=int, default=0)
    ap.add_argument("--max-enum", type=int, default=200000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--values", default="auto",
                    choices=["auto", "ascii", "u32le", "u32be", "u64le", "u64be", "none"])
    ap.add_argument("--no-raw", action="store_true")
    ap.add_argument("--readahead", type=int, default=8 << 20)
    ap.add_argument("--pin", action="store_true")
    ap.add_argument("--async-io", action="store_true")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    opts = Options() if args.no_raw else Options(raw_mode=True)
    db = Rdict(args.db, options=opts)
    cf = db if args.cf is None else db.get_column_family(args.cf)

    if args.k is None:
        it = iter(cf.keys())
        try:
            first_key = next(it)
        except StopIteration:
            print("DB empty", file=sys.stderr)
            sys.exit(2)
        k = len(first_key)
    else:
        k = args.k

    if str(args.alphabet).lower() == "auto":
        alphabet = guess_alphabet(cf, k)
        sys.stderr.write(f"[kgrep] alphabet=auto → inferred {alphabet.decode('ascii', 'ignore')}\n")
    else:
        alphabet = args.alphabet.encode("ascii", "ignore").upper()
        if not alphabet:
            alphabet = b"ACGT"

    rule_groups = parse_rules_or(args)
    sanitized_groups: List[Dict[int, bytes]] = []
    for g in rule_groups:
        sg = {p: allowed for p, allowed in g.items() if 1 <= p <= k and allowed}
        if sg:
            sanitized_groups.append(sg)
    rule_groups = sanitized_groups

    allowed_groups = [build_allowed_arrays(k, g) for g in rule_groups]
    user_provided_rules = bool(args.mask or args.rule)

    prefixes_set: Set[bytes] = set()
    if not rule_groups:
        prefixes = [b""]
    else:
        for g in rule_groups:
            P_g = choose_prefix_len(g, k, args.max_enum, alphabet)
            for pref in enumerated_prefixes(P_g, g, alphabet):
                prefixes_set.add(pref)
        prefixes = sorted(prefixes_set) if prefixes_set else [b""]

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

    out = sys.stdout
    _fh = None
    if args.output:
        _fh = open(args.output, "w")
        out = _fh

    if args.values == "auto":
        decoder = ValueDecoder.infer_from_cf(cf)
        sys.stderr.write(f"[kgrep] values=auto → inferred {decoder.mode}\n")
    else:
        decoder = ValueDecoder(args.values)

    printed = 0
    seen_keys: Set[bytes] = set()  # deduplication

    for pref in prefixes:
        lo, hi = seek_range_for_prefix(pref, k)
        _set_bound(ro, "lower", lo)
        _set_bound(ro, "upper", hi)

        it = cf.iter(ro)
        it.seek(lo)

        while it.valid():
            key = it.key()
            if len(key) != k or not key.startswith(pref):
                break
            if key in seen_keys:
                it.next()
                continue
            if match_kmer_bytes_any(key, allowed_groups, user_provided_rules=user_provided_rules):
                val = it.value()
                cnt = None if decoder.mode == "none" else decoder.decode(val)
                if cnt is None or cnt >= args.min_count:
                    seen_keys.add(key)
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

    db.close()
    if _fh:
        _fh.close()


if __name__ == "__main__":
    main()
