#!/usr/bin/env python3

from __future__ import annotations

import argparse
import io
import os
import re
import sys
from collections import OrderedDict, deque
from typing import Iterable, Iterator, List, Dict, Tuple, Optional, Set

from rocksdict import Rdict, Options, ReadOptions, AccessType


# ---------- rocksdict helpers ----------
def _try_call(obj, candidates: Iterable[str], *args) -> bool:
    for name in candidates:
        fn = getattr(obj, name, None)
        if callable(fn):
            try:
                fn(*args)
                return True
            except TypeError:
                continue
    return False


def _set_bound(ro: ReadOptions, kind: str, key: bytes) -> None:
    if kind == "lower":
        if _try_call(ro, ["set_iterate_lower_bound", "iterate_lower_bound"], key):
            return
    else:
        if _try_call(ro, ["set_iterate_upper_bound", "iterate_upper_bound"], key):
            return


# ---------- rule parsing (scan mode) ----------
_TERM_RE = re.compile(r"^\s*\d+\s*:\s*[!-~]+\s*$")


def _parse_rule_group(group: str) -> Dict[int, bytes]:
    rules: Dict[int, bytes] = {}
    group = group.strip()
    if not group:
        return rules
    for term in group.split(";"):
        term = term.strip()
        if not term or not _TERM_RE.match(term):
            continue
        p_str, vals = term.split(":", 1)
        try:
            p = int(p_str)
        except ValueError:
            continue
        allowed = bytes(sorted(set(vals.upper().encode("ascii", "ignore"))))
        if allowed:
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


# ---------- alphabet inference (scan mode) ----------
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


# ---------- scan helpers ----------
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
    import itertools as _it
    for tup in _it.product(*slots):
        yield bytes(tup)


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


# ---------- bounds ----------
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


# ---------- value decoding ----------
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


# ---------- chunk input streaming ----------
def _iter_chunks_from_path(path: str) -> Iterator[str]:
    with open(path, "r", buffering=1024 * 1024) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            yield s


def _iter_chunks_from_stdin() -> Iterator[str]:
    data_stream = io.TextIOWrapper(sys.stdin.buffer, encoding="utf-8", errors="ignore", line_buffering=False)
    for line in data_stream:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        yield s


def _iter_ascii_chunks(args: argparse.Namespace) -> Iterator[bytes]:
    # inline --chunks first (if any), then file/stdin
    if args.chunks:
        for item in args.chunks:
            parts = (p.strip() for p in item.split(","))
            for s in parts:
                if s:
                    try:
                        yield s.encode("ascii", "strict")
                    except UnicodeEncodeError:
                        sys.stderr.write(f"[kgrep] warning: non-ASCII kmer '{s}' – skipping\n")
    if args.chunks_file:
        if args.chunks_file == "-":
            for s in _iter_chunks_from_stdin():
                try:
                    yield s.encode("ascii", "strict")
                except UnicodeEncodeError:
                    sys.stderr.write(f"[kgrep] warning: non-ASCII kmer '{s}' – skipping\n")
        else:
            for s in _iter_chunks_from_path(args.chunks_file):
                try:
                    yield s.encode("ascii", "strict")
                except UnicodeEncodeError:
                    sys.stderr.write(f"[kgrep] warning: non-ASCII kmer '{s}' – skipping\n")


def _iter_hex_chunks(args: argparse.Namespace) -> Iterator[bytes]:
    if args.chunks:
        for item in args.chunks:
            parts = (p.strip() for p in item.split(","))
            for s in parts:
                if not s:
                    continue
                s2 = s[2:] if s.lower().startswith("0x") else s
                try:
                    yield bytes.fromhex(s2)
                except ValueError:
                    sys.stderr.write(f"[kgrep] warning: could not parse hex kmer '{s}' – skipping\n")
    if args.chunks_file:
        src = _iter_chunks_from_stdin() if args.chunks_file == "-" else _iter_chunks_from_path(args.chunks_file)
        for s in src:
            s2 = s[2:] if s.lower().startswith("0x") else s
            try:
                yield bytes.fromhex(s2)
            except ValueError:
                sys.stderr.write(f"[kgrep] warning: could not parse hex kmer '{s}' – skipping\n")


def iter_chunk_keys(args: argparse.Namespace) -> Iterator[bytes]:
    if args.chunks_hex:
        yield from _iter_hex_chunks(args)
    else:
        yield from _iter_ascii_chunks(args)


def count_lines_in_file(path: str) -> int:
    # Count non-empty, non-comment lines for progress %
    n = 0
    with open(path, "rb", buffering=1024 * 1024) as fh:
        for raw in fh:
            if not raw:
                continue
            # quick check: skip blank or '#'
            if raw[:1] == b"\n" or raw[:1] == b"#":
                # still need to skip lines that are just spaces
                # (rare) – fall back to strip for those small cases
                s = raw.strip()
                if not s or s.startswith(b"#"):
                    continue
            n += 1
    return n


# ---------- reverse complement ----------
_RC_TABLE = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


def reverse_complement_bytes(seq: bytes) -> Optional[bytes]:
    if not seq:
        return b""
    for bch in seq:
        if bch not in b"ACGTacgt":
            return None
    return seq.translate(_RC_TABLE)[::-1]


# ---------- robust point get with positive-only cache ----------
def _is_byteslike(x) -> bool:
    return isinstance(x, (bytes, bytearray, memoryview))


class Prober:
    def __init__(self, cf, ro: ReadOptions, allow_rc: bool, chunks_hex: bool, k: int):
        self.cf = cf
        self.ro = ro
        self.rc_enabled = bool(allow_rc and not chunks_hex)
        self.k = k
        self._get_fn = getattr(cf, "get", None)
        self._cache: Dict[bytes, Tuple[bool, Optional[bytes], Optional[bytes]]] = {}

        if allow_rc and chunks_hex:
            sys.stderr.write("[kgrep] --rc ignored for hex inputs; cannot compute reverse complement.\n")

    def _point_get(self, key: bytes) -> Optional[bytes]:
        fn = self._get_fn
        if callable(fn):
            for args in ((key,), (key, self.ro), (self.ro, key)):
                try:
                    val = fn(*args)
                except TypeError:
                    continue
                if _is_byteslike(val) or val is None:
                    if isinstance(val, bytearray):
                        return bytes(val)
                    if isinstance(val, memoryview):
                        return val.tobytes()
                    return val
        it = self.cf.iter(self.ro)
        it.seek(key)
        if it.valid() and it.key() == key:
            return it.value()
        return None

    def probe(self, kmer: bytes) -> Tuple[bool, Optional[bytes], Optional[bytes]]:
        cached = self._cache.get(kmer)
        if cached is not None:
            return cached

        val = self._point_get(kmer)
        if val is not None:
            res = (True, kmer, val)
            self._cache[kmer] = res
            return res

        if self.rc_enabled:
            rc = reverse_complement_bytes(kmer)
            if rc is not None and len(rc) == self.k:
                cached_rc = self._cache.get(rc)
                if cached_rc is not None:
                    self._cache[kmer] = cached_rc
                    return cached_rc
                val_rc = self._point_get(rc)
                if val_rc is not None:
                    res = (True, rc, val_rc)
                    self._cache[kmer] = res
                    self._cache[rc] = res
                    return res

        return (False, None, None)  # no negative caching


# ---------- bounded output de-dup ----------
class LRUSet:
    def __init__(self, capacity: int):
        self.cap = max(1, capacity)
        self.od = OrderedDict()

    def add(self, key: bytes) -> bool:
        if key in self.od:
            self.od.move_to_end(key)
            return False
        self.od[key] = None
        if len(self.od) > self.cap:
            self.od.popitem(last=False)
        return True


# ---------- main ----------
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
    ap.add_argument("--readahead", type=int, default=0,  # 0 is best for point lookups
                    help="Readahead bytes for iterators (0 is typical for point lookups).")
    ap.add_argument("--pin", action="store_true")
    ap.add_argument("--async-io", action="store_true")
    ap.add_argument("--output", default=None)
    ap.add_argument("--read-only", action="store_true",
                    help="Open DB in read-only mode (allows parallel readers without taking the lock).")

    # chunk retrieval options
    ap.add_argument("--chunks", action="append", default=None,
                    help="ASCII kmers to retrieve; can be repeated or comma-separated.")
    ap.add_argument("--chunks-file", default=None,
                    help="Path to file with one kmer per line (use '-' for stdin).")
    ap.add_argument("--chunks-hex", action="store_true",
                    help="Interpret provided kmers as hex-encoded bytes (0x prefix optional).")
    ap.add_argument("--progress", action="store_true",
                    help="Show progress on stderr during chunk retrieval.")
    ap.add_argument("--progress-step", type=int, default=100,
                    help="Update progress every N percent when --progress is set (default 100).")
    ap.add_argument("--rc", action="store_true",
                    help="If exact kmer not found, try reverse-complement (ACGT only).")

    # performance / memory knobs
    ap.add_argument("--assume-unique", action="store_true",
                    help="Assume inputs are unique; do not allocate input dedup structures.")
    ap.add_argument("--no-output-dedup", action="store_true",
                    help="Do not de-duplicate outputs (print every hit).")
    ap.add_argument("--buffer-size", type=int, default=10_000,
                    help="Batch size for buffered stdout writes in chunk retrieval.")
    ap.add_argument("--output-lru", type=int, default=200_000,
                    help="LRU size for output de-dup when enabled.")

    args = ap.parse_args()

    # open DB
    opts = Options() if args.no_raw else Options(raw_mode=True)
    if args.read_only:
        db = Rdict(args.db, options=opts, access_type=AccessType.read_only())
    else:
        db = Rdict(args.db, options=opts)
    cf = db if args.cf is None else db.get_column_family(args.cf)

    # value decoder
    if args.values == "auto":
        decoder = ValueDecoder.infer_from_cf(cf)
        sys.stderr.write(f"[kgrep] values=auto → inferred {decoder.mode}\n")
    else:
        decoder = ValueDecoder(args.values)

    # output
    out = sys.stdout
    _fh = None
    if args.output:
        _fh = open(args.output, "w", buffering=1024 * 1024)
        out = _fh

    # ReadOptions for point lookups
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

    # --------------- CHUNK RETRIEVAL (streaming) ---------------
    if args.chunks or args.chunks_file:
        # Determine k:
        # prefer user-specified; otherwise peek the first kmer from the stream;
        # if nothing provided, fall back to DB's first key length.
        stream = iter_chunk_keys(args)

        first: Optional[bytes] = None
        try:
            first = next(stream)
        except StopIteration:
            # no inputs; fall back to DB-based k (for consistency with original behavior)
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
            # nothing to do
            db.close()
            if _fh:
                _fh.close()
            return

        # compute k
        if args.k is not None:
            k = args.k
        else:
            k = len(first)

        # progress: line count only if file path (not stdin)
        total_inputs: Optional[int] = None
        if args.progress:
            if args.chunks_file and args.chunks_file != "-":
                try:
                    total_inputs = count_lines_in_file(args.chunks_file)
                except Exception:
                    total_inputs = None  # fallback to count-only progress
            else:
                total_inputs = None

        prober = Prober(cf, ro, allow_rc=args.rc, chunks_hex=args.chunks_hex, k=k)

        # output de-dup (optional)
        use_output_dedup = not args.no_output_dedup
        lru_outputs = LRUSet(capacity=max(1, int(args.output_lru))) if use_output_dedup else None

        # counters & buffers
        show_progress = bool(args.progress)
        step = args.progress_step if 1 <= args.progress_step <= 100 else 100
        next_threshold = step
        seen = 0
        found = 0
        printed = 0
        bufsize = max(1, int(args.buffer_size))
        outbuf = deque()

        # process the first item then the rest
        def process_one(kmer: bytes) -> None:
            nonlocal found, printed
            if len(kmer) != k:
                # strict length check (keep behavior predictable)
                return
            ok, db_key, db_val = prober.probe(kmer)
            if not ok:
                return
            found += 1
            # optional de-dup by DB key
            if lru_outputs is not None:
                if db_key is None or not lru_outputs.add(db_key):
                    return
            # decode + buffer print
            vb: Optional[bytes]
            if isinstance(db_val, bytearray):
                vb = bytes(db_val)
            elif isinstance(db_val, memoryview):
                vb = db_val.tobytes()
            else:
                vb = db_val
            cnt = None if decoder.mode == "none" else decoder.decode(vb)  # type: ignore[arg-type]
            try:
                s_key = db_key.decode("ascii") if db_key is not None else ""  # type: ignore[union-attr]
            except UnicodeDecodeError:
                s_key = db_key.hex() if db_key is not None else ""  # type: ignore[union-attr]
            if cnt is None:
                outbuf.append(f"{s_key}\n")
            else:
                if cnt >= args.min_count:
                    outbuf.append(f"{s_key}\t{cnt}\n")
                else:
                    return
            if len(outbuf) >= bufsize:
                out.writelines(outbuf)
                outbuf.clear()
                if _fh:
                    out.flush()
            printed += 1

        # first element
        seen += 1
        process_one(first)
        # rest of stream
        for kmer in stream:
            seen += 1
            process_one(kmer)

            if args.limit and printed >= args.limit:
                break

            if show_progress:
                if total_inputs:
                    pct = int(seen * 100 / total_inputs)
                    if pct >= next_threshold or seen == total_inputs:
                        sys.stderr.write(f"\r[kgrep] chunks: {seen}/{total_inputs} ({pct}%)")
                        sys.stderr.flush()
                        while next_threshold <= pct:
                            next_threshold += step
                else:
                    # count-only progress (no %)
                    if seen % (100_000) == 0:
                        sys.stderr.write(f"\r[kgrep] chunks: {seen}")
                        sys.stderr.flush()

        # final flush
        if outbuf:
            out.writelines(outbuf)
            outbuf.clear()
            if _fh:
                out.flush()
        if show_progress:
            sys.stderr.write("\n")
            sys.stderr.flush()

        # not-found report (per-line, since we assume unique)
        not_found = seen - found
        if total_inputs is not None and total_inputs != seen:
            # in case of early --limit or parse skips
            total_display = seen
        else:
            total_display = seen
        sys.stderr.write(f"[kgrep] not found: {not_found}/{total_display}\n")
        sys.stderr.flush()

        db.close()
        if _fh:
            _fh.close()
        return

    # --------------- SCAN / ENUMERATION MODE (unchanged logic) ---------------
    # determine k (for scan mode)
    if args.k is not None:
        k = args.k
    else:
        it = iter(cf.keys())
        try:
            first_key = next(it)
        except StopIteration:
            print("DB empty", file=sys.stderr)
            sys.exit(2)
        k = len(first_key)

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

    out = sys.stdout
    _fh = None
    if args.output:
        _fh = open(args.output, "w", buffering=1024 * 1024)
        out = _fh

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

    printed = 0
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
            if match_kmer_bytes_any(key, allowed_groups, user_provided_rules=user_provided_rules):
                val = it.value()
                cnt = None if decoder.mode == "none" else decoder.decode(val)
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
