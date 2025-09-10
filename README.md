# kgrep

**Fast positional k‑mer grep over RocksDict / RocksDB**

`kgrep` scans and retrieves fixed‑length k‑mers from a RocksDict/RocksDB database. It supports position‑specific filters, exact lookups for large k‑mer lists, and high‑throughput streaming designed to scale linearly to tens of millions of queries.

---

## Features

- **Mask syntax** with wildcards and character sets (see *Mask Syntax*).
- **Rule syntax** with 1‑based positional constraints; supports **OR** logic across rule groups.
- Tight prefix‑bounded scans using RocksDB lower/upper iterate bounds.
- Automatic alphabet inference (default DNA `ACGT`), configurable via `--alphabet`.
- Optional output to file via `--output`.
- **Chunk retrieval** (exact lookups for N kmers) with:
  - Streaming, O(1) memory per line (no preloading of large files).
  - Optional reverse‑complement lookup via `--rc` (ACGT ASCII input).
  - Progress reporting to `stderr`.
  - Value decoding (`--values`) and count filtering (`--min-count`).
  - Buffered output for high throughput.
  - Optional output de‑duplication to avoid double prints when different inputs map to the same DB key (e.g., via RC).

---

## Installation

```bash
pip install "git+https://github.com/vogulcan/kgrep@master"
```

---

## Example Usage

### 1) Search with mask syntax and count filter

> **Note:** `--mask` overrides `--rule` and is applied from the beginning of the k‑mer.

```bash
kgrep --db /path/dbfolder --mask "A[CT]?AA" --min-count 10
```

### 2) Search with rule syntax and output to file

> **Note:** `--rule` positions are **1‑based**.

```bash
kgrep --db /path/dbfolder --rule "1:A;2:CT" --min-count 1 --output output.txt
```

### 3) List all k‑mers without constraints

```bash
kgrep --db /path/dbfolder
```

### 4) Complex multi‑position rule

```bash
kgrep --db /path/dbfolder --rule "3:A;5:CG;10:T" --limit 1000
```

### 5) OR logic between rule groups

```bash
# (7:A AND 8:A) OR (17:T AND 18:T)
kgrep --db /path/dbfolder --rule "7:A;8:A|17:T;18:T"

# Same as above with multiple --rule flags:
kgrep --db /path/dbfolder --rule "7:A;8:A" --rule "17:T;18:T"
```

### 6) Chunk retrieval (exact lookups for N kmers)

Retrieve values for an explicit list of kmers (all must have the same length):

```bash
# Inline ASCII kmers (comma-separated or multiple flags)
kgrep --db /path/dbfolder --chunks AAAAC --chunks TTTTG,CCCCA

# From a file (one kmer per line; lines starting with '#' are ignored)
kgrep --db /path/dbfolder --chunks-file kmers.txt

# Reverse-complement fallback (ACGT ASCII inputs)
kgrep --db /path/dbfolder --chunks-file kmers.txt --rc

# From STDIN
cat kmers.txt | kgrep --db /path/dbfolder --chunks-file -

# Hex input (for non-ASCII keys); 0x prefix optional
kgrep --db /path/dbfolder --chunks-hex --chunks 0x4141414143 --chunks 7474745447
```

#### Notes for chunk retrieval

- If `--k` is not provided, *k* is inferred from the first provided k‑mer; if none are provided, it falls back to the length of the first key in the DB (original behavior).
- Values are decoded per `--values` (`auto`, `ascii`, `u32le`, `u32be`, `u64le`, `u64be`, `none`), and filtered with `--min-count` when numeric.
- Output for each found key is either:
  - `<key>\t<count>` when a numeric value is decoded, or
  - `<key>` if `--values none` is used or the value cannot be decoded.
- Non‑ASCII keys are printed in hex.
- Missing keys are skipped silently (no output line).
- `--limit` still applies to the number of printed lines.
- Optional progress reporting to `stderr`:
  ```bash
  kgrep --db /path/dbfolder --chunks-file kmers.txt --progress --progress-step 100
  ```
  When `--chunks-file` is a normal file, a percentage is shown; for stdin, a running count is shown.

---

## Tuning Flags for Chunk Retrieval (Scales to 10M–50M lines)

These flags complement existing ones to keep throughput **linear** on huge input files.

- `--read-only`  
  Open DB in read‑only mode so multiple jobs can run in parallel without DB locks.

- `--values <mode>`  
  Prefer a fixed mode (e.g., `u32be`) to skip the sampling phase of `--values auto`.

- `--rc`  
  If an exact kmer isn’t found, try its reverse‑complement (ACGT ASCII only). If found, prints the **DB key that exists** and its value. Not found reporting counts a k‑mer as missing **only if both the k‑mer and its RC are absent**.

- `--progress`, `--progress-step N`  
  Show progress on `stderr`. Larger `N` (e.g., `1000`) = fewer TTY writes.

- `--readahead 0`  
  Recommended for point lookups (default for chunk retrieval in the latest version).

- `--pin`  
  Enable RocksDB pinning to reduce copies on reads.

- `--buffer-size N` (default: `10000`)  
  Buffer N output lines before flushing, reducing per‑line I/O overhead.

- `--no-output-dedup`  
  By default, kgrep de‑dups **outputs** by DB key using a bounded LRU to avoid double prints when different inputs map to the same DB key (e.g., via `--rc`). Use this flag to **disable** de‑dup and print every hit.

- `--output-lru N` (default: `200000`)  
  Size of the bounded LRU used for output de‑duplication (ignored when `--no-output-dedup` is set).

- `--assume-unique`  
  Explicitly indicate inputs are unique. The streaming mode already assumes uniqueness and does **not** allocate large input de‑dup structures; this flag simply documents the assumption and keeps behavior explicit.

### Fast recipes

**Fast single-process, minimal overhead**
```bash
uv run kgrep --db /path/dbfolder \
  --read-only --rc --chunks-file kmers.txt \
  --values u32be --readahead 0 --pin \
  --progress --progress-step 1000 \
  --buffer-size 50000 \
  --output hits.txt
```

**Parallel sharding (8 jobs), safe with `--read-only`**
```bash
split -n l:8 --additional-suffix=.txt kmers.txt kmers.part.

parallel -j8 --linebuffer \
  'uv run kgrep --db /path/dbfolder \
     --read-only --rc --chunks-file {} \
     --values u32be --readahead 0 --pin \
     --progress-step 1000 \
     --output {}.out' ::: kmers.part.*.txt

cat kmers.part.*.txt.out | LC_ALL=C sort -u > hits.txt
```

> Tip: If you only need presence/absence (not counts), use `--values none` for maximum speed.

---

## Mask Syntax

| Symbol  | Meaning                                                        |
|---------|----------------------------------------------------------------|
| `?`     | Any character from the alphabet (wildcard)                     |
| `[AC]`  | Allowed set at this position (only listed characters allowed)  |
| `G`     | Exact match for the specified character                        |

Positions in masks are **1‑based**.

---

## Performance Notes

- The scanner chooses an optimal prefix length and bounds each scan range tightly with `(lo, hi)` to minimize iteration cost.
- Chunk retrieval uses robust point lookups (`cf.get` when available; iterator seek otherwise), buffered output, and no negative caching for RC to keep memory usage small and speed consistent.
- Works best when the DB is built with a fixed prefix extractor and prefix Bloom filters; `prefix_same_as_start` is enabled at read time.
- For pure point lookups, prefer `--readahead 0`, `--pin`, and fixed `--values` mode.

---

## Caveats

- If keys are not raw bytes or ASCII k‑mers, use `--no-raw` or provide inputs with `--chunks-hex`.
- DB should contain fixed‑length k‑mers as keys.
- Reverse‑complement lookup (`--rc`) is only supported for ASCII A/C/G/T inputs (not for hex mode).

---

## Command‑Line Options

Run:

```bash
kgrep --help
```

to see all options.

---

## License

MIT License — see repository for details.
