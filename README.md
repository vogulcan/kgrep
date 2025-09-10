# kgrep

**Fast positional k-mer grep over RocksDict / RocksDB**

`kgrep` is a tool for scanning k-mers in a RocksDict/RocksDB database with position-specific filtering rules.  
It uses tight prefix-bounded scans and in-memory constraint checks to minimize I/O and maximize speed.

## Features

- Supports **mask** syntax (wildcards, allowed sets, exact matches).
- Supports **rule** syntax for position-specific constraints.
- `--rule` supports **OR** logic (multiple rule groups, matched if any group passes).
- Optimized RocksDB iteration with prefix bounds (`iterate_lower_bound`, `iterate_upper_bound`).
- Configurable enumeration limit to balance performance and completeness.
- DNA alphabet by default (`ACGT`), but can be overridden.
- Optional output to file.

## Installation

```bash
pip install "git+https://github.com/vogulcan/kgrep@master"
```

## Example usage

### 1. Search with mask syntax and count filter 

```bash
# Note: --mask overrides --rule
# Note: --mask is applied from the beginning of the k-mer, see Mask Syntax.
kgrep --db /path/dbfolder --mask "A[CT]?AA" --min-count 10
```

### 2. Search with rule syntax and output to file

```bash
# Note: --rule is 1-based.
kgrep --db /path/dbfolder --rule "1:A;2:CT" --min-count 1 --output output.txt
```

### 3. List all k-mers without constraints

```bash
kgrep --db /path/dbfolder
```

### 4. Complex multi-position rule

```bash
kgrep --db /path/dbfolder --rule "3:A;5:CG;10:T" --limit 1000
```

### 5. OR logic between rule groups

```bash
# Example: match (7:A AND 8:A) OR (17:T AND 18:T)
kgrep --db /path/dbfolder --rule "7:A;8:A|17:T;18:T"

# Same as above, using multiple --rule flags:
kgrep --db /path/dbfolder --rule "7:A;8:A" --rule "17:T;18:T"
```

### 6. Chunk retrieval (exact lookups for N kmers)
Retrieve values for an explicit list of kmers (all must have the same length):
```bash
# Inline ASCII kmers (comma-separated or multiple flags)
kgrep --db /path/dbfolder --chunks AAAAC --chunks TTTTG,CCCCA

# From a file (one kmer per line; lines starting with '#' are ignored)
kgrep --db /path/dbfolder --chunks-file kmers.txt

# From STDIN
cat kmers.txt | kgrep --db /path/dbfolder --chunks-file -

# Hex input (for non-ASCII keys); 0x prefix optional
kgrep --db /path/dbfolder --chunks-hex --chunks 0x4141414143 --chunks 7474745447
```
#### Notes for chunk retrieval:

- If `--k` is not provided, k is inferred from the provided kmers.
- If no kmers are provided, it falls back to the length of the first key in the DB (original behavior).
- Values are decoded per `--values` (default auto), and filtered with `--min-count` if a numeric value is available.
- Output for each found key is either:
    - `<key>\t<count>` when a numeric value is decoded, or
    - `<key>` if `--values none` is used or the value cannot be decoded.
- Non-ASCII keys are printed in hex.
- Missing keys are skipped silently (no output line).
- Duplicates are de-duplicated.
- `--limit` still applies to the number of printed lines.

## Mask Syntax

| Symbol  | Meaning |
| ------- | ---------------------------------------------------------------- |
| `?`     | Any character from the alphabet (wildcard)                       |
| `[AC]`  | Allowed set at this position (only listed characters allowed)    |
| `G`     | Exact match for the specified character                          |

Positions in masks are **1-based**.

## Performance Notes

- The script chooses an optimal prefix length to enumerate (`--max-enum` controls limit).
- Prefix scanning is bounded tightly (`lo, hi`) to minimize iteration cost.
- Works best with RocksDB prefix extractor enabled and `prefix_same_as_start` set.
- Chunk retrieval uses direct iterator seeks for each provided key (fast exact lookups) and respects the same decoding/filters.

## Caveats

- If keys are not raw bytes or ASCII k-mers, use `--no-raw` or provide inputs with `--chunks-hex`.
- DB should contain fixed-length k-mers as keys.
- Performance degrades without bound setters.

## Command-Line Options

Run:

```bash
kgrep --help
```

to see all options.

## License

MIT License — see repository for details.



