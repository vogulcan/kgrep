# kgrep

**Fast positional k-mer grep over RocksDict / RocksDB**

`kgrep` is a high-performance tool for scanning k-mers in a RocksDict/RocksDB database with position-specific filtering rules.  
It uses tight prefix-bounded scans and in-memory constraint checks to minimize I/O and maximize speed.

## Features

- Supports **mask** syntax (wildcards, allowed sets, exact matches).
- Supports **rule** syntax for position-specific constraints.
- Optimized RocksDB iteration with prefix bounds (`iterate_lower_bound`, `iterate_upper_bound`).
- Configurable enumeration limit to balance performance and completeness.
- Multiple value decoding modes:
  - Big-endian unsigned 32-bit
  - Little-endian unsigned 64-bit
  - ASCII integer fallback
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

## Caveats

- If keys are not raw bytes or ASCII k-mers, use `--no-raw`.
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
