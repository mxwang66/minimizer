# Minimal btllib-style `indexlr` for minimizers

This repository provides a compact C++/Python implementation of an `indexlr`-like workflow:

```python
from btllib import indexlr
```

It reuses the existing ntHash/minimizer code in this repo and adds:
- Python bindings via `pybind11`
- FASTA parsing for plain `.fa/.fasta` and gzipped `.gz` inputs via `zlib`

## API

```python
def indexlr(
    assembly_path: Path,
    kmerlen: int,
    windowsize: int,
    assembly_idx: int,
    is_target: bool,
) -> tuple[bytes, tuple[str, ...]]:
    ...
```

The returned `bytes` buffer is laid out as records of:

- `hash` (`uint64`)
- `pos` (`uint32`)
- `record_idx` (`uint16`)
- `assembly_idx` (`uint16`)
- `is_target` (`bool` / 1 byte)

in exactly that order, with no struct padding.

## Conda setup (recommended)

```bash
conda create -n minimizer python=3.12 cmake ninja pybind11 zlib numpy -c conda-forge
conda activate minimizer
```

## Build

From repo root:

```bash
cmake -S . -B build -G Ninja \
  -DPython_EXECUTABLE="$(which python)" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

After build, the extension is created at:

- `build/btllib/_core*.so`

The Python wrapper package is in:

- `btllib/__init__.py`

Use either:

```bash
export PYTHONPATH="$PWD:$PWD/build"
```

or in Python:

```python
import sys
sys.path[:0] = ["/path/to/repo", "/path/to/repo/build"]
```

## Minimal Python example

```python
from pathlib import Path
import numpy as np
from btllib import indexlr

KMER_DTYPE = np.dtype([
    ("hash", np.uint64),
    ("pos", np.uint32),
    ("record_idx", np.uint16),
    ("assembly_idx", np.uint16),
    ("is_target", np.bool_),
])

kmers, idx_to_id = indexlr(
    assembly_path=Path("example.fa.gz"),
    kmerlen=31,
    windowsize=10,
    assembly_idx=0,
    is_target=True,
)

arr = np.frombuffer(kmers, dtype=KMER_DTYPE)
print(arr.shape, idx_to_id)
```

## Notes

- Record IDs are emitted for every FASTA record in order, even when a record contributes no minimizers.
- Minimizers are serialized record-by-record, preserving record and within-record minimizer order.
