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
    assembly_path: list[Path],
    kmerlen: int,
    windowsize: int,
    assembly_idx: list[int],
    is_target: list[bool],
) -> tuple[
    np.ndarray,
    list[tuple[str, ...]],
    list[int],
    list[int],
]:
    ...
```

The returned `numpy.ndarray` (dtype `uint8`) buffer is laid out as records of:

- `hash` (`uint64`)
- `pos` (`uint32`)
- `record_idx` (`uint16`)
- `assembly_idx` (`uint16`)
- `is_target` (`bool` / 1 byte)

in exactly that order, with no struct padding.

`assembly_path`, `assembly_idx`, and `is_target` are parallel lists. Minimizers for each
assembly are computed and serialized exactly as before, and each assembly's bytes are appended
to one shared output byte buffer in input order. The returned NumPy buffer is writable.

The function also returns:
- `record_offsets`: global minimizer indices where a new FASTA record starts contributing
  minimizers (flattened across assemblies). Records with zero minimizers are omitted.
- `assembly_offsets`: global minimizer indices where a new assembly starts contributing
  minimizers. Assemblies with zero minimizers are omitted.

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

kmers_u8, idx_to_id, record_offsets, assembly_offsets = indexlr(
    assembly_path=[Path("example.fa.gz"), Path("example2.fa.gz")],
    kmerlen=31,
    windowsize=10,
    assembly_idx=[0, 1],
    is_target=[True, False],
)

arr = kmers_u8.view(KMER_DTYPE)
print(arr.shape, idx_to_id, record_offsets, assembly_offsets)
```

## Notes

- Record IDs are emitted for every FASTA record in order for each assembly, even when a record contributes no minimizers.
- Minimizers are serialized record-by-record, preserving record and within-record minimizer order.
- `record_offsets` and `assembly_offsets` are minimizer-index offsets (not byte offsets). Byte offset can be computed as `offset * 17`.
