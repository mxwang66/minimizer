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
    assembly_paths: Iterable[Path],
    kmerlen: int,
    windowsize: int,
    assembly_idx: Iterable[int],
    is_target: Iterable[bool],
) -> tuple[
    np.ndarray,
    list[tuple[str, ...]],
    np.ndarray[np.uint64],
    np.ndarray[np.uint64],
]:
    ...
```

The returned first `numpy.ndarray` is a 1-D structured array with dtype:

- `hash` (`uint64`)
- `pos` (`uint32`)
- `record_idx` (`uint16`)
- `assembly_idx` (`uint16`)
- `is_target` (`bool`)

in exactly that order (equivalent to `btllib.KMER_DTYPE`).

`assembly_paths`, `assembly_idx`, and `is_target` are parallel iterables.

The function also returns:
- `record_offsets` (`np.ndarray[np.uint64]`): global minimizer indices where a new FASTA record starts contributing
  minimizers (flattened across assemblies). Records with zero minimizers are omitted.
- `assembly_offsets` (`np.ndarray[np.uint64]`): global minimizer indices where a new assembly starts contributing
  minimizers. Assemblies with zero minimizers are omitted.

## Conda setup (recommended)

```bash
conda create -n minimizer python cmake ninja cxx-compiler pybind11 zlib numpy -c conda-forge
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
from btllib import KMER_DTYPE, indexlr

kmers, idx_to_id, record_offsets, assembly_offsets = indexlr(
    assembly_paths=[Path("example.fa.gz"), Path("example2.fa.gz")],
    kmerlen=31,
    windowsize=10,
    assembly_idx=[0, 1],
    is_target=[True, False],
)

print(kmers.dtype == KMER_DTYPE, kmers.shape)
print(idx_to_id, record_offsets, assembly_offsets)
```

## Notes

- Record IDs are emitted for every FASTA record in order for each assembly, even when a record contributes no minimizers.
- Minimizers are serialized record-by-record, preserving record and within-record minimizer order.
- `record_offsets` and `assembly_offsets` are minimizer-index offsets.
- The native `indexlr` compute path releases the Python GIL, so concurrent Python threads can call `btllib.indexlr(...)` at the same time.
