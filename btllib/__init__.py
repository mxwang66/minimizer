from __future__ import annotations

from os import PathLike
from pathlib import Path
from pkgutil import extend_path

import numpy as np
from numpy.typing import NDArray

__path__ = extend_path(__path__, __name__)

from ._core import indexlr_native

KMER_DTYPE = np.dtype([
    ('hash', np.uint64), 
    ('pos', np.uint32), 
    ('record_idx', np.uint16), 
    ('assembly_idx', np.uint16), 
    ('is_target', np.bool_), 
])


def indexlr(
    assembly_path: list[Path],
    kmerlen: int,
    windowsize: int,
    assembly_idx: list[int],
    is_target: list[bool],
) -> tuple[
    NDArray[np.uint8],
    list[tuple[str, ...]],
    NDArray[np.uint64],
    NDArray[np.uint64],
]:
    """Compute minimizers for all FASTA records in each assembly in order."""

    if len(assembly_path) != len(assembly_idx) or len(assembly_path) != len(is_target):
        raise ValueError("assembly_path, assembly_idx, and is_target must have the same length")

    path_strs = [
        str(Path(path)) if isinstance(path, PathLike) else str(path) for path in assembly_path
    ]
    kmers, idx_to_id, record_offsets, assembly_offsets = indexlr_native(
        path_strs,
        int(kmerlen),
        int(windowsize),
        [int(idx) for idx in assembly_idx],
        [bool(target) for target in is_target],
    )
    return kmers.view(KMER_DTYPE), [tuple(ids) for ids in idx_to_id], record_offsets, assembly_offsets


__all__ = ["indexlr"]
