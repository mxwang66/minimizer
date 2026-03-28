from __future__ import annotations

from os import PathLike
from pathlib import Path
from pkgutil import extend_path

__path__ = extend_path(__path__, __name__)

from ._core import indexlr_native


def indexlr(
    assembly_path: list[Path],
    kmerlen: int,
    windowsize: int,
    assembly_idx: list[int],
    is_target: list[bool],
) -> tuple[
    bytes,
    list[tuple[str, ...]],
]:
    """Compute minimizers for all FASTA records in each assembly in order."""

    if len(assembly_path) != len(assembly_idx) or len(assembly_path) != len(is_target):
        raise ValueError("assembly_path, assembly_idx, and is_target must have the same length")

    path_strs = [
        str(Path(path)) if isinstance(path, PathLike) else str(path) for path in assembly_path
    ]
    kmers, idx_to_id = indexlr_native(
        path_strs,
        int(kmerlen),
        int(windowsize),
        [int(idx) for idx in assembly_idx],
        [bool(target) for target in is_target],
    )
    return kmers, [tuple(ids) for ids in idx_to_id]


__all__ = ["indexlr"]
