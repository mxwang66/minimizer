from __future__ import annotations

from os import PathLike
from pathlib import Path
from pkgutil import extend_path

__path__ = extend_path(__path__, __name__)

from ._core import indexlr_native


def indexlr(
    assembly_path: Path,
    kmerlen: int,
    windowsize: int,
    assembly_idx: int,
    is_target: bool,
) -> tuple[bytes, tuple[str, ...]]:
    """Compute minimizers for all FASTA records in order."""

    path_str = str(Path(assembly_path)) if isinstance(assembly_path, PathLike) else str(assembly_path)
    kmers, idx_to_id = indexlr_native(
        path_str,
        int(kmerlen),
        int(windowsize),
        int(assembly_idx),
        bool(is_target),
    )
    return kmers, tuple(idx_to_id)


__all__ = ["indexlr"]
