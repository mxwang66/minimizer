from __future__ import annotations

import hashlib
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from btllib import indexlr


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes", "y", "t"}:
        return True
    if normalized in {"false", "0", "no", "n", "f"}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value!r}")


def load_assemblies(csv_path: Path) -> tuple[list[str], list[int], list[bool]]:
    lines = csv_path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError("assemblies.csv must include a header and at least one data row")

    rows = [line.split(",") for line in lines[1:] if line.strip()]
    if not rows:
        raise ValueError("assemblies.csv has no data rows")

    idx = [int(row[0].strip()) for row in rows]
    paths = [row[1].strip() for row in rows]
    is_tar = [parse_bool(row[2]) for row in rows]
    return paths, idx, is_tar


def run_indexlr(paths: list[str], idx: list[int], is_tar: list[bool]):
    return indexlr(
        assembly_paths=paths,
        kmerlen=21,
        windowsize=200,
        assembly_idx=idx,
        is_target=is_tar,
    )


def kmers_digest(kmers) -> str:
    return hashlib.sha256(kmers.view("u1").tobytes()).hexdigest()


def assert_matches_baseline(name: str, baseline, candidate) -> None:
    base_kmers, base_idx_to_id, base_record_offsets, base_assembly_offsets = baseline
    cand_kmers, cand_idx_to_id, cand_record_offsets, cand_assembly_offsets = candidate

    assert cand_kmers.shape == base_kmers.shape, f"{name}: kmers.shape mismatch"
    assert kmers_digest(cand_kmers) == kmers_digest(base_kmers), f"{name}: kmers digest mismatch"
    assert cand_idx_to_id == base_idx_to_id, f"{name}: idx_to_id mismatch"
    assert cand_record_offsets.tolist() == base_record_offsets.tolist(), f"{name}: record_offsets mismatch"
    assert cand_assembly_offsets.tolist() == base_assembly_offsets.tolist(), f"{name}: assembly_offsets mismatch"


def main() -> None:
    csv_path = Path("assemblies.csv")
    paths, idx, is_tar = load_assemblies(csv_path)

    print(f"Loaded {len(paths)} assemblies from {csv_path}")

    t0 = time.perf_counter()
    baseline = run_indexlr(paths, idx, is_tar)
    serial_once_s = time.perf_counter() - t0

    base_kmers, base_idx_to_id, base_record_offsets, base_assembly_offsets = baseline
    print(
        "Baseline:",
        base_kmers.shape,
        len(base_idx_to_id),
        base_record_offsets.shape,
        base_assembly_offsets.shape,
        f"sha256={kmers_digest(base_kmers)}",
    )

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(run_indexlr, paths, idx, is_tar) for _ in range(2)]
        threaded_results = [f.result() for f in futures]

    for i, threaded_result in enumerate(threaded_results, start=1):
        assert_matches_baseline(f"threaded_run_{i}", baseline, threaded_result)
    print("Two concurrent runs matched baseline output")

    repeat_factor = 6
    heavy_paths = paths * repeat_factor
    heavy_idx = idx * repeat_factor
    heavy_is_tar = is_tar * repeat_factor
    print(f"Timing smoke check with repeat_factor={repeat_factor} ({len(heavy_paths)} assemblies)")

    t1 = time.perf_counter()
    _ = run_indexlr(heavy_paths, heavy_idx, heavy_is_tar)
    _ = run_indexlr(heavy_paths, heavy_idx, heavy_is_tar)
    serial_s = time.perf_counter() - t1

    t2 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_a = executor.submit(run_indexlr, heavy_paths, heavy_idx, heavy_is_tar)
        future_b = executor.submit(run_indexlr, heavy_paths, heavy_idx, heavy_is_tar)
        _ = future_a.result()
        _ = future_b.result()
    parallel_s = time.perf_counter() - t2

    print(f"Single baseline call: {serial_once_s:.3f}s")
    print(f"Serial heavy (2 calls): {serial_s:.3f}s")
    print(f"Parallel heavy (2 threads): {parallel_s:.3f}s")

    if parallel_s >= serial_s * 0.95:
        print("WARNING: Limited overlap observed; workload or environment may be too small.")
    else:
        print("Observed overlap between concurrent indexlr calls.")


if __name__ == "__main__":
    main()
