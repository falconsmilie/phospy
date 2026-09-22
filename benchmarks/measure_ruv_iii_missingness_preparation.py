"""Benchmark RUV-III missingness preparation stages and kernel context."""
# ruff: noqa: E402

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
DEFAULT_REPORT_PATH = ROOT / ("benchmarks/reports/ruv-iii-missingness-preparation.json")
for _path in (ROOT, SRC):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from phospy.science.batch_correction.ruv_iii import (
    RuvIIIKernel,
    RuvIIIReplicateStructure,
)
from phospy.science.batch_correction.ruv_iii_executor import (
    _materialize_observation_mask,
    _prepare_matrix,
    _temporary_row_median_completion,
    _validate_actual_vs_governed_missingness,
)
from phospy.science.configs.preprocessing import TemporaryImputationMethod

TierName = Literal["small", "representative", "stress"]
CaseName = Literal["complete", "actual_missing", "upstream_imputed", "mixed"]


@dataclass(frozen=True, slots=True)
class BenchmarkTier:
    rows: int
    columns: int
    governed_missing_fraction: float
    default_repeats: int


TIERS: dict[TierName, BenchmarkTier] = {
    "small": BenchmarkTier(1_000, 20, 0.05, 5),
    "representative": BenchmarkTier(10_000, 50, 0.15, 3),
    "stress": BenchmarkTier(30_000, 100, 0.30, 1),
}
CASES: tuple[CaseName, ...] = (
    "complete",
    "actual_missing",
    "upstream_imputed",
    "mixed",
)


@dataclass(frozen=True, slots=True)
class _BenchmarkObservationMask:
    feature_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    originally_missing_cells: tuple[tuple[str, str], ...]

    def to_payload(self) -> dict[str, object]:
        return {
            "feature_ids": list(self.feature_ids),
            "sample_ids": list(self.sample_ids),
            "originally_missing_cells": [
                list(cell) for cell in self.originally_missing_cells
            ],
        }


@dataclass(frozen=True, slots=True)
class _BenchmarkTemporaryImputationPolicy:
    allowed: bool
    method: TemporaryImputationMethod
    method_parameters: tuple[tuple[str, int], ...]

    def to_payload(self) -> dict[str, object]:
        return {
            "allowed": self.allowed,
            "method": self.method.value,
            "method_parameters": dict(self.method_parameters),
        }


@dataclass(frozen=True, slots=True)
class _BenchmarkPlan:
    observation_mask: _BenchmarkObservationMask
    temporary_imputation_policy: _BenchmarkTemporaryImputationPolicy


@dataclass(frozen=True, slots=True)
class BenchmarkInput:
    matrix: pd.DataFrame
    plan: _BenchmarkPlan


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tier",
        choices=(*TIERS, "all"),
        default="small",
        help="run one tier or all tiers; stress is never selected implicitly",
    )
    parser.add_argument(
        "--case",
        choices=(*CASES, "all"),
        default="all",
        help="run one missingness case or all cases",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=None,
        help="override the tier repeat count",
    )
    parser.add_argument(
        "--skip-kernel",
        action="store_true",
        help="skip the separately reported RUV-III kernel context measurement",
    )
    parser.add_argument(
        "--phase",
        choices=("baseline", "final", "ad-hoc"),
        default="ad-hoc",
        help="label the implementation phase represented by the evidence",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "optional JSON result path (use the default scratch location "
            f"{DEFAULT_REPORT_PATH}); JSON is always emitted to stdout"
        ),
    )
    return parser.parse_args()


def _site_key(position: int) -> str:
    return (
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        f"protein_identifier=P{position:06d}|residue=S|position={position}"
    )


def build_benchmark_input(
    *, tier_name: TierName, case_name: CaseName
) -> BenchmarkInput:
    tier = TIERS[tier_name]
    seed = 72_001 + tuple(TIERS).index(tier_name) * 100 + CASES.index(case_name)
    rng = np.random.default_rng(seed)
    feature_ids = tuple(_site_key(position + 1) for position in range(tier.rows))
    sample_ids = tuple(f"sample_{position + 1:03d}" for position in range(tier.columns))
    values = rng.normal(loc=16.0, scale=2.5, size=(tier.rows, tier.columns))
    complete = pd.DataFrame(
        values,
        index=pd.Index(feature_ids, name="site_key"),
        columns=pd.Index(sample_ids, name="sample_id"),
    )

    governed_count = (
        0
        if case_name == "complete"
        else int(round(values.size * tier.governed_missing_fraction))
    )
    flat_positions = (
        np.empty(0, dtype=np.int64)
        if governed_count == 0
        else np.sort(rng.choice(values.size, size=governed_count, replace=False))
    )
    row_positions, column_positions = np.unravel_index(flat_positions, values.shape)
    governed_cells = tuple(
        (feature_ids[int(row)], sample_ids[int(column)])
        for row, column in zip(row_positions, column_positions, strict=True)
    )

    matrix = complete.copy(deep=True)
    if case_name == "actual_missing":
        matrix_values = matrix.to_numpy(copy=True)
        matrix_values[row_positions, column_positions] = np.nan
        matrix = pd.DataFrame(
            matrix_values,
            index=complete.index.copy(),
            columns=complete.columns.copy(),
        )
    elif case_name == "mixed":
        actual_count = governed_count // 2
        matrix_values = matrix.to_numpy(copy=True)
        matrix_values[row_positions[:actual_count], column_positions[:actual_count]] = (
            np.nan
        )
        matrix = pd.DataFrame(
            matrix_values,
            index=complete.index.copy(),
            columns=complete.columns.copy(),
        )
    observed_counts = matrix.notna().to_numpy(dtype=bool).sum(axis=1)
    if bool((observed_counts < 2).any()):
        raise AssertionError("benchmark generation produced an ineligible row")

    mask = _BenchmarkObservationMask(feature_ids, sample_ids, governed_cells)
    allowed_policy = _BenchmarkTemporaryImputationPolicy(
        allowed=True,
        method=TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY,
        method_parameters=(("min_observed_values", 2),),
    )
    return BenchmarkInput(
        matrix=matrix,
        plan=_BenchmarkPlan(mask, allowed_policy),
    )


def _median_runtime(callable_: object, *, repeats: int) -> float:
    timings: list[float] = []
    for _ in range(repeats):
        started = time.perf_counter()
        callable_()  # type: ignore[operator]
        timings.append(time.perf_counter() - started)
    return float(statistics.median(timings))


def _median_runtime_with_setup(
    setup: object,
    callable_: object,
    *,
    repeats: int,
) -> float:
    timings: list[float] = []
    for _ in range(repeats):
        argument = setup()  # type: ignore[operator]
        started = time.perf_counter()
        callable_(argument)  # type: ignore[operator]
        timings.append(time.perf_counter() - started)
    return float(statistics.median(timings))


def _current_process_rss_mib() -> float | None:
    if os.name == "posix":
        try:
            resident_pages = int(Path("/proc/self/statm").read_text().split()[1])
            return resident_pages * os.sysconf("SC_PAGE_SIZE") / (1024 * 1024)
        except (OSError, ValueError, IndexError):
            return None
    if os.name != "nt":
        return None
    try:
        import ctypes
        from ctypes import wintypes

        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(ProcessMemoryCounters)
        kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
        psapi = ctypes.windll.psapi  # type: ignore[attr-defined]
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = (
            wintypes.HANDLE,
            ctypes.POINTER(ProcessMemoryCounters),
            wintypes.DWORD,
        )
        psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
        handle = kernel32.GetCurrentProcess()
        ok = psapi.GetProcessMemoryInfo(
            handle,
            ctypes.byref(counters),
            counters.cb,
        )
        if not ok:
            return None
        return float(int(counters.WorkingSetSize)) / (1024 * 1024)
    except (AttributeError, OSError):
        return None


def _preparation_rss_mib(
    benchmark_input: BenchmarkInput,
) -> tuple[float | None, float | None]:
    samples: list[float] = []
    stop = threading.Event()
    starting_rss = _current_process_rss_mib()
    if starting_rss is not None:
        samples.append(starting_rss)

    def sample_memory() -> None:
        while not stop.is_set():
            value = _current_process_rss_mib()
            if value is not None:
                samples.append(value)
            stop.wait(0.002)

    sampler = threading.Thread(target=sample_memory, daemon=True)
    sampler.start()
    try:
        _prepare_matrix(phospho=benchmark_input.matrix, plan=benchmark_input.plan)
    finally:
        stop.set()
        sampler.join()
    final_value = _current_process_rss_mib()
    if final_value is not None:
        samples.append(final_value)
    if not samples:
        return None, None
    peak_rss = max(samples)
    return peak_rss, (
        None if starting_rss is None else max(0.0, peak_rss - starting_rss)
    )


def _kernel_seconds(benchmark_input: BenchmarkInput, *, repeats: int) -> float:
    prepared = _prepare_matrix(
        phospho=benchmark_input.matrix,
        plan=benchmark_input.plan,
    )
    sample_ids = tuple(str(value) for value in prepared.working.columns)
    replicate_by_sample = {
        sample_id: f"replicate_{position // 2 + 1:03d}"
        for position, sample_id in enumerate(sample_ids)
    }
    replicate_structure = RuvIIIReplicateStructure.from_assignments(
        sample_order=sample_ids,
        replicate_by_sample=replicate_by_sample,
    )
    control_count = min(64, prepared.working.shape[0])
    controls = tuple(str(value) for value in prepared.working.index[:control_count])
    kernel = RuvIIIKernel()
    return _median_runtime(
        lambda: kernel.run(
            phospho=prepared.working,
            control_site_keys=controls,
            replicate_structure=replicate_structure,
            k=1,
        ),
        repeats=repeats,
    )


def run_case(
    *,
    tier_name: TierName,
    case_name: CaseName,
    repeats: int,
    include_kernel: bool,
) -> dict[str, object]:
    benchmark_input = build_benchmark_input(
        tier_name=tier_name,
        case_name=case_name,
    )
    matrix = benchmark_input.matrix
    plan = benchmark_input.plan
    values = matrix.to_numpy(dtype="float64", copy=True)
    originally_missing = _materialize_observation_mask(
        phospho=matrix,
        plan=plan,
    )
    actual_missing = _validate_actual_vs_governed_missingness(
        values=values,
        originally_missing=originally_missing,
    )
    mask_seconds = _median_runtime(
        lambda: _materialize_observation_mask(phospho=matrix, plan=plan),
        repeats=repeats,
    )
    validation_seconds = _median_runtime(
        lambda: _validate_actual_vs_governed_missingness(
            values=values,
            originally_missing=originally_missing,
        ),
        repeats=repeats,
    )
    completion_seconds = (
        _median_runtime_with_setup(
            values.copy,
            lambda fresh_values: _temporary_row_median_completion(
                values=fresh_values,
                actual_missing=actual_missing,
                index=matrix.index,
                policy=plan.temporary_imputation_policy,
            ),
            repeats=repeats,
        )
        if bool(actual_missing.any())
        else 0.0
    )
    complete_preparation_seconds = _median_runtime(
        lambda: _prepare_matrix(
            phospho=matrix,
            plan=plan,
        ),
        repeats=repeats,
    )

    tier = TIERS[tier_name]
    governed_count = len(benchmark_input.plan.observation_mask.originally_missing_cells)
    actual_count = int(benchmark_input.matrix.isna().to_numpy().sum())
    peak_rss_mib, peak_rss_delta_mib = _preparation_rss_mib(benchmark_input)
    result: dict[str, object] = {
        "tier": tier_name,
        "case": case_name,
        "rows": tier.rows,
        "columns": tier.columns,
        "governed_missing_fraction": tier.governed_missing_fraction,
        "governed_missing_cell_count": governed_count,
        "actual_missing_cell_count": actual_count,
        "upstream_imputed_cell_count": governed_count - actual_count,
        "repeats": repeats,
        "observation_mask_materialisation_seconds": mask_seconds,
        "actual_vs_governed_validation_seconds": validation_seconds,
        "temporary_row_median_completion_seconds": completion_seconds,
        "complete_matrix_preparation_seconds": complete_preparation_seconds,
        "peak_process_rss_mib": peak_rss_mib,
        "preparation_peak_rss_delta_mib": peak_rss_delta_mib,
        "ruv_iii_kernel_seconds": (
            _kernel_seconds(benchmark_input, repeats=repeats)
            if include_kernel
            else None
        ),
    }
    result["mask_fraction_of_preparation"] = mask_seconds / complete_preparation_seconds
    result["completion_fraction_of_preparation"] = (
        completion_seconds / complete_preparation_seconds
    )
    return result


def _dependency_version(distribution: str) -> str:
    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        return "unavailable"


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git_output(*arguments: str) -> bytes | None:
    try:
        completed = subprocess.run(
            ("git", *arguments),
            cwd=ROOT,
            check=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return completed.stdout


def _source_provenance() -> dict[str, object]:
    relative_sources = (
        "src/phospy/science/batch_correction/_missingness_preparation.py",
        "src/phospy/science/batch_correction/executor.py",
        "src/phospy/science/batch_correction/ruv_iii_executor.py",
    )
    source_hashes = {
        relative_path: _sha256((ROOT / relative_path).read_bytes())
        for relative_path in relative_sources
        if (ROOT / relative_path).is_file()
    }
    head = _git_output("rev-parse", "HEAD")
    status = _git_output("status", "--porcelain", "--untracked-files=no")
    diff = _git_output("diff", "--binary", "HEAD", "--", *relative_sources)
    return {
        "command": [sys.executable, *sys.argv],
        "git_head": None if head is None else head.decode().strip(),
        "git_worktree_dirty": None if status is None else bool(status.strip()),
        "tracked_source_diff_sha256": None if diff is None else _sha256(diff),
        "production_source_sha256": source_hashes,
    }


def main() -> int:
    args = _parse_args()
    tier_names: tuple[TierName, ...] = (
        tuple(TIERS) if args.tier == "all" else (args.tier,)
    )
    case_names: tuple[CaseName, ...] = CASES if args.case == "all" else (args.case,)
    results: list[dict[str, object]] = []
    for tier_name in tier_names:
        repeats = (
            TIERS[tier_name].default_repeats if args.repeats is None else args.repeats
        )
        if repeats < 1:
            raise SystemExit("--repeats must be at least 1")
        for case_name in case_names:
            results.append(
                run_case(
                    tier_name=tier_name,
                    case_name=case_name,
                    repeats=repeats,
                    include_kernel=not args.skip_kernel,
                )
            )

    payload = {
        "schema_version": 1,
        "benchmark": "ruv_iii_missingness_preparation",
        "phase": args.phase,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "stage_measurement_method": (
            "direct repeated calls to the private production materialisation, "
            "validation, and completion stages; completion setup copies are "
            "performed before timing"
        ),
        "kernel_context_excluded_from_preparation": True,
        "source_provenance": _source_provenance(),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": _dependency_version("numpy"),
            "pandas": _dependency_version("pandas"),
        },
        "results": results,
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True)
    print(f"benchmark_result_json={json.dumps(payload, sort_keys=True)}")
    if args.output is not None:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
