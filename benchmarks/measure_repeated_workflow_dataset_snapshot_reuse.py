#!/usr/bin/env python3
"""Benchmark repeated workflow reuse of dataset-owned immutable snapshots.

Targets:
- `phospy.workflows.differential.DifferentialAnalysisWorkflow.run`
- `phospy.workflows.kinase.KinaseWorkflow.run`

This is an explicitly invoked local benchmark. It records machine-specific
runtime, tracemalloc peak memory, full-frame pandas deep-copy counts, and
dataset internal snapshot construction counts for repeated workflow use of the
same unchanged dataset.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
import tracemalloc
from collections import Counter
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

REPORTS_DIRECTORY = Path("benchmarks/reports")
REPORT_PATH = REPORTS_DIRECTORY / "repeated-workflow-dataset-snapshot-reuse.json"
DOCUMENTED_MAIN_METRIC_KEYS = (
    "dataset_phospho_rows",
    "dataset_phospho_columns",
    "dataset_site_metadata_columns",
    "frame_dtypes_json",
    "differential_first_run_seconds",
    "differential_repeated_run_seconds",
    "differential_first_peak_tracemalloc_mib",
    "differential_repeated_peak_tracemalloc_mib",
    "kinase_first_run_seconds",
    "kinase_repeated_run_seconds",
    "kinase_first_peak_tracemalloc_mib",
    "kinase_repeated_peak_tracemalloc_mib",
    "full_frame_deep_copy_counts_json",
    "snapshot_construction_counts_json",
    "environment_json",
    "dependency_versions_json",
)


@dataclass(frozen=True, slots=True)
class BenchmarkConfig:
    n_sites: int = 240
    n_samples: int = 8
    n_kinases: int = 8
    substrates_per_kinase: int = 24


@dataclass(frozen=True, slots=True)
class RunMeasurement:
    runtime_seconds: float
    peak_tracemalloc_mib: float
    full_frame_deep_copies: dict[str, int]
    snapshot_constructions: dict[str, int]


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    config: BenchmarkConfig
    dataset_dimensions: dict[str, int]
    frame_dtypes: dict[str, dict[str, str]]
    differential_first: RunMeasurement
    differential_repeated: RunMeasurement
    kinase_first: RunMeasurement
    kinase_repeated: RunMeasurement
    total_full_frame_deep_copies: dict[str, int]
    total_snapshot_constructions: dict[str, int]
    environment: dict[str, str]
    dependencies: dict[str, str]


@dataclass(frozen=True, slots=True)
class _DatasetFrameSignature:
    label: str
    shape: tuple[int, int]
    index: tuple[object, ...]
    columns: tuple[object, ...]


@dataclass(slots=True)
class _InstrumentationCounts:
    full_frame_deep_copies: Counter[str]
    snapshot_constructions: Counter[str]


def default_config() -> BenchmarkConfig:
    return BenchmarkConfig()


def default_report_path() -> Path:
    return ROOT / REPORT_PATH


def _dependency_version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "unavailable"


def _environment() -> dict[str, str]:
    return {
        "python_version": sys.version.replace("\n", " "),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
    }


def _dependencies() -> dict[str, str]:
    return {
        "phospy": _dependency_version("phospy"),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": _dependency_version("scipy"),
    }


def _json_metric(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _site_sequence(site: str) -> str:
    residue = str(site).strip().upper()[0]
    return ("A" * 15) + residue + ("A" * 15)


def _build_dataset(config: BenchmarkConfig):
    from phospy.api import Organism
    from tests.support.analysis_ready_dataset_factories import (
        trusted_analysis_ready_dataset_from_tables,
    )
    from tests.support.intensity_scale_states import (
        supported_log2_intensity_scale_state,
        supported_log2_processing_state,
    )
    from tests.support.site_keys import (
        protein_site_key_index,
        site_key_context_columns,
    )

    genes = [f"GENE{i:04d}" for i in range(int(config.n_sites))]
    residues = ("S", "T", "Y")
    sites = [
        f"{residues[position % len(residues)]}{position + 1}"
        for position in range(int(config.n_sites))
    ]
    site_index = protein_site_key_index(protein_identifiers=genes, sites=sites)
    sample_ids = [
        f"{condition}_{replicate + 1}"
        for condition in ("A", "B")
        for replicate in range(int(config.n_samples) // 2)
    ]
    if len(sample_ids) != int(config.n_samples):
        raise ValueError("n_samples must be even")
    row_effect = np.linspace(0.0, 1.0, int(config.n_sites), dtype=np.float64)
    columns: dict[str, np.ndarray] = {}
    for sample_position, sample_id in enumerate(sample_ids):
        condition_effect = 0.0 if sample_id.startswith("A_") else 0.35
        replicate_offset = float(sample_position % (int(config.n_samples) // 2)) * 0.03
        columns[sample_id] = (
            10.0
            + row_effect
            + condition_effect
            + replicate_offset
            + (np.arange(int(config.n_sites), dtype=np.float64) % 7.0) * 0.01
        )
    phospho = pd.DataFrame(columns, index=site_index.copy())
    display_ids = [f"{gene};{site};" for gene, site in zip(genes, sites, strict=True)]
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": genes,
            "site": sites,
            "site_sequence": [_site_sequence(site) for site in sites],
            "protein_id": genes,
            "localisation_confidence": np.full(int(config.n_sites), 0.95),
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )


def _differential_request(dataset: Any):
    from phospy.api import (
        Contrast,
        DifferentialAnalysisRequest,
        ExperimentalDesign,
        SampleDesignRecord,
    )

    sample_ids = [str(sample_id) for sample_id in dataset.phospho.columns]
    half = len(sample_ids) // 2
    records = [
        SampleDesignRecord(
            sample_id=sample_id,
            condition="A" if position < half else "B",
            biological_replicate_id=f"{sample_id}_bio",
        )
        for position, sample_id in enumerate(sample_ids)
    ]
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=ExperimentalDesign(samples=tuple(records)),
        contrasts=(
            Contrast(
                name="B_vs_A",
                numerator_condition="B",
                denominator_condition="A",
            ),
        ),
    )


def _kinase_request(dataset: Any, config: BenchmarkConfig):
    from phospy.advanced import (
        KinasePredictionConfig,
        KinaseScoringConfig,
        ReferenceContextCompatibilityPolicy,
    )
    from phospy.api import (
        KinaseWorkflowRequest,
        ReferenceBundle,
    )

    site_metadata = dataset.site_metadata
    display_ids = [str(value) for value in site_metadata["display_id"].tolist()]
    rows: list[dict[str, str]] = []
    site_count = len(display_ids)
    for kinase_position in range(int(config.n_kinases)):
        kinase = f"KINASE_{kinase_position + 1:03d}"
        for offset in range(int(config.substrates_per_kinase)):
            site_id = display_ids[(kinase_position * 7 + offset) % site_count]
            rows.append({"kinase": kinase, "substrate_site": site_id})
    reference_sequences = pd.DataFrame(
        {"site_sequence": site_metadata["site_sequence"].to_numpy(dtype=str)},
        index=pd.Index(display_ids, name="site_id"),
    )
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=ReferenceBundle(
            organism=dataset.organism,
            kinase_substrate_map=pd.DataFrame(rows),
            site_sequences=reference_sequences,
        ),
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=4,
            deterministic_max_selected_kinases=min(8, int(config.n_kinases)),
            adaptive_ensemble_runs=4,
        ),
        activity_config=None,
    )


def _frame_signatures(dataset: Any) -> tuple[_DatasetFrameSignature, ...]:
    frames = (
        ("dataset.phospho", dataset._phospho),
        ("dataset.site_metadata", dataset._site_metadata),
        ("dataset.sample_metadata", dataset._sample_metadata),
        ("dataset.total", dataset._total),
        ("dataset.comparisons", dataset._comparisons),
    )
    return tuple(
        _DatasetFrameSignature(
            label=label,
            shape=frame.shape,
            index=tuple(frame.index.tolist()),
            columns=tuple(frame.columns.tolist()),
        )
        for label, frame in frames
        if frame is not None
    )


def _matching_frame_label(
    frame: pd.DataFrame,
    signatures: tuple[_DatasetFrameSignature, ...],
) -> str | None:
    frame_shape = frame.shape
    frame_columns = tuple(frame.columns.tolist())
    for signature in signatures:
        if frame_shape != signature.shape or frame_columns != signature.columns:
            continue
        if tuple(frame.index.tolist()) == signature.index:
            return signature.label
    return None


@contextmanager
def _instrument_dataset(dataset: Any) -> Iterator[_InstrumentationCounts]:
    import phospy.science.datasets.internal_frame_store as internal_frame_store_module

    signatures = _frame_signatures(dataset)
    counts = _InstrumentationCounts(
        full_frame_deep_copies=Counter(),
        snapshot_constructions=Counter(),
    )
    original_copy = pd.DataFrame.copy
    original_snapshot = internal_frame_store_module.immutable_dataframe_snapshot

    def wrapped_copy(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        deep = kwargs.get("deep", args[0] if args else True)
        if bool(deep):
            label = _matching_frame_label(self, signatures)
            if label is not None:
                counts.full_frame_deep_copies.update((label,))
        return original_copy(self, *args, **kwargs)

    def counting_snapshot(
        value: pd.DataFrame,
        *,
        field_name: str,
        error_type: type[Exception] = TypeError,
    ):
        counts.snapshot_constructions.update((field_name,))
        return original_snapshot(
            value,
            field_name=field_name,
            error_type=error_type,
        )

    pd.DataFrame.copy = wrapped_copy  # pyright: ignore[reportAttributeAccessIssue] - benchmark instrumentation monkeypatches pandas copy counts.
    internal_frame_store_module.immutable_dataframe_snapshot = counting_snapshot
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_copy  # pyright: ignore[reportAttributeAccessIssue] - restore benchmark instrumentation.
        internal_frame_store_module.immutable_dataframe_snapshot = original_snapshot


def _counter_delta(after: Counter[str], before: Counter[str]) -> dict[str, int]:
    keys = set(after) | set(before)
    return {
        key: int(after[key] - before[key])
        for key in sorted(keys)
        if int(after[key] - before[key]) != 0
    }


def _measure_run(
    run: Callable[[], object],
    counts: _InstrumentationCounts,
) -> RunMeasurement:
    full_copy_before = counts.full_frame_deep_copies.copy()
    snapshot_before = counts.snapshot_constructions.copy()
    tracemalloc.start()
    started = time.perf_counter()
    run()
    runtime_seconds = time.perf_counter() - started
    _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return RunMeasurement(
        runtime_seconds=float(runtime_seconds),
        peak_tracemalloc_mib=float(peak_bytes / (1024.0 * 1024.0)),
        full_frame_deep_copies=_counter_delta(
            counts.full_frame_deep_copies,
            full_copy_before,
        ),
        snapshot_constructions=_counter_delta(
            counts.snapshot_constructions,
            snapshot_before,
        ),
    )


def _frame_dtypes(dataset: Any) -> dict[str, dict[str, str]]:
    frames = {
        "phospho": dataset._phospho,
        "site_metadata": dataset._site_metadata,
        "sample_metadata": dataset._sample_metadata,
        "total": dataset._total,
        "comparisons": dataset._comparisons,
    }
    return {
        name: {str(column): str(dtype) for column, dtype in frame.dtypes.items()}
        for name, frame in frames.items()
        if frame is not None
    }


def run_benchmark(config: BenchmarkConfig | None = None) -> BenchmarkResult:
    from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow

    resolved_config = config or default_config()
    dataset = _build_dataset(resolved_config)
    differential_request = _differential_request(dataset)
    kinase_request = _kinase_request(dataset, resolved_config)

    with _instrument_dataset(dataset) as counts:
        differential_first = _measure_run(
            lambda: DifferentialAnalysisWorkflow().run(differential_request),
            counts,
        )
        differential_repeated = _measure_run(
            lambda: DifferentialAnalysisWorkflow().run(differential_request),
            counts,
        )
        kinase_first = _measure_run(
            lambda: KinaseWorkflow().run(kinase_request),
            counts,
        )
        kinase_repeated = _measure_run(
            lambda: KinaseWorkflow().run(kinase_request),
            counts,
        )

    return BenchmarkResult(
        config=resolved_config,
        dataset_dimensions={
            "phospho_rows": int(dataset._phospho.shape[0]),
            "phospho_columns": int(dataset._phospho.shape[1]),
            "site_metadata_rows": int(dataset._site_metadata.shape[0]),
            "site_metadata_columns": int(dataset._site_metadata.shape[1]),
        },
        frame_dtypes=_frame_dtypes(dataset),
        differential_first=differential_first,
        differential_repeated=differential_repeated,
        kinase_first=kinase_first,
        kinase_repeated=kinase_repeated,
        total_full_frame_deep_copies=dict(
            sorted(counts.full_frame_deep_copies.items())
        ),
        total_snapshot_constructions=dict(
            sorted(counts.snapshot_constructions.items())
        ),
        environment=_environment(),
        dependencies=_dependencies(),
    )


def _measurement_payload(measurement: RunMeasurement) -> dict[str, Any]:
    return {
        "runtime_seconds": measurement.runtime_seconds,
        "peak_tracemalloc_mib": measurement.peak_tracemalloc_mib,
        "full_frame_deep_copies": measurement.full_frame_deep_copies,
        "snapshot_constructions": measurement.snapshot_constructions,
    }


def report_payload(result: BenchmarkResult) -> dict[str, Any]:
    return {
        "benchmark": "repeated_workflow_dataset_snapshot_reuse",
        "config": {
            "n_sites": int(result.config.n_sites),
            "n_samples": int(result.config.n_samples),
            "n_kinases": int(result.config.n_kinases),
            "substrates_per_kinase": int(result.config.substrates_per_kinase),
        },
        "dataset": {
            "dimensions": result.dataset_dimensions,
            "frame_dtypes": result.frame_dtypes,
        },
        "runs": {
            "differential": {
                "first": _measurement_payload(result.differential_first),
                "repeated": _measurement_payload(result.differential_repeated),
            },
            "kinase": {
                "first": _measurement_payload(result.kinase_first),
                "repeated": _measurement_payload(result.kinase_repeated),
            },
        },
        "totals": {
            "full_frame_deep_copies": result.total_full_frame_deep_copies,
            "snapshot_constructions": result.total_snapshot_constructions,
        },
        "environment": result.environment,
        "dependencies": result.dependencies,
        "observation_scope": (
            "Machine-specific local benchmark observation; not a portable "
            "runtime or memory guarantee."
        ),
    }


def _print_metrics(result: BenchmarkResult) -> None:
    payload = report_payload(result)
    print("benchmark_suite=repeated_workflow_dataset_snapshot_reuse_v1")
    print(f"dataset_phospho_rows={result.dataset_dimensions['phospho_rows']}")
    print(f"dataset_phospho_columns={result.dataset_dimensions['phospho_columns']}")
    print(
        "dataset_site_metadata_columns="
        f"{result.dataset_dimensions['site_metadata_columns']}"
    )
    print(f"frame_dtypes_json={_json_metric(result.frame_dtypes)}")
    print(
        "differential_first_run_seconds="
        f"{result.differential_first.runtime_seconds:.6f}"
    )
    print(
        "differential_repeated_run_seconds="
        f"{result.differential_repeated.runtime_seconds:.6f}"
    )
    print(
        "differential_first_peak_tracemalloc_mib="
        f"{result.differential_first.peak_tracemalloc_mib:.3f}"
    )
    print(
        "differential_repeated_peak_tracemalloc_mib="
        f"{result.differential_repeated.peak_tracemalloc_mib:.3f}"
    )
    print(f"kinase_first_run_seconds={result.kinase_first.runtime_seconds:.6f}")
    print(f"kinase_repeated_run_seconds={result.kinase_repeated.runtime_seconds:.6f}")
    print(
        "kinase_first_peak_tracemalloc_mib="
        f"{result.kinase_first.peak_tracemalloc_mib:.3f}"
    )
    print(
        "kinase_repeated_peak_tracemalloc_mib="
        f"{result.kinase_repeated.peak_tracemalloc_mib:.3f}"
    )
    print(
        "full_frame_deep_copy_counts_json="
        f"{_json_metric(result.total_full_frame_deep_copies)}"
    )
    print(
        "snapshot_construction_counts_json="
        f"{_json_metric(result.total_snapshot_constructions)}"
    )
    print(f"environment_json={_json_metric(payload['environment'])}")
    print(f"dependency_versions_json={_json_metric(payload['dependencies'])}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Measure repeated differential and kinase workflow reuse of the same "
            "dataset-owned immutable snapshots."
        )
    )
    parser.add_argument("--n-sites", type=int, default=default_config().n_sites)
    parser.add_argument("--n-samples", type=int, default=default_config().n_samples)
    parser.add_argument("--n-kinases", type=int, default=default_config().n_kinases)
    parser.add_argument(
        "--substrates-per-kinase",
        type=int,
        default=default_config().substrates_per_kinase,
    )
    parser.add_argument(
        "--write-report",
        action="store_true",
        help="Write a JSON report under benchmarks/reports/.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    result = run_benchmark(
        BenchmarkConfig(
            n_sites=int(args.n_sites),
            n_samples=int(args.n_samples),
            n_kinases=int(args.n_kinases),
            substrates_per_kinase=int(args.substrates_per_kinase),
        )
    )
    _print_metrics(result)
    if bool(args.write_report):
        report_path = default_report_path()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            json.dumps(report_payload(result), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"report_path={report_path.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
