#!/usr/bin/env python3
"""Benchmark repeated workflow reuse of dataset-owned immutable snapshots.

Targets:
- `phospy.workflows.differential.DifferentialAnalysisWorkflow.run`
- `phospy.workflows.kinase.KinaseWorkflow.run`

This is an explicitly invoked local benchmark. It records machine-specific
runtime, tracemalloc peak memory, full-frame and projected pandas deep-copy
counts, and dataset internal snapshot construction counts for setup and
repeated workflow use of the same unchanged dataset.
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
from dataclasses import dataclass, field
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
    "setup_run_seconds",
    "setup_peak_tracemalloc_mib",
    "differential_first_run_seconds",
    "differential_repeated_run_seconds",
    "differential_first_peak_tracemalloc_mib",
    "differential_repeated_peak_tracemalloc_mib",
    "kinase_first_run_seconds",
    "kinase_repeated_run_seconds",
    "kinase_first_peak_tracemalloc_mib",
    "kinase_repeated_peak_tracemalloc_mib",
    "setup_full_frame_deep_copy_counts_json",
    "per_run_full_frame_deep_copy_counts_json",
    "full_frame_deep_copy_counts_json",
    "projected_frame_deep_copy_counts_json",
    "per_run_projected_frame_deep_copy_counts_json",
    "setup_snapshot_construction_counts_json",
    "per_run_snapshot_construction_counts_json",
    "snapshot_construction_counts_json",
    "environment_json",
    "dependency_versions_json",
)
_FRAME_SOURCE_ATTR = "_phospy_repeated_workflow_frame_source"


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
    projected_frame_deep_copies: dict[str, int]
    snapshot_constructions: dict[str, int]


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    config: BenchmarkConfig
    dataset_dimensions: dict[str, int]
    frame_dtypes: dict[str, dict[str, str]]
    setup: RunMeasurement
    differential_first: RunMeasurement
    differential_repeated: RunMeasurement
    kinase_first: RunMeasurement
    kinase_repeated: RunMeasurement
    total_full_frame_deep_copies: dict[str, int]
    total_projected_frame_deep_copies: dict[str, int]
    total_snapshot_constructions: dict[str, int]
    environment: dict[str, str]
    dependencies: dict[str, str]


@dataclass(frozen=True, slots=True)
class _DatasetFrameSource:
    label: str
    index: tuple[object, ...]
    columns: tuple[object, ...]


@dataclass(slots=True)
class _InstrumentationCounts:
    full_frame_deep_copies: Counter[str]
    projected_frame_deep_copies: Counter[str]
    snapshot_constructions: Counter[str]
    frame_sources_by_id: dict[int, _DatasetFrameSource] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class _BenchmarkSetup:
    dataset: Any
    differential_request: Any
    kinase_request: Any


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


def _mark_frame_source(
    counts: _InstrumentationCounts,
    frame: pd.DataFrame,
    *,
    label: str,
) -> None:
    source = _DatasetFrameSource(
        label=label,
        index=tuple(frame.index.tolist()),
        columns=tuple(frame.columns.tolist()),
    )
    counts.frame_sources_by_id[id(frame)] = source
    frame.attrs[_FRAME_SOURCE_ATTR] = source


def _mark_dataset_frame_sources(
    counts: _InstrumentationCounts,
    dataset: Any,
) -> None:
    for label, frame in (
        ("dataset.phospho", dataset._phospho),
        ("dataset.site_metadata", dataset._site_metadata),
        ("dataset.sample_metadata", dataset._sample_metadata),
        ("dataset.total", dataset._total),
        ("dataset.comparisons", dataset._comparisons),
    ):
        if isinstance(frame, pd.DataFrame):
            _mark_frame_source(counts, frame, label=label)


def _frame_source_for(
    counts: _InstrumentationCounts,
    frame: pd.DataFrame,
) -> _DatasetFrameSource | None:
    source = counts.frame_sources_by_id.get(id(frame))
    if source is not None:
        return source
    attr_source = frame.attrs.get(_FRAME_SOURCE_ATTR)
    if isinstance(attr_source, _DatasetFrameSource):
        return attr_source
    return None


def _is_full_source_frame(
    frame: pd.DataFrame,
    source: _DatasetFrameSource,
) -> bool:
    return (
        tuple(frame.index.tolist()) == source.index
        and tuple(frame.columns.tolist()) == source.columns
    )


def _snapshot_source_label(field_name: str) -> str | None:
    suffix = " internal snapshot"
    if not field_name.endswith(suffix):
        return None
    return field_name[: -len(suffix)]


def _mark_snapshot_frame_source(
    counts: _InstrumentationCounts,
    snapshot: object,
    *,
    field_name: str,
    source_value: pd.DataFrame,
) -> None:
    snapshot_frame = getattr(snapshot, "_frame", None)
    if not isinstance(snapshot_frame, pd.DataFrame):
        return
    source = _frame_source_for(counts, source_value)
    label = source.label if source is not None else _snapshot_source_label(field_name)
    if label is None:
        return
    _mark_frame_source(counts, snapshot_frame, label=label)


def _build_dataset(
    config: BenchmarkConfig,
    *,
    counts: _InstrumentationCounts | None = None,
):
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
    sample_metadata = pd.DataFrame(
        {
            "condition": [
                "A" if sample_id.startswith("A_") else "B" for sample_id in sample_ids
            ],
            "batch": [
                f"batch_{(position % 2) + 1}"
                for position, _sample_id in enumerate(sample_ids)
            ],
        },
        index=pd.Index(sample_ids, name="sample_id"),
    )
    comparisons = pd.DataFrame(
        {
            "B_vs_A_expected_shift": np.full(
                int(config.n_sites),
                0.35,
                dtype=np.float64,
            )
        },
        index=site_index.copy(),
    )
    if counts is not None:
        _mark_frame_source(counts, phospho, label="dataset.phospho")
        _mark_frame_source(counts, site_metadata, label="dataset.site_metadata")
        _mark_frame_source(counts, sample_metadata, label="dataset.sample_metadata")
        _mark_frame_source(counts, comparisons, label="dataset.comparisons")
    dataset = trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=sample_metadata,
        comparisons=comparisons,
        organism=Organism.RAT,
        intensity_scale_state=supported_log2_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_log2_processing_state(has_total_matrix=False),
    )
    if counts is not None:
        _mark_dataset_frame_sources(counts, dataset)
    return dataset


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


@contextmanager
def _instrument_copy_accounting() -> Iterator[_InstrumentationCounts]:
    import phospy.science.datasets.internal_frame_store as internal_frame_store_module

    counts = _InstrumentationCounts(
        full_frame_deep_copies=Counter(),
        projected_frame_deep_copies=Counter(),
        snapshot_constructions=Counter(),
    )
    original_copy = pd.DataFrame.copy
    original_snapshot = internal_frame_store_module.immutable_dataframe_snapshot
    original_optional_snapshot = (
        internal_frame_store_module.immutable_optional_dataframe_snapshot
    )

    def wrapped_copy(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        deep = kwargs.get("deep", args[0] if args else True)
        if bool(deep):
            source = _frame_source_for(counts, self)
            if source is not None:
                if _is_full_source_frame(self, source):
                    counts.full_frame_deep_copies.update((source.label,))
                else:
                    counts.projected_frame_deep_copies.update(
                        (f"{source.label} projected",)
                    )
        return original_copy(self, *args, **kwargs)

    def counting_snapshot(
        value: pd.DataFrame,
        *,
        field_name: str,
        error_type: type[Exception] = TypeError,
    ):
        counts.snapshot_constructions.update((field_name,))
        snapshot = original_snapshot(
            value,
            field_name=field_name,
            error_type=error_type,
        )
        _mark_snapshot_frame_source(
            counts,
            snapshot,
            field_name=field_name,
            source_value=value,
        )
        return snapshot

    def counting_optional_snapshot(
        value: pd.DataFrame | None,
        *,
        field_name: str,
        error_type: type[Exception] = TypeError,
    ):
        if isinstance(value, pd.DataFrame):
            counts.snapshot_constructions.update((field_name,))
        snapshot = original_optional_snapshot(
            value,
            field_name=field_name,
            error_type=error_type,
        )
        if isinstance(value, pd.DataFrame):
            _mark_snapshot_frame_source(
                counts,
                snapshot,
                field_name=field_name,
                source_value=value,
            )
        return snapshot

    pd.DataFrame.copy = wrapped_copy  # pyright: ignore[reportAttributeAccessIssue] - benchmark instrumentation monkeypatches pandas copy counts.
    internal_frame_store_module.immutable_dataframe_snapshot = counting_snapshot
    internal_frame_store_module.immutable_optional_dataframe_snapshot = (
        counting_optional_snapshot
    )
    try:
        yield counts
    finally:
        pd.DataFrame.copy = original_copy  # pyright: ignore[reportAttributeAccessIssue] - restore benchmark instrumentation.
        internal_frame_store_module.immutable_dataframe_snapshot = original_snapshot
        internal_frame_store_module.immutable_optional_dataframe_snapshot = (
            original_optional_snapshot
        )


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
    _result, measurement = _measure_operation(run, counts)
    return measurement


def _measure_operation(
    run: Callable[[], object],
    counts: _InstrumentationCounts,
) -> tuple[object, RunMeasurement]:
    full_copy_before = counts.full_frame_deep_copies.copy()
    projected_copy_before = counts.projected_frame_deep_copies.copy()
    snapshot_before = counts.snapshot_constructions.copy()
    tracemalloc.start()
    started = time.perf_counter()
    result = run()
    runtime_seconds = time.perf_counter() - started
    _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return (
        result,
        RunMeasurement(
            runtime_seconds=float(runtime_seconds),
            peak_tracemalloc_mib=float(peak_bytes / (1024.0 * 1024.0)),
            full_frame_deep_copies=_counter_delta(
                counts.full_frame_deep_copies,
                full_copy_before,
            ),
            projected_frame_deep_copies=_counter_delta(
                counts.projected_frame_deep_copies,
                projected_copy_before,
            ),
            snapshot_constructions=_counter_delta(
                counts.snapshot_constructions,
                snapshot_before,
            ),
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


def _build_setup(
    config: BenchmarkConfig,
    counts: _InstrumentationCounts,
) -> _BenchmarkSetup:
    dataset = _build_dataset(config, counts=counts)
    differential_request = _differential_request(dataset)
    kinase_request = _kinase_request(dataset, config)
    _mark_dataset_frame_sources(counts, dataset)
    return _BenchmarkSetup(
        dataset=dataset,
        differential_request=differential_request,
        kinase_request=kinase_request,
    )


def run_benchmark(config: BenchmarkConfig | None = None) -> BenchmarkResult:
    from phospy import DifferentialAnalysisWorkflow, KinaseWorkflow

    resolved_config = config or default_config()
    with _instrument_copy_accounting() as counts:
        setup_payload, setup = _measure_operation(
            lambda: _build_setup(resolved_config, counts),
            counts,
        )
        if not isinstance(setup_payload, _BenchmarkSetup):
            raise TypeError("benchmark setup did not return _BenchmarkSetup")
        dataset = setup_payload.dataset
        differential_request = setup_payload.differential_request
        kinase_request = setup_payload.kinase_request
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
        setup=setup,
        differential_first=differential_first,
        differential_repeated=differential_repeated,
        kinase_first=kinase_first,
        kinase_repeated=kinase_repeated,
        total_full_frame_deep_copies=dict(
            sorted(counts.full_frame_deep_copies.items())
        ),
        total_projected_frame_deep_copies=dict(
            sorted(counts.projected_frame_deep_copies.items())
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
        "projected_frame_deep_copies": measurement.projected_frame_deep_copies,
        "snapshot_constructions": measurement.snapshot_constructions,
    }


def _per_phase_metric(
    result: BenchmarkResult,
    attribute_name: str,
) -> dict[str, dict[str, int]]:
    return {
        "setup": getattr(result.setup, attribute_name),
        "differential_first": getattr(result.differential_first, attribute_name),
        "differential_repeated": getattr(
            result.differential_repeated,
            attribute_name,
        ),
        "kinase_first": getattr(result.kinase_first, attribute_name),
        "kinase_repeated": getattr(result.kinase_repeated, attribute_name),
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
        "setup": _measurement_payload(result.setup),
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
            "projected_frame_deep_copies": result.total_projected_frame_deep_copies,
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
    print(f"setup_run_seconds={result.setup.runtime_seconds:.6f}")
    print(f"setup_peak_tracemalloc_mib={result.setup.peak_tracemalloc_mib:.3f}")
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
        "setup_full_frame_deep_copy_counts_json="
        f"{_json_metric(result.setup.full_frame_deep_copies)}"
    )
    print(
        "per_run_full_frame_deep_copy_counts_json="
        f"{_json_metric(_per_phase_metric(result, 'full_frame_deep_copies'))}"
    )
    print(
        "full_frame_deep_copy_counts_json="
        f"{_json_metric(result.total_full_frame_deep_copies)}"
    )
    print(
        "projected_frame_deep_copy_counts_json="
        f"{_json_metric(result.total_projected_frame_deep_copies)}"
    )
    print(
        "per_run_projected_frame_deep_copy_counts_json="
        f"{_json_metric(_per_phase_metric(result, 'projected_frame_deep_copies'))}"
    )
    print(
        "setup_snapshot_construction_counts_json="
        f"{_json_metric(result.setup.snapshot_constructions)}"
    )
    print(
        "per_run_snapshot_construction_counts_json="
        f"{_json_metric(_per_phase_metric(result, 'snapshot_constructions'))}"
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
