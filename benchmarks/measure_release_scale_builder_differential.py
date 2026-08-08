"""Measure the optional 50,000-site x 48-sample builder+differential workload.

This is an explicitly invoked local benchmark. It is intentionally excluded
from pytest, CI, and release gates because runtime and process memory are
machine-dependent observations rather than portable release-blocking facts.
"""
# ruff: noqa: E402

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
for _path in (ROOT, SRC):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
)
from phospy.advanced import (
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetSiteMatrixConfig,
    DifferentialAnalysisConfig,
    EmpiricalBayesConfig,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.contracts.configs.differential import (
    IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
)
from phospy.provenance.hashing import (
    DEFAULT_TOLERANCE_TABLE_HASH_ALGORITHM,
    hash_json_payload,
    hash_table_tolerance,
)
from phospy.provenance.models import (
    JsonValue,
    PreprocessingStageProvenance,
    TableFingerprint,
)
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from tests.support.performance_contracts import (
    DEFAULT_PERFORMANCE_SEED,
    deterministic_analysis_ready_site_keys,
    deterministic_analysis_ready_site_metadata,
    deterministic_matrix,
    deterministic_sample_columns,
    measure_wall_clock,
    with_missing_fraction,
)

REPORTS_DIRECTORY = Path("benchmarks/reports")
REPORT_FILENAME = "release-scale-builder-differential.json"
RSS_UNAVAILABLE = "unavailable"
CONTRAST_NAME = "C2_vs_C1"
EXPECTED_PREPROCESSING_STAGE_SEQUENCE = (
    "localisation_confidence",
    "missing_data",
    "intensity_transform",
    "normalisation",
)
WORKLOAD_SEGMENT_KEYS = (
    "request_preparation_seconds",
    "builder_execution_seconds",
    "differential_analysis_seconds",
    "result_table_assembly_seconds",
)
DOCUMENTED_MAIN_METRIC_KEYS = (
    "sites",
    "samples",
    "missing_fraction",
    "original_missing_cell_count",
    "final_missing_cell_count",
    "output_rows",
    "output_columns",
    "tested_feature_count",
    "total_runtime_seconds",
    "request_preparation_seconds",
    "builder_execution_seconds",
    "preprocessing_execution_seconds",
    "preprocessing_report_assembly_seconds",
    "provenance_fingerprinting_seconds",
    "differential_analysis_seconds",
    "result_table_assembly_seconds",
    "result_table_fingerprinting_seconds",
    "scientific_summary_assembly_seconds",
    "process_rss_peak_mib",
    "scientific_summary_digest",
    "report_path",
)


@dataclass(frozen=True)
class ReleaseScaleBenchmarkConfig:
    n_sites: int = 50_000
    n_samples: int = 48
    missing_fraction: float = 0.03
    seed: int = DEFAULT_PERFORMANCE_SEED
    site_start: int = 700_000
    gene_prefix: str = "RELGENE"
    sequence_width: int = 31
    shifted_row_fraction: float = 0.08
    condition_shift: float = 2.5


@dataclass
class _TimingRecorder:
    timings: dict[str, float] = field(default_factory=dict)

    def add(self, key: str, seconds: float) -> None:
        self.timings[key] = self.timings.get(key, 0.0) + float(seconds)

    def measure(self, key: str, func: Callable[[], object]) -> object:
        result, elapsed_seconds = measure_wall_clock(func, warmup=False)
        self.add(key, elapsed_seconds)
        return result


@dataclass(frozen=True)
class _ReleaseScaleWorkflowMeasurement:
    dataset: AnalysisReadyPhosphoDataset
    result: DifferentialAnalysisResult
    table: pd.DataFrame
    timings: Mapping[str, float]
    original_missing_cell_count: int


@dataclass(frozen=True)
class _ReleaseScaleScientificSummary:
    canonical_json: str

    @classmethod
    def from_payload(
        cls, payload: Mapping[str, object]
    ) -> _ReleaseScaleScientificSummary:
        json_ready = _json_compatible_payload(payload)
        if not isinstance(json_ready, Mapping):
            raise AssertionError(
                "release-scale scientific summary must be a JSON object"
            )
        return cls(_canonical_json(json_ready))

    def to_payload(self) -> dict[str, object]:
        loaded = json.loads(self.canonical_json)
        if not isinstance(loaded, dict):
            raise AssertionError(
                "release-scale scientific summary JSON is not an object"
            )
        return loaded

    @property
    def digest(self) -> str:
        return hash_json_payload(cast(JsonValue, self.to_payload()))


@dataclass(frozen=True)
class _ScientificSummaryMeasurement:
    summary: _ReleaseScaleScientificSummary
    timings: Mapping[str, float]


@dataclass(frozen=True)
class ReleaseScaleBenchmarkResult:
    config: ReleaseScaleBenchmarkConfig
    timings: Mapping[str, float]
    metrics: Mapping[str, object]
    scientific_summary: Mapping[str, object]


def default_config() -> ReleaseScaleBenchmarkConfig:
    """Return the declared default release-scale benchmark configuration."""

    return ReleaseScaleBenchmarkConfig()


def run_benchmark(
    config: ReleaseScaleBenchmarkConfig | None = None,
) -> ReleaseScaleBenchmarkResult:
    resolved_config = default_config() if config is None else config
    started = time.perf_counter()

    workflow = _run_release_scale_workflow(resolved_config)
    _assert_release_scale_outputs(
        config=resolved_config,
        dataset=workflow.dataset,
        result=workflow.result,
        table=workflow.table,
    )

    summary_started = time.perf_counter()
    summary_measurement = _build_release_scale_scientific_summary(
        config=resolved_config,
        dataset=workflow.dataset,
        result=workflow.result,
        table=workflow.table,
        original_missing_cell_count=workflow.original_missing_cell_count,
    )
    summary_seconds = time.perf_counter() - summary_started

    timings = dict(workflow.timings)
    for key, value in summary_measurement.timings.items():
        timings[key] = float(value)
    timings["scientific_summary_assembly_seconds"] = float(summary_seconds)
    timings["single_workload_runtime_seconds"] = sum(
        timings[key] for key in WORKLOAD_SEGMENT_KEYS
    )
    timings["total_runtime_seconds"] = float(time.perf_counter() - started)

    summary_payload = summary_measurement.summary.to_payload()
    row_count, column_count = _summary_dimensions(summary_payload)
    metrics: dict[str, object] = {
        "sites": int(resolved_config.n_sites),
        "samples": int(resolved_config.n_samples),
        "missing_fraction": float(resolved_config.missing_fraction),
        "original_missing_cell_count": _required_int(
            summary_payload,
            "original_missing_cell_count",
        ),
        "final_missing_cell_count": _required_int(
            summary_payload,
            "final_missing_cell_count",
        ),
        "output_rows": row_count,
        "output_columns": column_count,
        "tested_feature_count": _required_int(
            summary_payload,
            "tested_feature_count",
        ),
        "process_rss_peak_mib": _rss_peak_report_value(_process_rss_peak_mib()),
        "scientific_summary_digest": summary_measurement.summary.digest,
    }
    metrics.update({key: float(value) for key, value in timings.items()})

    return ReleaseScaleBenchmarkResult(
        config=resolved_config,
        timings=timings,
        metrics=metrics,
        scientific_summary=summary_payload,
    )


def write_report(result: ReleaseScaleBenchmarkResult) -> Path:
    report_path = default_report_path()
    report_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "benchmark": "release_scale_builder_differential",
        "description": (
            "Optional local 50,000-site x 48-sample builder, preprocessing, "
            "provenance/fingerprinting, and one-contrast differential benchmark."
        ),
        "config": asdict(result.config),
        "metrics": dict(result.metrics),
        "timings": dict(result.timings),
        "scientific_summary": dict(result.scientific_summary),
        "python": {
            "version": sys.version,
            "executable": sys.executable,
        },
    }
    report_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return report_path


def default_report_path() -> Path:
    return Path(__file__).resolve().parents[1] / REPORTS_DIRECTORY / REPORT_FILENAME


def _run_release_scale_workflow(
    config: ReleaseScaleBenchmarkConfig,
) -> _ReleaseScaleWorkflowMeasurement:
    recorder = _TimingRecorder()
    request = recorder.measure(
        "request_preparation_seconds",
        lambda: _build_release_scale_dataset_request(config),
    )
    if not isinstance(request, DatasetBuildRequest):
        raise AssertionError("release-scale request preparation returned wrong type")
    request = cast(DatasetBuildRequest, request)
    if not isinstance(request.phospho, pd.DataFrame):
        raise AssertionError("release-scale request phospho table is not a DataFrame")
    original_missing_cell_count = int(request.phospho.isna().sum().sum())

    with _record_builder_phase_timings(recorder):
        dataset = recorder.measure(
            "builder_execution_seconds",
            lambda: AnalysisReadyDatasetBuilder().run(request),
        )
    if not isinstance(dataset, AnalysisReadyPhosphoDataset):
        raise AssertionError("release-scale builder returned wrong type")

    differential_request = _differential_request(dataset)
    result = recorder.measure(
        "differential_analysis_seconds",
        lambda: DifferentialAnalysisWorkflow().run(differential_request),
    )
    if not isinstance(result, DifferentialAnalysisResult):
        raise AssertionError("release-scale differential returned wrong type")

    table = recorder.measure(
        "result_table_assembly_seconds",
        lambda: result.table_for(CONTRAST_NAME),
    )
    if not isinstance(table, pd.DataFrame):
        raise AssertionError("release-scale result table export returned wrong type")

    return _ReleaseScaleWorkflowMeasurement(
        dataset=dataset,
        result=result,
        table=table,
        timings=dict(recorder.timings),
        original_missing_cell_count=original_missing_cell_count,
    )


@contextmanager
def _record_builder_phase_timings(recorder: _TimingRecorder) -> Iterator[None]:
    originals = {
        "_run_preprocessing": DatasetBuildExecutor._run_preprocessing,
        "_assemble_preprocessing_report": (
            DatasetBuildExecutor._assemble_preprocessing_report
        ),
        "_assemble_run_provenance": DatasetBuildExecutor._assemble_run_provenance,
    }

    def timed_method(
        timing_key: str,
        method: Callable[..., object],
    ) -> Callable[..., object]:
        def wrapper(self: object, *args: object, **kwargs: object) -> object:
            started = time.perf_counter()
            try:
                return method(self, *args, **kwargs)
            finally:
                recorder.add(timing_key, time.perf_counter() - started)

        return wrapper

    DatasetBuildExecutor._run_preprocessing = cast(
        Any,
        timed_method(
            "preprocessing_execution_seconds",
            originals["_run_preprocessing"],
        ),
    )
    DatasetBuildExecutor._assemble_preprocessing_report = cast(
        Any,
        timed_method(
            "preprocessing_report_assembly_seconds",
            originals["_assemble_preprocessing_report"],
        ),
    )
    DatasetBuildExecutor._assemble_run_provenance = cast(
        Any,
        timed_method(
            "provenance_fingerprinting_seconds",
            originals["_assemble_run_provenance"],
        ),
    )
    try:
        yield
    finally:
        for name, method in originals.items():
            setattr(DatasetBuildExecutor, name, method)


def _build_release_scale_scientific_summary(
    *,
    config: ReleaseScaleBenchmarkConfig,
    dataset: AnalysisReadyPhosphoDataset,
    result: DifferentialAnalysisResult,
    table: pd.DataFrame,
    original_missing_cell_count: int,
) -> _ScientificSummaryMeasurement:
    if dataset.provenance is None:
        raise AssertionError("release-scale dataset provenance is missing")
    if result.policy_provenance is None:
        raise AssertionError("release-scale differential policy provenance is missing")
    if result.workflow_provenance is None:
        raise AssertionError(
            "release-scale differential workflow provenance is missing"
        )

    stage_sequence = tuple(
        stage.stage for stage in dataset.provenance.preprocessing_stages
    )
    if stage_sequence != EXPECTED_PREPROCESSING_STAGE_SEQUENCE:
        raise AssertionError(
            "release-scale preprocessing stages changed: "
            f"expected={EXPECTED_PREPROCESSING_STAGE_SEQUENCE!r}; "
            f"observed={stage_sequence!r}"
        )

    expected_original_missing_cell_count = int(
        round(
            int(config.n_sites) * int(config.n_samples) * float(config.missing_fraction)
        )
    )
    if int(original_missing_cell_count) != expected_original_missing_cell_count:
        raise AssertionError(
            "release-scale original missing-cell count changed: "
            f"expected={expected_original_missing_cell_count}; "
            f"observed={original_missing_cell_count}"
        )

    result_fingerprint_started = time.perf_counter()
    result_table_fingerprint = _compact_dataframe_fingerprint_payload(
        table,
        name=f"differential.{CONTRAST_NAME}.result_table",
    )
    result_table_fingerprinting_seconds = (
        time.perf_counter() - result_fingerprint_started
    )

    final_missing_cell_count = int(dataset.phospho.isna().sum().sum())
    input_table_fingerprints = [
        _compact_table_fingerprint_payload(fingerprint)
        for fingerprint in dataset.provenance.input_tables
    ]
    output_table_fingerprints = [
        _compact_table_fingerprint_payload(fingerprint)
        for fingerprint in dataset.provenance.output_tables
    ]
    preprocessing_trace = [
        _compact_preprocessing_stage_payload(stage)
        for stage in dataset.provenance.preprocessing_stages
    ]
    policy_payload = asdict(result.policy_provenance)
    workflow_payload = dict(result.workflow_provenance)
    dataset_workflow_payload = {
        "workflow_name": dataset.provenance.workflow_name,
        "workflow_parameters": dict(dataset.provenance.workflow_parameters),
        "random_state": dataset.provenance.random_state,
        "random_seed_policy": dataset.provenance.random_seed_policy,
        "scientific_policies": [
            policy.to_payload() for policy in dataset.provenance.scientific_policies
        ],
    }
    payload = {
        "matrix_dimensions": {
            "rows": int(dataset.phospho.shape[0]),
            "columns": int(dataset.phospho.shape[1]),
        },
        "original_missing_cell_count": int(original_missing_cell_count),
        "expected_original_missing_cell_count": expected_original_missing_cell_count,
        "final_missing_cell_count": final_missing_cell_count,
        "expected_preprocessing_stage_sequence": list(
            EXPECTED_PREPROCESSING_STAGE_SEQUENCE
        ),
        "preprocessing_stage_sequence": list(stage_sequence),
        "input_table_fingerprints": input_table_fingerprints,
        "output_table_fingerprints": output_table_fingerprints,
        "input_table_fingerprints_digest": _stable_payload_digest(
            input_table_fingerprints
        ),
        "output_table_fingerprints_digest": _stable_payload_digest(
            output_table_fingerprints
        ),
        "processing_state_completeness_flags": _processing_state_completeness_flags(
            dataset
        ),
        "preprocessing_trace": preprocessing_trace,
        "preprocessing_trace_digest": _stable_payload_digest(preprocessing_trace),
        "contrast_names": [CONTRAST_NAME],
        "tested_feature_count": int(
            (table["result_status"].astype(str) == "tested").sum()
        ),
        "differential_result_table_fingerprint": result_table_fingerprint,
        "differential_result_table_digest": result_table_fingerprint["hash_value"],
        "differential_policy_provenance_digest": _stable_payload_digest(policy_payload),
        "differential_workflow_provenance_digest": _stable_payload_digest(
            workflow_payload
        ),
        "dataset_workflow_provenance_digest": _stable_payload_digest(
            dataset_workflow_payload
        ),
        "provenance_complete": (
            dataset.provenance is not None
            and result.policy_provenance is not None
            and result.workflow_provenance is not None
        ),
    }
    return _ScientificSummaryMeasurement(
        summary=_ReleaseScaleScientificSummary.from_payload(payload),
        timings={
            "result_table_fingerprinting_seconds": (
                result_table_fingerprinting_seconds
            ),
        },
    )


def _compact_table_fingerprint_payload(
    fingerprint: TableFingerprint,
) -> dict[str, object]:
    return {
        "name": fingerprint.name,
        "rows": int(fingerprint.rows),
        "columns": int(fingerprint.columns),
        "tolerance_hash_algorithm": fingerprint.tolerance_hash_algorithm,
        "tolerance_hash_value": fingerprint.tolerance_hash_value,
    }


def _compact_dataframe_fingerprint_payload(
    table: pd.DataFrame,
    *,
    name: str,
) -> dict[str, object]:
    return {
        "name": name,
        "rows": int(table.shape[0]),
        "columns": int(table.shape[1]),
        "hash_algorithm": DEFAULT_TOLERANCE_TABLE_HASH_ALGORITHM,
        "hash_value": hash_table_tolerance(table, name=name),
    }


def _compact_preprocessing_stage_payload(
    stage: PreprocessingStageProvenance,
) -> dict[str, object]:
    return {
        "stage": str(stage.stage),
        "operation": str(stage.operation),
        "input_shape": list(stage.input_shape),
        "output_shape": list(stage.output_shape),
        "input_hash": stage.input_hash,
        "output_hash": stage.output_hash,
        "phospho_input_hash": stage.phospho_input_hash,
        "phospho_output_hash": stage.phospho_output_hash,
        "parameters_digest": _stable_payload_digest(stage.parameters),
        "dropped_row_count": int(stage.dropped_row_count),
        "imputed_cell_count": int(stage.imputed_cell_count),
    }


def _processing_state_completeness_flags(
    dataset: AnalysisReadyPhosphoDataset,
) -> dict[str, object]:
    state = dataset.processing_state
    missing_value_count = state.missing_data.missing_value_count
    return {
        "missing_data_complete_matrix": bool(state.missing_data.complete_matrix),
        "missing_data_imputed": bool(state.missing_data.imputed),
        "missing_data_has_missing_values": bool(state.missing_data.has_missing_values),
        "missing_data_missing_value_count": (
            0 if missing_value_count is None else int(missing_value_count)
        ),
        "ruv_enabled": bool(state.ruv_readiness.enabled),
        "ruv_ready": bool(state.ruv_readiness.ready),
        "ruv_requires_complete_matrix": bool(
            state.ruv_readiness.requires_complete_matrix
        ),
        "ruv_matrix_complete": bool(state.ruv_readiness.matrix_complete),
        "ruv_missingness_mask_preserved": bool(
            state.ruv_readiness.missingness_mask_preserved
        ),
    }


def _stable_payload_digest(payload: object) -> str:
    return hash_json_payload(cast(JsonValue, _json_compatible_payload(payload)))


def _json_compatible_payload(value: object) -> object:
    if is_dataclass(value) and not isinstance(value, type):
        return _json_compatible_payload(asdict(value))
    if isinstance(value, Mapping):
        return {
            str(key): _json_compatible_payload(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, tuple | list):
        return [_json_compatible_payload(item) for item in value]
    if isinstance(value, np.generic):
        return _json_compatible_payload(value.item())
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, Enum):
        enum_value = value.value
        if isinstance(enum_value, str | int | float | bool):
            return enum_value
    return value


def _canonical_json(payload: Mapping[str, object]) -> str:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _required_int(payload: Mapping[str, object], key: str) -> int:
    value = payload[key]
    if isinstance(value, bool) or not isinstance(value, int):
        raise AssertionError(f"release-scale summary field {key!r} must be an int")
    return int(value)


def _summary_dimensions(payload: Mapping[str, object]) -> tuple[int, int]:
    raw_dimensions = payload["matrix_dimensions"]
    if not isinstance(raw_dimensions, Mapping):
        raise AssertionError("release-scale matrix_dimensions must be an object")
    dimensions = cast(Mapping[str, object], raw_dimensions)
    return (
        _required_int(dimensions, "rows"),
        _required_int(dimensions, "columns"),
    )


def _assert_release_scale_outputs(
    *,
    config: ReleaseScaleBenchmarkConfig,
    dataset: AnalysisReadyPhosphoDataset,
    result: DifferentialAnalysisResult,
    table: pd.DataFrame,
) -> None:
    expected_shape = (int(config.n_sites), int(config.n_samples))
    if dataset.phospho.shape != expected_shape:
        raise AssertionError(
            f"release-scale final matrix shape changed: {dataset.phospho.shape!r}"
        )
    if dataset.sample_metadata is None:
        raise AssertionError("release-scale sample metadata is missing")
    if dataset.sample_metadata.shape[0] != int(config.n_samples):
        raise AssertionError("release-scale sample metadata row count is wrong")
    if dataset.site_metadata.shape[0] != int(config.n_sites):
        raise AssertionError("release-scale site metadata row count is wrong")
    if "site_sequence" not in dataset.site_metadata.columns:
        raise AssertionError("release-scale site metadata lost required site_sequence")
    if int(dataset.phospho.isna().sum().sum()) != 0:
        raise AssertionError("release-scale final matrix still contains missing values")
    if dataset.preprocessing_report is None:
        raise AssertionError("release-scale preprocessing report is missing")
    if dataset.preprocessing_report.row_counts.empty:
        raise AssertionError("release-scale preprocessing row counts are missing")
    if dataset.processing_state.missing_data.complete_matrix is not True:
        raise AssertionError("release-scale missing-data state is not complete")
    if dataset.processing_state.ruv_readiness.missingness_mask_preserved is not True:
        raise AssertionError("release-scale missingness mask was not preserved")
    if dataset.provenance is None:
        raise AssertionError("release-scale dataset provenance is missing")
    stage_names = {stage.stage for stage in dataset.provenance.preprocessing_stages}
    if not set(EXPECTED_PREPROCESSING_STAGE_SEQUENCE).issubset(stage_names):
        raise AssertionError("release-scale preprocessing provenance is incomplete")
    if not dataset.provenance.input_tables:
        raise AssertionError("release-scale input table fingerprints are missing")
    if not dataset.provenance.output_tables:
        raise AssertionError("release-scale output table fingerprints are missing")

    if set(result.contrast_tables) != {CONTRAST_NAME}:
        raise AssertionError("release-scale differential contrast output changed")
    if table.shape[0] != int(config.n_sites):
        raise AssertionError("release-scale differential result row count is wrong")
    if table.index.name != "site_key":
        raise AssertionError("release-scale differential result index is not site_key")
    if (
        table.loc[:, "site_key"].astype(str).tolist()
        != table.index.astype(str).tolist()
    ):
        raise AssertionError("release-scale result site_key column and index differ")
    required_columns = {
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val",
        "imputed_fraction",
        "result_status",
        "result_status_reason",
    }
    if not required_columns.issubset(table.columns):
        raise AssertionError("release-scale result table columns are incomplete")
    if not (table["result_status"].astype(str) == "tested").all():
        raise AssertionError("release-scale differential did not test every feature")
    if not table["imputed_fraction"].between(0.0, 0.25, inclusive="both").all():
        raise AssertionError("release-scale imputed_fraction output is out of bounds")
    if float(table["imputed_fraction"].max()) <= 0.0:
        raise AssertionError("release-scale result did not retain missingness evidence")
    if not np.isfinite(table["logFC"].to_numpy(dtype=float)).all():
        raise AssertionError("release-scale logFC contains non-finite values")
    if result.policy_provenance is None:
        raise AssertionError("release-scale differential policy provenance is missing")
    if result.workflow_provenance is None:
        raise AssertionError(
            "release-scale differential workflow provenance is missing"
        )
    if result.input_dataset_preprocessing_report is None:
        raise AssertionError("release-scale differential did not retain preprocessing")
    if result.empirical_bayes_trend is not True:
        raise AssertionError("release-scale differential trend policy changed")
    if result.mean_variance_trend_diagnostics is None:
        raise AssertionError("release-scale trend diagnostics are missing")
    if result.prior_diagnostics.prior_variance.shape != (int(config.n_sites),):
        raise AssertionError("release-scale prior diagnostics feature count is wrong")


def _differential_request(
    dataset: AnalysisReadyPhosphoDataset,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=_two_condition_design(sample_ids=dataset.phospho.columns),
        contrasts=(
            Contrast(
                name=CONTRAST_NAME,
                numerator_condition="C2",
                denominator_condition="C1",
            ),
        ),
        config=DifferentialAnalysisConfig(
            imputed_value_policy=IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
            imputed_value_max_fraction=0.25,
            empirical_bayes=EmpiricalBayesConfig(
                method="standard",
                trend=True,
            ),
        ),
    )


def _build_release_scale_dataset_request(
    config: ReleaseScaleBenchmarkConfig,
) -> DatasetBuildRequest:
    sample_columns = deterministic_sample_columns(
        int(config.n_samples),
        prefix="release_sample",
    )
    site_keys = deterministic_analysis_ready_site_keys(
        int(config.n_sites),
        start=int(config.site_start),
        gene_prefix=config.gene_prefix,
    )
    phospho = deterministic_matrix(
        n_sites=int(config.n_sites),
        n_samples=int(config.n_samples),
        seed=int(config.seed) + 50_000,
        site_ids=site_keys,
        sample_columns=sample_columns,
    )
    phospho = phospho + 40.0
    shifted_rows = int(int(config.n_sites) * float(config.shifted_row_fraction))
    phospho.iloc[:shifted_rows, int(config.n_samples) // 2 :] += float(
        config.condition_shift
    )
    phospho = with_missing_fraction(
        phospho,
        missing_fraction=float(config.missing_fraction),
        seed=int(config.seed) + 50_001,
    )
    return DatasetBuildRequest(
        phospho=phospho,
        site_metadata=_release_scale_site_metadata(site_keys, config=config),
        sample_metadata=_release_scale_sample_metadata(sample_columns),
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig(
            intensity_transform=DatasetIntensityTransformConfig(policy="log2"),
            normalisation=DatasetNormalisationConfig(policy="median_center"),
            missing_data=DatasetMissingDataConfig(
                policy="impute_row_median",
                min_observed_values=1,
            ),
            site_matrix=DatasetSiteMatrixConfig(policy="as_input"),
        ),
    )


def _release_scale_site_metadata(
    site_keys: pd.Index,
    *,
    config: ReleaseScaleBenchmarkConfig,
) -> pd.DataFrame:
    metadata = deterministic_analysis_ready_site_metadata(
        site_keys,
        start=int(config.site_start),
        gene_prefix=config.gene_prefix,
        sequence_width=int(config.sequence_width),
    )
    row_count = metadata.shape[0]
    row_positions = np.arange(row_count, dtype=int)
    metadata = metadata.assign(
        protein_accession=[
            f"UPI{position:09d}" for position in range(1, row_count + 1)
        ],
        isoform_label=[
            f"{config.gene_prefix}{position:05d}-{(position % 3) + 1}"
            for position in range(1, row_count + 1)
        ],
        evidence_count=(2 + (row_positions % 5)).astype(int),
        peptide_count=(1 + (row_positions % 4)).astype(int),
    )
    return metadata


def _release_scale_sample_metadata(sample_ids: pd.Index) -> pd.DataFrame:
    count = int(sample_ids.size)
    condition = np.asarray(["C1"] * (count // 2) + ["C2"] * (count // 2), dtype=object)
    batch = np.asarray([f"batch_{(index % 6) + 1}" for index in range(count)])
    return pd.DataFrame(
        {
            "sample_id": sample_ids.astype(str).tolist(),
            "condition": condition.tolist(),
            "batch": batch.tolist(),
            "instrument": [f"orbitrap_{(index % 4) + 1}" for index in range(count)],
            "donor_id": [f"donor_{(index % 24) + 1:02d}" for index in range(count)],
            "acquisition_order": list(range(1, count + 1)),
            "injection_volume_ul": [1.0 + (index % 3) * 0.1 for index in range(count)],
        },
        index=sample_ids.copy(),
    )


def _two_condition_design(*, sample_ids: Iterable[str]) -> ExperimentalDesign:
    sample_list = [str(sample_id) for sample_id in sample_ids]
    midpoint = len(sample_list) // 2
    records: list[SampleDesignRecord] = []
    for sample_index, sample_id in enumerate(sample_list):
        condition = "C1" if sample_index < midpoint else "C2"
        replicate_index = (
            sample_index + 1 if condition == "C1" else sample_index - midpoint + 1
        )
        records.append(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=condition,
                biological_replicate_id=f"{condition}_rep{replicate_index:02d}",
            )
        )
    return ExperimentalDesign(samples=tuple(records))


def _process_rss_peak_mib() -> float | None:
    if os.name == "posix":
        return _posix_process_rss_peak_mib()
    if os.name == "nt":
        return _windows_process_rss_peak_mib()
    return None


def _posix_process_rss_peak_mib() -> float | None:
    status_path = Path(f"/proc/{os.getpid()}/status")
    try:
        for line in status_path.read_text(encoding="utf-8").splitlines():
            if line.startswith("VmHWM:"):
                parts = line.split()
                if len(parts) >= 2:
                    return float(int(parts[1]) * 1024) / (1024.0 * 1024.0)
    except OSError:
        return None
    return None


def _windows_process_rss_peak_mib() -> float | None:
    try:
        import ctypes
        from ctypes import wintypes
    except ImportError:
        return None

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

    process_query_limited_information = 0x1000
    process_vm_read = 0x0010
    handle = ctypes.windll.kernel32.OpenProcess(  # type: ignore[attr-defined]
        process_query_limited_information | process_vm_read,
        False,
        os.getpid(),
    )
    if not handle:
        return None
    try:
        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(ProcessMemoryCounters)
        ok = ctypes.windll.psapi.GetProcessMemoryInfo(  # type: ignore[attr-defined]
            handle,
            ctypes.byref(counters),
            counters.cb,
        )
        if not ok:
            return None
        return float(int(counters.PeakWorkingSetSize)) / (1024.0 * 1024.0)
    finally:
        ctypes.windll.kernel32.CloseHandle(handle)  # type: ignore[attr-defined]


def _rss_peak_report_value(value: float | None) -> float | str:
    return RSS_UNAVAILABLE if value is None else round(float(value), 6)


def _format_metric_value(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def format_metric_values(metrics: Mapping[str, object]) -> str:
    return " ".join(
        f"{key}={_format_metric_value(metrics[key])}"
        for key in DOCUMENTED_MAIN_METRIC_KEYS
        if key in metrics
    )


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the optional local 50,000-site x 48-sample PhosPy "
            "builder+differential benchmark once."
        ),
    )
    parser.add_argument(
        "--write-report",
        action="store_true",
        help=(
            "Write a JSON report under benchmarks/reports/. The benchmark always "
            "prints key=value metrics to stdout."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    result = run_benchmark()
    metrics = dict(result.metrics)
    metrics["report_path"] = "not_written"
    if args.write_report:
        metrics["report_path"] = str(write_report(result))
    formatted_metrics = format_metric_values(metrics)
    print(f"benchmark=release_scale_builder_differential {formatted_metrics}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
