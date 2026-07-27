from __future__ import annotations

import json
import os
import sys
import time
import tracemalloc
from collections.abc import Callable, Iterable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    AnalysisReadyPhosphoDataset,
    DifferentialAnalysisWorkflow,
)
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetMissingDataConfig,
    DatasetNormalisationConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DifferentialAnalysisConfig,
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
    EmpiricalBayesConfig,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from phospy.contracts.configs.differential import (
    IMPUTED_VALUE_POLICY_WITHHOLD_IMPUTED_FEATURES,
)
from phospy.science.datasets.builders.executor import DatasetBuildExecutor
from tests.support.performance_contracts import (
    DEFAULT_PERFORMANCE_SEED,
    END_TO_END_RELEASE_SCALE_INSTRUMENTED_TIMEOUT_SECONDS,
    END_TO_END_RELEASE_SCALE_MISSING_FRACTION,
    END_TO_END_RELEASE_SCALE_N_SAMPLES,
    END_TO_END_RELEASE_SCALE_N_SITES,
    END_TO_END_RELEASE_SCALE_PEAK_MIB_MAX,
    END_TO_END_RELEASE_SCALE_RUNTIME_SECONDS_MAX,
    deterministic_analysis_ready_site_keys,
    deterministic_analysis_ready_site_metadata,
    deterministic_matrix,
    deterministic_sample_columns,
    measure_wall_clock,
    run_subprocess_with_peak_rss,
    with_missing_fraction,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]

_ROOT = Path(__file__).resolve().parents[2]
_MEMORY_PROBE_PREFIX = "release_scale_memory_probe_json="
_REPORT_FILENAME = "release-scale-performance-contract.json"
_PRODUCTION_SEGMENT_KEYS = (
    "dataset_request_preparation_seconds",
    "builder_execution_seconds",
    "differential_execution_seconds",
    "serialization_report_assembly_seconds",
)


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
class _ReleaseScaleProductionMeasurement:
    dataset: AnalysisReadyPhosphoDataset
    result: DifferentialAnalysisResult
    table: pd.DataFrame
    timings: Mapping[str, float]


@dataclass(frozen=True)
class _ReleaseScaleMemoryMeasurement:
    instrumented_runtime_seconds: float
    tracemalloc_peak_mib: float
    subprocess_elapsed_seconds: float
    process_rss_peak_mib: float | None
    metrics: Mapping[str, object]


def test_end_to_end_release_scale_builder_and_differential_contract(
    record_property: Any,
) -> None:
    production = _run_release_scale_production_measurement()

    _assert_release_scale_outputs(
        dataset=production.dataset,
        result=production.result,
        table=production.table,
    )
    assert (
        production.timings["production_runtime_seconds"]
        < END_TO_END_RELEASE_SCALE_RUNTIME_SECONDS_MAX
    )

    memory = _run_release_scale_memory_probe_subprocess()

    _record_release_scale_properties(
        record_property=record_property,
        production=production,
        memory=memory,
    )
    _write_release_scale_report(production=production, memory=memory)
    _print_release_scale_summary(production=production, memory=memory)

    assert memory.tracemalloc_peak_mib < END_TO_END_RELEASE_SCALE_PEAK_MIB_MAX


def _run_release_scale_production_measurement() -> _ReleaseScaleProductionMeasurement:
    recorder = _TimingRecorder()
    request = recorder.measure(
        "dataset_request_preparation_seconds",
        _build_release_scale_dataset_request,
    )
    if not isinstance(request, DatasetBuildRequest):
        raise AssertionError("release-scale request preparation returned wrong type")

    with _record_builder_phase_timings(recorder):
        dataset = recorder.measure(
            "builder_execution_seconds",
            lambda: AnalysisReadyDatasetBuilder().run(request),
        )
    if not isinstance(dataset, AnalysisReadyPhosphoDataset):
        raise AssertionError("release-scale builder returned wrong type")

    differential_request = _differential_request(dataset)
    result = recorder.measure(
        "differential_execution_seconds",
        lambda: DifferentialAnalysisWorkflow().run(differential_request),
    )
    if not isinstance(result, DifferentialAnalysisResult):
        raise AssertionError("release-scale differential returned wrong type")

    table = recorder.measure(
        "serialization_report_assembly_seconds",
        lambda: result.table_for("C2_vs_C1"),
    )
    if not isinstance(table, pd.DataFrame):
        raise AssertionError("release-scale result table export returned wrong type")

    recorder.timings["production_runtime_seconds"] = sum(
        recorder.timings[key] for key in _PRODUCTION_SEGMENT_KEYS
    )
    return _ReleaseScaleProductionMeasurement(
        dataset=dataset,
        result=result,
        table=table,
        timings=dict(recorder.timings),
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

    DatasetBuildExecutor._run_preprocessing = timed_method(
        "preprocessing_execution_seconds",
        originals["_run_preprocessing"],
    )
    DatasetBuildExecutor._assemble_preprocessing_report = timed_method(
        "preprocessing_report_assembly_seconds",
        originals["_assemble_preprocessing_report"],
    )
    DatasetBuildExecutor._assemble_run_provenance = timed_method(
        "provenance_fingerprinting_seconds",
        originals["_assemble_run_provenance"],
    )
    try:
        yield
    finally:
        for name, method in originals.items():
            setattr(DatasetBuildExecutor, name, method)


def _run_release_scale_memory_probe_subprocess() -> _ReleaseScaleMemoryMeasurement:
    timeout_seconds = _instrumented_timeout_seconds()
    helper = (
        "from tests.performance.test_end_to_end_release_scale_contract "
        "import _release_scale_memory_probe_main; "
        "_release_scale_memory_probe_main()"
    )
    env = dict(os.environ)
    pythonpath_parts = [str(_ROOT), str(_ROOT / "src")]
    existing_pythonpath = env.get("PYTHONPATH")
    if existing_pythonpath:
        pythonpath_parts.append(existing_pythonpath)
    env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

    try:
        completed = run_subprocess_with_peak_rss(
            [sys.executable, "-c", helper],
            cwd=_ROOT,
            env=env,
            timeout_seconds=timeout_seconds,
        )
    except TimeoutError as exc:
        pytest.fail(
            "release-scale tracemalloc memory probe exceeded the explicit "
            f"{timeout_seconds:.1f}s completion timeout: {exc}"
        )

    if completed.returncode != 0:
        pytest.fail(
            "release-scale tracemalloc memory probe failed with return code "
            f"{completed.returncode}\nstdout:\n{completed.stdout}\nstderr:\n"
            f"{completed.stderr}"
        )

    metrics = _parse_memory_probe_metrics(completed.stdout)
    return _ReleaseScaleMemoryMeasurement(
        instrumented_runtime_seconds=float(metrics["instrumented_runtime_seconds"]),
        tracemalloc_peak_mib=float(metrics["tracemalloc_peak_mib"]),
        subprocess_elapsed_seconds=float(completed.elapsed_seconds),
        process_rss_peak_mib=completed.peak_rss_mib,
        metrics=metrics,
    )


def _release_scale_memory_probe_main() -> None:
    metrics = _run_release_scale_tracemalloc_measurement()
    print(
        f"{_MEMORY_PROBE_PREFIX}{json.dumps(metrics, sort_keys=True)}",
        flush=True,
    )


def _run_release_scale_tracemalloc_measurement() -> dict[str, object]:
    tracemalloc.start()
    started = time.perf_counter()
    try:
        production = _run_release_scale_production_measurement()
        instrumented_runtime_seconds = time.perf_counter() - started
        _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    _assert_release_scale_outputs(
        dataset=production.dataset,
        result=production.result,
        table=production.table,
    )
    return {
        "instrumentation": "tracemalloc",
        "instrumented_runtime_seconds": float(instrumented_runtime_seconds),
        "tracemalloc_peak_mib": float(peak_bytes) / (1024.0 * 1024.0),
        "timings": {
            str(key): float(value) for key, value in production.timings.items()
        },
        "final_matrix_shape": (
            f"{production.dataset.phospho.shape[0]}x"
            f"{production.dataset.phospho.shape[1]}"
        ),
        "tested_feature_count": int(
            (production.table["result_status"].astype(str) == "tested").sum()
        ),
    }


def _parse_memory_probe_metrics(stdout: str) -> dict[str, object]:
    for line in reversed(stdout.splitlines()):
        if line.startswith(_MEMORY_PROBE_PREFIX):
            loaded = json.loads(line.removeprefix(_MEMORY_PROBE_PREFIX))
            if not isinstance(loaded, dict):
                raise AssertionError("release-scale memory probe JSON is not an object")
            return loaded
    raise AssertionError(
        f"release-scale memory probe did not emit metrics JSON; stdout was:\n{stdout}"
    )


def _instrumented_timeout_seconds() -> float:
    raw = os.environ.get("PHOSPY_RELEASE_SCALE_INSTRUMENTED_TIMEOUT_SECONDS")
    if raw is None:
        return END_TO_END_RELEASE_SCALE_INSTRUMENTED_TIMEOUT_SECONDS
    try:
        timeout_seconds = float(raw)
    except ValueError as exc:
        raise AssertionError(
            "PHOSPY_RELEASE_SCALE_INSTRUMENTED_TIMEOUT_SECONDS must be numeric"
        ) from exc
    if timeout_seconds <= 0.0:
        raise AssertionError(
            "PHOSPY_RELEASE_SCALE_INSTRUMENTED_TIMEOUT_SECONDS must be > 0"
        )
    return timeout_seconds


def _record_release_scale_properties(
    *,
    record_property: Any,
    production: _ReleaseScaleProductionMeasurement,
    memory: _ReleaseScaleMemoryMeasurement,
) -> None:
    table = production.table
    dataset = production.dataset
    record_property("release_scale_sites", END_TO_END_RELEASE_SCALE_N_SITES)
    record_property("release_scale_samples", END_TO_END_RELEASE_SCALE_N_SAMPLES)
    record_property(
        "release_scale_missing_fraction",
        END_TO_END_RELEASE_SCALE_MISSING_FRACTION,
    )
    record_property("python_version", sys.version.split()[0])
    record_property(
        "production_runtime_seconds",
        round(production.timings["production_runtime_seconds"], 6),
    )
    record_property(
        "production_runtime_budget_seconds",
        END_TO_END_RELEASE_SCALE_RUNTIME_SECONDS_MAX,
    )
    for key, value in sorted(production.timings.items()):
        record_property(key, round(float(value), 6))
    record_property(
        "instrumented_runtime_seconds",
        round(memory.instrumented_runtime_seconds, 6),
    )
    record_property(
        "instrumented_completion_timeout_seconds",
        round(_instrumented_timeout_seconds(), 6),
    )
    record_property("tracemalloc_peak_mib", round(memory.tracemalloc_peak_mib, 6))
    record_property(
        "tracemalloc_peak_budget_mib",
        END_TO_END_RELEASE_SCALE_PEAK_MIB_MAX,
    )
    record_property(
        "process_rss_peak_mib",
        (
            "unavailable_no_portable_project_helper"
            if memory.process_rss_peak_mib is None
            else round(memory.process_rss_peak_mib, 6)
        ),
    )
    record_property(
        "final_matrix_shape",
        f"{dataset.phospho.shape[0]}x{dataset.phospho.shape[1]}",
    )
    record_property(
        "tested_feature_count",
        int((table["result_status"].astype(str) == "tested").sum()),
    )


def _write_release_scale_report(
    *,
    production: _ReleaseScaleProductionMeasurement,
    memory: _ReleaseScaleMemoryMeasurement,
) -> None:
    report_dir = Path(os.environ.get("PHOSPY_PERFORMANCE_REPORT_DIR", "build/reports"))
    report_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "contract": "end_to_end_release_scale_builder_plus_differential",
        "dimensions": {
            "sites": END_TO_END_RELEASE_SCALE_N_SITES,
            "samples": END_TO_END_RELEASE_SCALE_N_SAMPLES,
            "missing_fraction": END_TO_END_RELEASE_SCALE_MISSING_FRACTION,
        },
        "budgets": {
            "production_runtime_seconds_max": (
                END_TO_END_RELEASE_SCALE_RUNTIME_SECONDS_MAX
            ),
            "tracemalloc_peak_mib_max": END_TO_END_RELEASE_SCALE_PEAK_MIB_MAX,
            "instrumented_completion_timeout_seconds": (
                _instrumented_timeout_seconds()
            ),
        },
        "production": {
            "runtime_seconds": production.timings["production_runtime_seconds"],
            "timings": dict(production.timings),
        },
        "instrumented_memory": {
            "runtime_seconds_reported_not_budgeted": (
                memory.instrumented_runtime_seconds
            ),
            "tracemalloc_peak_mib": memory.tracemalloc_peak_mib,
            "subprocess_elapsed_seconds": memory.subprocess_elapsed_seconds,
            "process_rss_peak_mib": memory.process_rss_peak_mib,
            "metrics": dict(memory.metrics),
        },
        "assertions": {
            "final_matrix_shape": (
                f"{production.dataset.phospho.shape[0]}x"
                f"{production.dataset.phospho.shape[1]}"
            ),
            "tested_feature_count": int(
                (production.table["result_status"].astype(str) == "tested").sum()
            ),
            "provenance_complete": production.dataset.provenance is not None
            and production.result.policy_provenance is not None
            and production.result.workflow_provenance is not None,
        },
        "python": {
            "version": sys.version,
            "executable": sys.executable,
        },
    }
    (report_dir / _REPORT_FILENAME).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _print_release_scale_summary(
    *,
    production: _ReleaseScaleProductionMeasurement,
    memory: _ReleaseScaleMemoryMeasurement,
) -> None:
    rss_value = (
        "unavailable_no_portable_project_helper"
        if memory.process_rss_peak_mib is None
        else f"{memory.process_rss_peak_mib:.3f}"
    )
    segments = " ".join(
        f"{key}={production.timings[key]:.3f}"
        for key in (
            "dataset_request_preparation_seconds",
            "builder_execution_seconds",
            "preprocessing_execution_seconds",
            "preprocessing_report_assembly_seconds",
            "provenance_fingerprinting_seconds",
            "differential_execution_seconds",
            "serialization_report_assembly_seconds",
        )
        if key in production.timings
    )
    print(
        "release_scale_e2e "
        f"production_runtime_seconds="
        f"{production.timings['production_runtime_seconds']:.3f} "
        f"instrumented_runtime_seconds={memory.instrumented_runtime_seconds:.3f} "
        f"tracemalloc_peak_mib={memory.tracemalloc_peak_mib:.3f} "
        f"process_rss_peak_mib={rss_value} "
        f"final_shape={production.dataset.phospho.shape} "
        f"{segments}"
    )


def _assert_release_scale_outputs(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    result: DifferentialAnalysisResult,
    table: pd.DataFrame,
) -> None:
    assert dataset.phospho.shape == (
        END_TO_END_RELEASE_SCALE_N_SITES,
        END_TO_END_RELEASE_SCALE_N_SAMPLES,
    )
    assert dataset.sample_metadata is not None
    assert dataset.sample_metadata.shape[0] == END_TO_END_RELEASE_SCALE_N_SAMPLES
    assert dataset.site_metadata.shape[0] == END_TO_END_RELEASE_SCALE_N_SITES
    assert int(dataset.phospho.isna().sum().sum()) == 0
    assert dataset.preprocessing_report is not None
    assert not dataset.preprocessing_report.row_counts.empty
    assert dataset.processing_state.missing_data.complete_matrix is True
    assert dataset.processing_state.ruv_readiness.missingness_mask_preserved is True
    assert dataset.provenance is not None
    stage_names = {stage.stage for stage in dataset.provenance.preprocessing_stages}
    assert {
        "localisation_confidence",
        "intensity_transform",
        "normalisation",
        "missing_data",
    }.issubset(stage_names)
    assert dataset.provenance.input_tables
    assert dataset.provenance.output_tables

    assert set(result.contrast_tables) == {"C2_vs_C1"}
    assert table.shape[0] == END_TO_END_RELEASE_SCALE_N_SITES
    assert table.index.name == "site_key"
    assert (
        table.loc[:, "site_key"].astype(str).tolist()
        == table.index.astype(str).tolist()
    )
    assert {
        "logFC",
        "t",
        "P.Value",
        "adj.P.Val",
        "imputed_fraction",
        "result_status",
        "result_status_reason",
    }.issubset(table.columns)
    assert (table["result_status"].astype(str) == "tested").all()
    assert table["imputed_fraction"].between(0.0, 0.25, inclusive="both").all()
    assert float(table["imputed_fraction"].max()) > 0.0
    assert np.isfinite(table["logFC"].to_numpy(dtype=float)).all()
    assert result.policy_provenance is not None
    assert result.workflow_provenance is not None
    assert result.input_dataset_preprocessing_report is not None
    assert result.empirical_bayes_trend is True
    assert result.mean_variance_trend_diagnostics is not None
    assert result.prior_diagnostics.prior_variance.shape == (
        END_TO_END_RELEASE_SCALE_N_SITES,
    )


def _run_release_scale_workflow(
    request: DatasetBuildRequest,
) -> tuple[AnalysisReadyPhosphoDataset, DifferentialAnalysisResult]:
    dataset = AnalysisReadyDatasetBuilder().run(request)
    result = DifferentialAnalysisWorkflow().run(_differential_request(dataset))
    return dataset, result


def _differential_request(
    dataset: AnalysisReadyPhosphoDataset,
) -> DifferentialAnalysisRequest:
    return DifferentialAnalysisRequest(
        dataset=dataset,
        design=_two_condition_design(sample_ids=dataset.phospho.columns),
        contrasts=(
            Contrast(
                name="C2_vs_C1",
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


def _build_release_scale_dataset_request() -> DatasetBuildRequest:
    sample_columns = deterministic_sample_columns(
        END_TO_END_RELEASE_SCALE_N_SAMPLES,
        prefix="release_sample",
    )
    site_keys = deterministic_analysis_ready_site_keys(
        END_TO_END_RELEASE_SCALE_N_SITES,
        start=700_000,
        gene_prefix="RELGENE",
    )
    phospho = deterministic_matrix(
        n_sites=END_TO_END_RELEASE_SCALE_N_SITES,
        n_samples=END_TO_END_RELEASE_SCALE_N_SAMPLES,
        seed=DEFAULT_PERFORMANCE_SEED + 50_000,
        site_ids=site_keys,
        sample_columns=sample_columns,
    )
    phospho = phospho + 40.0
    shifted_rows = int(END_TO_END_RELEASE_SCALE_N_SITES * 0.08)
    phospho.iloc[:shifted_rows, END_TO_END_RELEASE_SCALE_N_SAMPLES // 2 :] += 2.5
    phospho = with_missing_fraction(
        phospho,
        missing_fraction=END_TO_END_RELEASE_SCALE_MISSING_FRACTION,
        seed=DEFAULT_PERFORMANCE_SEED + 50_001,
    )
    return DatasetBuildRequest(
        phospho=phospho,
        site_metadata=_release_scale_site_metadata(site_keys),
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


def _release_scale_site_metadata(site_keys: pd.Index) -> pd.DataFrame:
    metadata = deterministic_analysis_ready_site_metadata(
        site_keys,
        start=700_000,
        gene_prefix="RELGENE",
        sequence_width=31,
    )
    row_count = metadata.shape[0]
    row_positions = np.arange(row_count, dtype=int)
    metadata = metadata.assign(
        protein_accession=[
            f"UPI{position:09d}" for position in range(1, row_count + 1)
        ],
        isoform_label=[
            f"RELGENE{position:05d}-{(position % 3) + 1}"
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


if __name__ == "__main__":
    _release_scale_memory_probe_main()
