from __future__ import annotations

from collections.abc import Iterable

import pytest

from phospy import AnalysisReadyDatasetBuilder, DifferentialAnalysisWorkflow
from phospy.api import (
    Contrast,
    DatasetBuildRequest,
    DatasetPreprocessingConfig,
    DifferentialAnalysisRequest,
    ExperimentalDesign,
    Organism,
    SampleDesignRecord,
)
from tests.support.performance_contracts import (
    DATASET_BUILD_MEDIUM_PEAK_MIB_MAX,
    DATASET_BUILD_MEDIUM_RUNTIME_SECONDS_MAX,
    DATASET_BUILD_SMOKE_PEAK_MIB_MAX,
    DATASET_BUILD_SMOKE_RUNTIME_SECONDS_MAX,
    DIFFERENTIAL_WORKFLOW_MEDIUM_PEAK_MIB_MAX,
    DIFFERENTIAL_WORKFLOW_MEDIUM_RUNTIME_SECONDS_MAX,
    DIFFERENTIAL_WORKFLOW_SMOKE_PEAK_MIB_MAX,
    DIFFERENTIAL_WORKFLOW_SMOKE_RUNTIME_SECONDS_MAX,
    WORKFLOW_MEDIUM_CONTRACT_MISSING_FRACTION,
    WORKFLOW_MEDIUM_CONTRACT_N_CONDITIONS,
    WORKFLOW_MEDIUM_CONTRACT_N_SAMPLES,
    WORKFLOW_MEDIUM_CONTRACT_N_SITES,
    WORKFLOW_SMOKE_CONTRACT_MISSING_FRACTION,
    WORKFLOW_SMOKE_CONTRACT_N_CONDITIONS,
    WORKFLOW_SMOKE_CONTRACT_N_SAMPLES,
    WORKFLOW_SMOKE_CONTRACT_N_SITES,
    deterministic_matrix,
    deterministic_site_ids,
    deterministic_site_metadata,
    measure_runtime_and_peak_mib,
    with_missing_fraction,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]

_SCALE_CASES = (
    pytest.param(
        "smoke",
        WORKFLOW_SMOKE_CONTRACT_N_SITES,
        WORKFLOW_SMOKE_CONTRACT_N_SAMPLES,
        WORKFLOW_SMOKE_CONTRACT_N_CONDITIONS,
        WORKFLOW_SMOKE_CONTRACT_MISSING_FRACTION,
        DATASET_BUILD_SMOKE_RUNTIME_SECONDS_MAX,
        DATASET_BUILD_SMOKE_PEAK_MIB_MAX,
        DIFFERENTIAL_WORKFLOW_SMOKE_RUNTIME_SECONDS_MAX,
        DIFFERENTIAL_WORKFLOW_SMOKE_PEAK_MIB_MAX,
        id="smoke",
    ),
    pytest.param(
        "medium",
        WORKFLOW_MEDIUM_CONTRACT_N_SITES,
        WORKFLOW_MEDIUM_CONTRACT_N_SAMPLES,
        WORKFLOW_MEDIUM_CONTRACT_N_CONDITIONS,
        WORKFLOW_MEDIUM_CONTRACT_MISSING_FRACTION,
        DATASET_BUILD_MEDIUM_RUNTIME_SECONDS_MAX,
        DATASET_BUILD_MEDIUM_PEAK_MIB_MAX,
        DIFFERENTIAL_WORKFLOW_MEDIUM_RUNTIME_SECONDS_MAX,
        DIFFERENTIAL_WORKFLOW_MEDIUM_PEAK_MIB_MAX,
        id="medium",
    ),
)


def _build_dataset_request(
    *,
    n_sites: int,
    n_samples: int,
    missing_fraction: float,
    seed: int,
) -> DatasetBuildRequest:
    site_ids = deterministic_site_ids(n_sites, start=220_000, gene_prefix="PERFSITE")
    phospho = deterministic_matrix(
        n_sites=n_sites,
        n_samples=n_samples,
        seed=seed,
        site_ids=site_ids,
    )
    phospho = phospho + 30.0
    phospho_with_missing = with_missing_fraction(
        phospho,
        missing_fraction=missing_fraction,
        seed=seed + 101,
    )
    site_metadata = deterministic_site_metadata(site_ids, include_protein_id=True)
    return DatasetBuildRequest(
        phospho=phospho_with_missing,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        preprocessing_config=DatasetPreprocessingConfig.from_raw_phosphosite_table(),
    )


def _build_design(
    *,
    sample_ids: Iterable[str],
    n_conditions: int,
) -> ExperimentalDesign:
    sample_list = [str(sample_id) for sample_id in sample_ids]
    if n_conditions < 2:
        raise ValueError("n_conditions must be >= 2")
    if len(sample_list) % n_conditions != 0:
        raise ValueError(
            "sample count must be evenly divisible by n_conditions for balanced design"
        )
    replicates_per_condition = len(sample_list) // n_conditions
    records: list[SampleDesignRecord] = []
    for sample_index, sample_id in enumerate(sample_list):
        condition_index = sample_index // replicates_per_condition
        condition = f"C{condition_index + 1}"
        replicate_index = (sample_index % replicates_per_condition) + 1
        records.append(
            SampleDesignRecord(
                sample_id=sample_id,
                condition=condition,
                biological_replicate_id=f"{condition}_rep{replicate_index}",
            )
        )
    return ExperimentalDesign(samples=tuple(records))


def _build_contrasts(*, n_conditions: int) -> tuple[Contrast, ...]:
    baseline = "C1"
    return tuple(
        Contrast(
            name=f"C{condition_index}_vs_{baseline}",
            numerator_condition=f"C{condition_index}",
            denominator_condition=baseline,
        )
        for condition_index in range(2, n_conditions + 1)
    )


@pytest.mark.parametrize(
    (
        "_scale_name",
        "n_sites",
        "n_samples",
        "_n_conditions",
        "missing_fraction",
        "dataset_runtime_max",
        "dataset_peak_max",
        "_differential_runtime_max",
        "_differential_peak_max",
    ),
    _SCALE_CASES,
)
def test_dataset_builder_performance_contract_for_smoke_and_medium_scales(
    _scale_name: str,
    n_sites: int,
    n_samples: int,
    _n_conditions: int,
    missing_fraction: float,
    dataset_runtime_max: float,
    dataset_peak_max: float,
    _differential_runtime_max: float,
    _differential_peak_max: float,
) -> None:
    request = _build_dataset_request(
        n_sites=n_sites,
        n_samples=n_samples,
        missing_fraction=missing_fraction,
        seed=18_013,
    )
    builder = AnalysisReadyDatasetBuilder()
    dataset, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: builder.run(request),
        warmup=True,
    )

    assert dataset.phospho.shape == (n_sites, n_samples)
    assert int(dataset.phospho.isna().sum().sum()) == 0
    assert dataset.preprocessing_report is not None
    assert not dataset.preprocessing_report.row_counts.empty
    assert dataset.processing_state.missing_data.complete_matrix is True
    assert runtime_seconds < dataset_runtime_max
    assert peak_mib < dataset_peak_max


@pytest.mark.parametrize(
    (
        "_scale_name",
        "n_sites",
        "n_samples",
        "n_conditions",
        "missing_fraction",
        "_dataset_runtime_max",
        "_dataset_peak_max",
        "differential_runtime_max",
        "differential_peak_max",
    ),
    _SCALE_CASES,
)
def test_differential_workflow_performance_contract_for_smoke_and_medium_scales(
    _scale_name: str,
    n_sites: int,
    n_samples: int,
    n_conditions: int,
    missing_fraction: float,
    _dataset_runtime_max: float,
    _dataset_peak_max: float,
    differential_runtime_max: float,
    differential_peak_max: float,
) -> None:
    dataset = AnalysisReadyDatasetBuilder().run(
        _build_dataset_request(
            n_sites=n_sites,
            n_samples=n_samples,
            missing_fraction=missing_fraction,
            seed=19_103,
        )
    )
    design = _build_design(
        sample_ids=dataset.phospho.columns,
        n_conditions=n_conditions,
    )
    contrasts = _build_contrasts(n_conditions=n_conditions)

    result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: DifferentialAnalysisWorkflow().run(
            DifferentialAnalysisRequest(
                dataset=dataset,
                design=design,
                contrasts=contrasts,
            )
        ),
        warmup=True,
    )

    assert set(result.contrast_tables) == {contrast.name for contrast in contrasts}
    for table in result.contrast_tables.values():
        assert table.shape == (n_sites, 4)
        assert list(table.columns) == ["logFC", "t", "P.Value", "adj.P.Val"]
    assert result.policy_provenance is not None
    assert result.input_dataset_preprocessing_report is not None
    assert runtime_seconds < differential_runtime_max
    assert peak_mib < differential_peak_max
