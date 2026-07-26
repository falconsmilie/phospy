from __future__ import annotations

from collections.abc import Iterable
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
from tests.support.performance_contracts import (
    DEFAULT_PERFORMANCE_SEED,
    END_TO_END_RELEASE_SCALE_MISSING_FRACTION,
    END_TO_END_RELEASE_SCALE_N_SAMPLES,
    END_TO_END_RELEASE_SCALE_N_SITES,
    END_TO_END_RELEASE_SCALE_PEAK_MIB_MAX,
    END_TO_END_RELEASE_SCALE_RUNTIME_SECONDS_MAX,
    deterministic_analysis_ready_site_keys,
    deterministic_analysis_ready_site_metadata,
    deterministic_matrix,
    deterministic_sample_columns,
    measure_runtime_and_peak_mib,
    with_missing_fraction,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]


def test_end_to_end_release_scale_builder_and_differential_contract(
    record_property: Any,
) -> None:
    request = _build_release_scale_dataset_request()

    dataset_and_result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        lambda: _run_release_scale_workflow(request),
        warmup=False,
    )
    dataset, result = dataset_and_result
    table = result.table_for("C2_vs_C1")
    rss_note = "unavailable_no_portable_project_helper"

    record_property("release_scale_sites", END_TO_END_RELEASE_SCALE_N_SITES)
    record_property("release_scale_samples", END_TO_END_RELEASE_SCALE_N_SAMPLES)
    record_property("runtime_seconds", round(runtime_seconds, 6))
    record_property("tracemalloc_peak_mib", round(peak_mib, 6))
    record_property("process_rss_peak_mib", rss_note)
    record_property(
        "final_matrix_shape",
        f"{dataset.phospho.shape[0]}x{dataset.phospho.shape[1]}",
    )
    record_property(
        "tested_feature_count",
        int((table["result_status"].astype(str) == "tested").sum()),
    )
    print(
        "release_scale_e2e "
        f"runtime_seconds={runtime_seconds:.3f} "
        f"tracemalloc_peak_mib={peak_mib:.3f} "
        f"process_rss_peak_mib={rss_note} "
        f"final_shape={dataset.phospho.shape}"
    )

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
    assert runtime_seconds < END_TO_END_RELEASE_SCALE_RUNTIME_SECONDS_MAX
    assert peak_mib < END_TO_END_RELEASE_SCALE_PEAK_MIB_MAX


def _run_release_scale_workflow(
    request: DatasetBuildRequest,
) -> tuple[AnalysisReadyPhosphoDataset, DifferentialAnalysisResult]:
    dataset = AnalysisReadyDatasetBuilder().run(request)
    result = DifferentialAnalysisWorkflow().run(
        DifferentialAnalysisRequest(
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
    )
    return dataset, result


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
