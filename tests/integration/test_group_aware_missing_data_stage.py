from __future__ import annotations

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.advanced import DatasetMissingDataConfig
from phospy.api import DatasetBuildRequest, DatasetPreprocessingConfig, Organism
from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import DeterminismKind
from phospy.science.datasets.preprocessing.stages.missing_data import (
    stage as missing_data_stage_module,
)
from tests.support.site_keys import protein_site_key_index


def _inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    display_ids = ["GENE1;S1;", "GENE2;S2;", "GENE3;S3;", "GENE4;S4;"]
    index = protein_site_key_index(
        protein_identifiers=["P1", "P2", "P3", "P4"],
        sites=["S1", "S2", "S3", "S4"],
        protein_namespace="protein_id",
        organism="rat",
    )
    columns = pd.Index(["A1", "A2", "A3", "B1", "B2", "B3"], name="sample")
    phospho = pd.DataFrame(
        [
            [1.0, 2.0, np.nan, 5.0, 6.0, 7.0],
            [np.nan, np.nan, np.nan, 8.0, 9.0, 10.0],
            [2.0, 3.0, 4.0, 6.0, 7.0, 8.0],
            [1.0, np.nan, np.nan, 4.0, 5.0, 6.0],
        ],
        index=index,
        columns=columns,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": ["GENE1", "GENE2", "GENE3", "GENE4"],
            "site": ["S1", "S2", "S3", "S4"],
            "site_sequence": ["AAAA", "CCCC", "DDDD", "EEEE"],
            "protein_id": ["P1", "P2", "P3", "P4"],
            "localisation_confidence": [0.9, 0.9, 0.9, 0.9],
        },
        index=index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"condition": ["A", "A", "A", "B", "B", "B"]},
        index=columns.copy(),
    )
    return phospho, site_metadata, sample_metadata


def _group_aware_config(seed: int) -> DatasetPreprocessingConfig:
    return DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_group_aware",
            group_column="condition",
            min_partial_observed_fraction=0.5,
            min_reference_observed_fraction=0.75,
            q=0.01,
            width=0.3,
            seed=seed,
            k=2,
            distance="nan_euclidean",
        )
    )


def _build(seed: int, *, sample_metadata: pd.DataFrame | None = None):
    phospho, site_metadata, default_sample_metadata = _inputs()
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=(
                default_sample_metadata if sample_metadata is None else sample_metadata
            ),
            organism=Organism.RAT,
            preprocessing_config=_group_aware_config(seed),
            input_intensity_scale="log2",
            allow_suspicious_declared_input_intensity_scale=True,
        )
    )


def _missing_data_stage(dataset):
    assert dataset.provenance is not None
    return next(
        stage
        for stage in dataset.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )


def test_group_aware_stage_produces_complete_mixed_mechanism_matrix(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    phospho, _, _ = _inputs()
    mechanism_inputs: list[pd.DataFrame] = []
    real_knn = missing_data_stage_module.impute_knn_targets
    real_minprob = missing_data_stage_module.impute_minprob_targets

    def _capture_knn_input(phospho: pd.DataFrame, **kwargs: object):
        mechanism_inputs.append(phospho.copy(deep=True))
        return real_knn(phospho, **kwargs)

    def _capture_minprob_input(phospho: pd.DataFrame, **kwargs: object):
        mechanism_inputs.append(phospho.copy(deep=True))
        return real_minprob(phospho, **kwargs)

    monkeypatch.setattr(
        missing_data_stage_module, "impute_knn_targets", _capture_knn_input
    )
    monkeypatch.setattr(
        missing_data_stage_module,
        "impute_minprob_targets",
        _capture_minprob_input,
    )
    dataset = _build(123)
    knn_row, minprob_row, _, dropped_row = phospho.index
    original_retained = phospho.drop(index=dropped_row)

    assert dropped_row not in dataset.phospho.index
    assert len(mechanism_inputs) == 2
    pdt.assert_frame_equal(mechanism_inputs[0], original_retained)
    pdt.assert_frame_equal(mechanism_inputs[1], original_retained)
    assert not dataset.phospho.isna().to_numpy().any()
    assert dataset.phospho.at[knn_row, "A3"] == 4.0
    assert dataset.phospho.loc[minprob_row, ["A1", "A2", "A3"]].notna().all()

    retained_input = phospho.loc[dataset.phospho.index]
    observed = retained_input.notna()
    np.testing.assert_array_equal(
        dataset.phospho.to_numpy()[observed.to_numpy()],
        retained_input.to_numpy()[observed.to_numpy()],
    )
    observation_mask = dataset.imputation_observed_mask_dataframe()
    assert observation_mask is not None
    pdt.assert_frame_equal(
        observation_mask,
        retained_input.notna(),
        check_names=False,
    )

    assert dataset.preprocessing_report is not None
    row_audit = dataset.preprocessing_report.row_audit
    assert set(row_audit["source_row_id"]) == {
        knn_row,
        minprob_row,
        dropped_row,
    }
    assert set(row_audit["action"]) == {"dropped", "imputed"}

    stage = _missing_data_stage(dataset)
    assert stage.determinism is DeterminismKind.SEEDED_STOCHASTIC
    assert stage.random_seed == 123
    assert {item.name for item in stage.consumed_input_tables} == {
        "dataset.phospho",
        "dataset.site_metadata",
        "dataset.sample_metadata",
    }
    diagnostics = stage.diagnostics or {}
    parameters = diagnostics.get("method_parameters", {})
    assert parameters["knn_target_cell_count"] == 1
    assert parameters["minprob_target_cell_count"] == 3
    assert parameters["mechanism_input"] == "original_retained_matrix"


def test_group_aware_seed_reproducibility_changes_only_minprob_cells() -> None:
    phospho, _, _ = _inputs()
    knn_row, minprob_row = phospho.index[:2]
    first = _build(123)
    repeated = _build(123)
    changed_seed = _build(456)

    pdt.assert_frame_equal(first.phospho, repeated.phospho)
    first_mask = first.imputation_observed_mask_dataframe()
    repeated_mask = repeated.imputation_observed_mask_dataframe()
    changed_mask = changed_seed.imputation_observed_mask_dataframe()
    pdt.assert_frame_equal(
        first_mask,
        repeated_mask,
    )
    pdt.assert_frame_equal(
        first_mask,
        changed_mask,
    )
    assert first.phospho.at[knn_row, "A3"] == changed_seed.phospho.at[knn_row, "A3"]
    assert not np.array_equal(
        first.phospho.loc[minprob_row, ["A1", "A2"]].to_numpy(),
        changed_seed.phospho.loc[minprob_row, ["A1", "A2"]].to_numpy(),
    )
    first_parameters = (_missing_data_stage(first).diagnostics or {})[
        "method_parameters"
    ]
    changed_parameters = (_missing_data_stage(changed_seed).diagnostics or {})[
        "method_parameters"
    ]
    assert (
        first_parameters["knn_target_mask_hash"]
        == changed_parameters["knn_target_mask_hash"]
    )
    assert (
        first_parameters["minprob_target_mask_hash"]
        == changed_parameters["minprob_target_mask_hash"]
    )


@pytest.mark.parametrize("metadata_case", ["missing", "misaligned"])
def test_group_aware_stage_requires_aligned_sample_metadata(metadata_case: str) -> None:
    phospho, site_metadata, sample_metadata = _inputs()
    resolved_metadata = None
    if metadata_case == "misaligned":
        resolved_metadata = sample_metadata.drop(index="B3")

    with pytest.raises(PhosPyInputError, match="sample_metadata"):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=resolved_metadata,
                organism=Organism.RAT,
                preprocessing_config=_group_aware_config(123),
                input_intensity_scale="log2",
                allow_suspicious_declared_input_intensity_scale=True,
            )
        )


def test_standalone_missing_data_policy_does_not_consume_sample_metadata() -> None:
    phospho, site_metadata, _ = _inputs()
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=None,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                    input_scale="log2",
                )
            ),
            input_intensity_scale="log2",
            allow_suspicious_declared_input_intensity_scale=True,
        )
    )

    stage = _missing_data_stage(dataset)
    assert "dataset.sample_metadata" not in {
        item.name for item in stage.consumed_input_tables
    }


def test_group_aware_knn_no_overlap_error_identifies_cell_context() -> None:
    phospho, site_metadata, sample_metadata = _inputs()
    phospho = phospho.iloc[:2].copy(deep=True)
    phospho.iloc[0] = [1.0, 2.0, np.nan, np.nan, np.nan, np.nan]
    phospho.iloc[1] = [np.nan, np.nan, 4.0, 5.0, 6.0, 7.0]
    site_metadata = site_metadata.loc[phospho.index]
    config = DatasetPreprocessingConfig(
        missing_data=DatasetMissingDataConfig(
            policy="impute_group_aware",
            group_column="condition",
            min_partial_observed_fraction=0.3,
            min_reference_observed_fraction=0.3,
            q=0.01,
            width=0.3,
            seed=123,
            k=1,
            distance="nan_euclidean",
        )
    )

    with pytest.raises(
        PhosPyInputError,
        match=r"no eligible donor.*affected row=.*affected column='A3'",
    ):
        AnalysisReadyDatasetBuilder().run(
            DatasetBuildRequest(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                organism=Organism.RAT,
                preprocessing_config=config,
                input_intensity_scale="log2",
                allow_suspicious_declared_input_intensity_scale=True,
            )
        )
