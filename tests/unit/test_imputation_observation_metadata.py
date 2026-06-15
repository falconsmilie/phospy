from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from phospy.api import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.errors.validation import DatasetValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

_PROTEINS = ["MAPK14", "AKT1", "GSK3B"]
_SITES = ["Y182", "T308", "S9"]
_DISPLAY_IDS = ["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"]


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=_PROTEINS,
        sites=_SITES,
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": _PROTEINS,
            "site": _SITES,
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
            "protein_id": _PROTEINS,
            "localisation_confidence": [0.95, 0.95, 0.95],
            "site_key": index.astype(str).tolist(),
            "display_id": _DISPLAY_IDS,
            **site_key_context_columns(index),
        },
        index=index.copy(),
    )


def _phospho_with_imputation_cases(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, np.nan],
            "sample_b": [2.0, np.nan, 5.0],
            "sample_c": [3.0, 4.0, np.nan],
            "sample_d": [4.0, 6.0, np.nan],
        },
        index=index.copy(),
    )


def _build_imputed_dataset() -> AnalysisReadyPhosphoDataset:
    index = _site_index()
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=_phospho_with_imputation_cases(index),
            site_metadata=_site_metadata(index),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                )
            ),
        )
    )


def _complete_phospho(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0],
            "sample_b": [1.5, 2.5, 3.5],
        },
        index=index.copy(),
    )


def _construct_dataset_with_mask(
    mask: pd.DataFrame,
) -> AnalysisReadyPhosphoDataset:
    index = _site_index()
    return AnalysisReadyPhosphoDataset(
        phospho=_complete_phospho(index),
        site_metadata=_site_metadata(index),
        imputation_observation_mask=mask,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def test_imputation_metadata_records_feature_imputed_counts() -> None:
    dataset = _build_imputed_dataset()

    summary = dataset.imputation_feature_metadata

    assert summary is not None
    assert summary.index.equals(dataset.phospho.index)
    assert summary["imputed_cell_count"].tolist() == [0, 1, 3]
    assert summary["imputed_fraction"].tolist() == [0.0, 0.25, 0.75]


def test_imputation_metadata_records_observed_counts() -> None:
    dataset = _build_imputed_dataset()

    summary = dataset.imputation_feature_metadata

    assert summary is not None
    assert summary["observed_cell_count"].tolist() == [4, 3, 1]


def test_dataset_imputation_summary_preserves_requested_feature_order() -> None:
    dataset = _build_imputed_dataset()
    requested_features = [dataset.phospho.index[2], dataset.phospho.index[0]]

    summary = dataset.imputation_observation_summary_dataframe(
        feature_ids=requested_features,
        sample_ids=["sample_b", "sample_a"],
    )

    assert summary is not None
    assert summary.index.tolist() == requested_features
    assert summary["feature_id"].tolist() == requested_features
    assert summary["observed_cell_count"].tolist() == [1, 2]
    assert summary["imputed_cell_count"].tolist() == [1, 0]
    assert summary["total_analysed_cell_count"].tolist() == [2, 2]
    assert summary["imputed_fraction"].tolist() == [0.5, 0.0]


def test_dataset_imputation_summary_validates_sample_labels() -> None:
    dataset = _build_imputed_dataset()

    with pytest.raises(
        DatasetValidationError,
        match="sample_ids.*unknown_sample",
    ):
        dataset.imputation_observation_summary_dataframe(
            feature_ids=[dataset.phospho.index[0]],
            sample_ids=["sample_a", "unknown_sample"],
        )


def test_dataset_imputation_summary_validates_feature_labels() -> None:
    dataset = _build_imputed_dataset()

    with pytest.raises(
        DatasetValidationError,
        match="feature_ids.*unknown_feature",
    ):
        dataset.imputation_observation_summary_dataframe(
            feature_ids=["unknown_feature"],
            sample_ids=["sample_a"],
        )


def test_dataset_imputation_summary_fails_when_imputed_state_lacks_mask() -> None:
    index = _site_index()
    processing_state = supported_linear_processing_state(has_total_matrix=False)
    imputed_processing_state = replace(
        processing_state,
        missing_data=replace(processing_state.missing_data, imputed=True),
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=_complete_phospho(index),
        site_metadata=_site_metadata(index),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=imputed_processing_state,
    )

    with pytest.raises(
        DatasetValidationError,
        match="imputation_observation_mask.*missing_data\\.imputed",
    ):
        dataset.imputation_observation_summary_dataframe(
            feature_ids=[index[0]],
            sample_ids=["sample_a"],
        )


def test_imputation_metadata_rejects_misaligned_rows() -> None:
    index = _site_index()
    mask = pd.DataFrame(
        True,
        index=pd.Index(list(reversed(index.tolist())), name="site_key"),
        columns=pd.Index(["sample_a", "sample_b"]),
    )

    with pytest.raises(
        DatasetValidationError,
        match="dataset.imputation_observation_mask.index",
    ):
        _construct_dataset_with_mask(mask)


def test_imputation_metadata_rejects_misaligned_columns() -> None:
    index = _site_index()
    mask = pd.DataFrame(
        True,
        index=index.copy(),
        columns=pd.Index(["sample_b", "sample_a"]),
    )

    with pytest.raises(
        DatasetValidationError,
        match="dataset.imputation_observation_mask.columns",
    ):
        _construct_dataset_with_mask(mask)


def test_imputation_metadata_public_export_is_defensive() -> None:
    dataset = _build_imputed_dataset()
    summary = dataset.imputation_feature_metadata
    mask = dataset.imputation_observed_mask_dataframe()

    assert summary is not None
    assert mask is not None
    summary.iloc[0, 0] = 999
    mask.iloc[0, 0] = False

    reread_summary = dataset.imputation_feature_metadata
    reread_mask = dataset.imputation_observed_mask_dataframe()

    assert reread_summary is not None
    assert reread_mask is not None
    assert int(reread_summary.iloc[0, 0]) == 0
    assert bool(reread_mask.iloc[0, 0]) is True
