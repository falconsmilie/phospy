from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from phospy.errors.validation import (
    DatasetValidationError,
    TransformationValidationError,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import Organism
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from phospy.validation.datasets.analysis_ready import (
    AnalysisReadyDatasetModelBoundaryValidator,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)

_BOUNDARY_VALIDATOR = AnalysisReadyDatasetModelBoundaryValidator()
_MODEL_BOUNDARY_ERRORS = (DatasetValidationError, TransformationValidationError)
_CENTRED_Y_SEQUENCE = "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"
_CENTRED_T_SEQUENCE = "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA"
_CENTRED_S_SEQUENCE = "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA"


def _site_key(*, protein_identifier: str, residue: str, position: int) -> str:
    key = build_protein_scoped_site_key(
        organism="rat",
        protein_namespace="protein_id",
        protein_identifier=protein_identifier,
        residue=residue,
        position=position,
        field_name="test.site_key",
        error_type=ValueError,
    )
    return encode_site_key(key)


def _valid_payload() -> dict[str, object]:
    site_keys = [
        _site_key(protein_identifier="MAPK14", residue="Y", position=182),
        _site_key(protein_identifier="AKT1", residue="T", position=308),
    ]
    site_index = pd.Index(site_keys, name="site_key")
    return {
        "phospho": pd.DataFrame(
            {
                "sample_a": [1.0, 2.0],
                "sample_b": [1.5, 2.5],
            },
            index=site_index.copy(),
        ),
        "site_metadata": pd.DataFrame(
            {
                "site_key": site_keys,
                "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
                "organism": ["rat", "rat"],
                "protein_namespace": ["protein_id", "protein_id"],
                "protein_identifier": ["MAPK14", "AKT1"],
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": [
                    _CENTRED_Y_SEQUENCE,
                    _CENTRED_T_SEQUENCE,
                ],
                "protein_id": ["MAPK14", "AKT1"],
            },
            index=site_index.copy(),
        ),
        "sample_metadata": pd.DataFrame(
            {"condition": ["a", "b"]},
            index=pd.Index(["sample_a", "sample_b"], name="sample_id"),
        ),
        "total": None,
        "comparisons": None,
        "organism": Organism.RAT,
        "intensity_scale_state": supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        "processing_state": supported_linear_processing_state(has_total_matrix=False),
        "preprocessing_report": None,
        "provenance": None,
    }


def _assert_constructor_and_adapter_reject(payload: dict[str, object]) -> None:
    with pytest.raises(_MODEL_BOUNDARY_ERRORS):
        AnalysisReadyPhosphoDataset(**payload)
    with pytest.raises(_MODEL_BOUNDARY_ERRORS):
        _BOUNDARY_VALIDATOR.run(**payload)


def test_model_boundary_validator_accepts_valid_payload() -> None:
    payload = _valid_payload()
    constructed = AnalysisReadyPhosphoDataset(**payload)
    validated = _BOUNDARY_VALIDATOR.run(**payload)
    assert isinstance(constructed, AnalysisReadyPhosphoDataset)
    assert isinstance(validated, AnalysisReadyPhosphoDataset)
    assert constructed.phospho.index.name == "site_key"
    assert constructed.site_metadata.index.name == "site_key"
    assert constructed.site_metadata.loc[:, "site_key"].tolist() == (
        constructed.site_metadata.index.tolist()
    )


@pytest.mark.parametrize(
    ("protein_id_values", "expected_present"),
    [
        pytest.param(None, False, id="absent"),
        pytest.param(["MAPK14", "AKT1"], True, id="complete"),
        pytest.param(["MAPK14", pd.NA], True, id="partially-missing"),
        pytest.param(["", "  "], True, id="blank"),
    ],
)
def test_model_boundary_validator_treats_protein_id_as_optional_signalome_metadata(
    protein_id_values: list[object] | None,
    expected_present: bool,
) -> None:
    payload = _valid_payload()
    site_metadata = payload["site_metadata"].copy(deep=True)
    if protein_id_values is None:
        site_metadata = site_metadata.drop(columns=["protein_id"])
    else:
        site_metadata.loc[:, "protein_id"] = protein_id_values
    payload["site_metadata"] = site_metadata

    constructed = AnalysisReadyPhosphoDataset(**payload)
    validated = _BOUNDARY_VALIDATOR.run(**payload)

    assert ("protein_id" in constructed.site_metadata.columns) is expected_present
    assert ("protein_id" in validated.site_metadata.columns) is expected_present


def test_model_boundary_validator_rejects_display_indexed_direct_constructor_payload() -> (
    None
):
    payload = _valid_payload()
    phospho = payload["phospho"].copy(deep=True)
    site_metadata = payload["site_metadata"].copy(deep=True)
    phospho.index = pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id")
    site_metadata.index = phospho.index.copy()
    payload["phospho"] = phospho
    payload["site_metadata"] = site_metadata

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="display-indexed"):
        AnalysisReadyPhosphoDataset(**payload)
    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="display-indexed"):
        _BOUNDARY_VALIDATOR.run(**payload)


def test_model_boundary_validator_rejects_missing_site_key_column() -> None:
    payload = _valid_payload()
    payload["site_metadata"] = payload["site_metadata"].drop(columns=["site_key"])
    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="site_key"):
        AnalysisReadyPhosphoDataset(**payload)


def test_model_boundary_validator_rejects_missing_display_id_column() -> None:
    payload = _valid_payload()
    payload["site_metadata"] = payload["site_metadata"].drop(columns=["display_id"])
    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="display_id"):
        AnalysisReadyPhosphoDataset(**payload)


@pytest.mark.parametrize(
    "column_name",
    ["organism", "protein_namespace", "protein_identifier"],
)
def test_model_boundary_validator_rejects_missing_site_key_context_column(
    column_name: str,
) -> None:
    payload = _valid_payload()
    payload["site_metadata"] = payload["site_metadata"].drop(columns=[column_name])

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match=column_name):
        AnalysisReadyPhosphoDataset(**payload)


def test_model_boundary_validator_rejects_site_key_residue_mismatch() -> None:
    payload = _valid_payload()
    site_metadata = payload["site_metadata"].copy(deep=True)
    row_key = site_metadata.index[1]
    site_metadata.loc[row_key, "display_id"] = "AKT1;S308;"
    site_metadata.loc[row_key, "site"] = "S308"
    site_metadata.loc[row_key, "site_sequence"] = _CENTRED_S_SEQUENCE
    payload["site_metadata"] = site_metadata

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="metadata-derived"):
        AnalysisReadyPhosphoDataset(**payload)


def test_model_boundary_validator_rejects_site_key_position_mismatch() -> None:
    payload = _valid_payload()
    site_metadata = payload["site_metadata"].copy(deep=True)
    row_key = site_metadata.index[1]
    site_metadata.loc[row_key, "display_id"] = "AKT1;T309;"
    site_metadata.loc[row_key, "site"] = "T309"
    payload["site_metadata"] = site_metadata

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="metadata-derived"):
        AnalysisReadyPhosphoDataset(**payload)


def test_model_boundary_validator_rejects_site_key_protein_identifier_mismatch() -> (
    None
):
    payload = _valid_payload()
    site_metadata = payload["site_metadata"].copy(deep=True)
    row_key = site_metadata.index[1]
    site_metadata.loc[row_key, "protein_identifier"] = "P31749"
    payload["site_metadata"] = site_metadata

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="metadata-derived"):
        AnalysisReadyPhosphoDataset(**payload)


def test_model_boundary_validator_allows_duplicate_display_id_with_distinct_site_key() -> (
    None
):
    first_key = _site_key(protein_identifier="P28482", residue="Y", position=182)
    second_key = _site_key(protein_identifier="Q9WVS8", residue="Y", position=182)
    site_index = pd.Index([first_key, second_key], name="site_key")
    payload = _valid_payload()
    payload["phospho"] = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=site_index.copy(),
    )
    payload["site_metadata"] = pd.DataFrame(
        {
            "site_key": [first_key, second_key],
            "display_id": ["MAPK14;Y182;", "MAPK14;Y182;"],
            "organism": ["rat", "rat"],
            "protein_namespace": ["protein_id", "protein_id"],
            "protein_identifier": ["P28482", "Q9WVS8"],
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": [_CENTRED_Y_SEQUENCE, _CENTRED_Y_SEQUENCE],
            "protein_id": ["P28482", "Q9WVS8"],
        },
        index=site_index.copy(),
    )

    dataset = AnalysisReadyPhosphoDataset(**payload)

    assert dataset.site_metadata.loc[:, "display_id"].nunique() == 1
    assert dataset.site_metadata.loc[:, "site_key"].nunique() == 2


def test_model_boundary_validator_rejects_duplicate_site_key() -> None:
    payload = _valid_payload()
    duplicate_key = payload["phospho"].index[0]  # type: ignore[index]
    duplicate_index = pd.Index([duplicate_key, duplicate_key], name="site_key")
    phospho = payload["phospho"].copy(deep=True)
    site_metadata = payload["site_metadata"].copy(deep=True)
    phospho.index = duplicate_index.copy()
    site_metadata.index = duplicate_index.copy()
    site_metadata.loc[:, "site_key"] = duplicate_index.tolist()
    payload["phospho"] = phospho
    payload["site_metadata"] = site_metadata

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="duplicate_site_key"):
        AnalysisReadyPhosphoDataset(**payload)


def test_model_boundary_validator_rejects_site_key_column_index_mismatch() -> None:
    payload = _valid_payload()
    site_metadata = payload["site_metadata"].copy(deep=True)
    site_metadata.iloc[1, site_metadata.columns.get_loc("site_key")] = _site_key(
        protein_identifier="AKT1",
        residue="T",
        position=309,
    )
    payload["site_metadata"] = site_metadata

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="site_key.*index"):
        AnalysisReadyPhosphoDataset(**payload)


def test_model_boundary_validator_rejects_phospho_site_metadata_index_mismatch() -> (
    None
):
    payload = _valid_payload()
    phospho = payload["phospho"].copy(deep=True)
    phospho.index = pd.Index(list(reversed(phospho.index.tolist())), name="site_key")
    payload["phospho"] = phospho

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="must exactly match"):
        AnalysisReadyPhosphoDataset(**payload)


def test_analysis_ready_from_owned_rejects_invalid_site_identity() -> None:
    payload = _valid_payload()
    phospho = payload["phospho"].copy(deep=True)
    site_metadata = payload["site_metadata"].copy(deep=True)
    phospho.index = pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id")
    site_metadata.index = phospho.index.copy()
    payload["phospho"] = phospho
    payload["site_metadata"] = site_metadata

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="display-indexed"):
        AnalysisReadyPhosphoDataset._from_owned(**payload)


def test_model_boundary_validator_parity_for_misaligned_site_metadata() -> None:
    payload = _valid_payload()
    mismatched_site_key = _site_key(
        protein_identifier="AKT1",
        residue="T",
        position=309,
    )
    payload["site_metadata"] = pd.DataFrame(
        {
            "site_key": [
                payload["site_metadata"].index[0],  # type: ignore[index]
                mismatched_site_key,
            ],
            "display_id": ["MAPK14;Y182;", "AKT1;T308;"],
            "organism": ["rat", "rat"],
            "protein_namespace": ["protein_id", "protein_id"],
            "protein_identifier": ["MAPK14", "AKT1"],
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                _CENTRED_Y_SEQUENCE,
                _CENTRED_T_SEQUENCE,
            ],
        },
        index=pd.Index(
            [
                payload["site_metadata"].index[0],  # type: ignore[index]
                mismatched_site_key,
            ],
            name="site_key",
        ),
    )
    _assert_constructor_and_adapter_reject(payload)


def test_model_boundary_validator_parity_for_missing_sample_metadata_rows() -> None:
    payload = _valid_payload()
    payload["sample_metadata"] = pd.DataFrame(
        {"condition": ["a"]},
        index=pd.Index(["sample_a"], name="sample_id"),
    )
    _assert_constructor_and_adapter_reject(payload)


def test_model_boundary_validator_parity_for_missing_or_invalid_scale_state() -> None:
    payload = _valid_payload()
    payload["intensity_scale_state"] = None
    _assert_constructor_and_adapter_reject(payload)


def test_model_boundary_validator_parity_for_invalid_processing_state() -> None:
    payload = _valid_payload()
    processing_state = payload["processing_state"]
    assert hasattr(processing_state, "missing_data")
    invalid_missing_data = replace(processing_state.missing_data, complete_matrix=False)
    payload["processing_state"] = replace(
        processing_state, missing_data=invalid_missing_data
    )
    _assert_constructor_and_adapter_reject(payload)


def test_model_boundary_validator_parity_for_missing_site_sequence() -> None:
    payload = _valid_payload()
    payload["site_metadata"] = payload["site_metadata"].drop(columns=["site_sequence"])
    _assert_constructor_and_adapter_reject(payload)


def test_model_boundary_validator_parity_for_incomplete_phospho_matrix() -> None:
    payload = _valid_payload()
    phospho = payload["phospho"].copy(deep=True)
    phospho.iloc[0, 0] = float("nan")
    payload["phospho"] = phospho
    _assert_constructor_and_adapter_reject(payload)


@pytest.mark.parametrize(
    ("field_name", "value"),
    [
        pytest.param("provenance", "invalid", id="invalid-provenance-type"),
        pytest.param(
            "preprocessing_report", "invalid", id="invalid-preprocessing-report-type"
        ),
    ],
)
def test_model_boundary_validator_parity_for_invalid_metadata_types(
    field_name: str,
    value: object,
) -> None:
    payload = _valid_payload()
    payload[field_name] = value
    _assert_constructor_and_adapter_reject(payload)
