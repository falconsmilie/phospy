from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

import phospy.science.datasets.models as dataset_models
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
from phospy.science.transformations.models import (
    IntensityScaleState,
    MatrixIntensityScaleState,
)
from phospy.validation.datasets.analysis_ready import (
    AnalysisReadyDatasetModelBoundaryValidator,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
    supported_log2_processing_state,
)

_BOUNDARY_VALIDATOR = AnalysisReadyDatasetModelBoundaryValidator()
_MODEL_BOUNDARY_ERRORS = (DatasetValidationError, TransformationValidationError)
_CENTRED_Y_SEQUENCE = "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"
_CENTRED_T_SEQUENCE = "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA"
_CENTRED_S_SEQUENCE = "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA"
_NO_MISSING_PROCESSING_STATE_ERROR = (
    "dataset.phospho must not contain missing values; "
    "dataset.processing_state.missing_data claims no missing values "
    "but dataset.phospho contains missing values"
)


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


def _valid_payload_with_explicit_positions(
    *,
    position_value: object,
    site_position_value: object,
) -> dict[str, object]:
    payload = _valid_payload()
    site_metadata = payload["site_metadata"].copy(deep=True)
    site_metadata.loc[:, "position"] = pd.Series(
        [182, position_value],
        index=site_metadata.index,
        dtype="object",
    )
    site_metadata.loc[:, "site_position"] = pd.Series(
        [182, site_position_value],
        index=site_metadata.index,
        dtype="object",
    )
    payload["site_metadata"] = site_metadata
    return payload


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


def test_analysis_ready_dataset_constructor_preserves_frame_ownership() -> None:
    payload = _valid_payload()
    phospho = payload["phospho"]
    site_metadata = payload["site_metadata"]
    sample_metadata = payload["sample_metadata"]
    assert isinstance(phospho, pd.DataFrame)
    assert isinstance(site_metadata, pd.DataFrame)
    assert isinstance(sample_metadata, pd.DataFrame)

    dataset = AnalysisReadyPhosphoDataset(**payload)

    assert dataset._phospho is not phospho
    assert dataset._site_metadata is not site_metadata
    assert dataset._sample_metadata is not sample_metadata

    phospho.iloc[0, 0] = 999.0
    site_metadata.loc[site_metadata.index[0], "gene_symbol"] = "CHANGED"
    sample_metadata.loc[sample_metadata.index[0], "condition"] = "changed"
    assert float(dataset._phospho.iloc[0, 0]) == 1.0
    assert str(dataset._site_metadata.loc[site_metadata.index[0], "gene_symbol"]) == (
        "MAPK14"
    )
    assert dataset._sample_metadata is not None
    assert str(dataset._sample_metadata.loc[sample_metadata.index[0], "condition"]) == (
        "a"
    )

    exported_phospho = dataset.phospho
    exported_phospho.iloc[0, 0] = 777.0
    assert exported_phospho is not dataset._phospho
    assert float(dataset._phospho.iloc[0, 0]) == 1.0


def test_analysis_ready_dataset_constructor_requires_site_sequence() -> None:
    payload = _valid_payload()
    payload["site_metadata"] = payload["site_metadata"].drop(columns=["site_sequence"])

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="site_sequence"):
        AnalysisReadyPhosphoDataset(**payload)


def test_analysis_ready_dataset_constructor_requires_established_log2_state() -> None:
    payload = _valid_payload()
    payload["intensity_scale_state"] = IntensityScaleState(
        phospho=MatrixIntensityScaleState.log2(established_by="test.unestablished")
    )
    payload["processing_state"] = supported_log2_processing_state(
        has_total_matrix=False
    )

    with pytest.raises(TransformationValidationError, match="must be established"):
        AnalysisReadyPhosphoDataset(**payload)


def test_analysis_ready_dataset_internal_construction_does_not_bypass_validation() -> (
    None
):
    payload = _valid_payload()
    payload["site_metadata"] = payload["site_metadata"].drop(columns=["site_sequence"])

    with pytest.raises(_MODEL_BOUNDARY_ERRORS, match="site_sequence"):
        AnalysisReadyPhosphoDataset._from_owned(**payload)


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


@pytest.mark.parametrize(
    ("column_name", "position_value"),
    [
        pytest.param("position", 308, id="position-python-int"),
        pytest.param("position", np.int64(308), id="position-numpy-int"),
        pytest.param("site_position", 308, id="site-position-python-int"),
        pytest.param("site_position", np.int64(308), id="site-position-numpy-int"),
    ],
)
def test_direct_dataset_boundary_accepts_strict_integer_position_metadata(
    column_name: str,
    position_value: object,
) -> None:
    payload = _valid_payload()
    site_metadata = payload["site_metadata"].copy(deep=True)
    site_metadata.loc[:, column_name] = pd.Series(
        [182, position_value],
        index=site_metadata.index,
        dtype="object",
    )
    payload["site_metadata"] = site_metadata

    dataset = AnalysisReadyPhosphoDataset(**payload)

    assert dataset.site_metadata.loc[site_metadata.index[1], column_name] == 308


def test_direct_dataset_boundary_accepts_equivalent_position_and_site_position_metadata() -> (
    None
):
    payload = _valid_payload_with_explicit_positions(
        position_value=308,
        site_position_value=np.int64(308),
    )

    dataset = AnalysisReadyPhosphoDataset(**payload)

    row_key = dataset.site_metadata.index[1]
    assert dataset.site_metadata.loc[row_key, "position"] == 308
    assert dataset.site_metadata.loc[row_key, "site_position"] == 308


@pytest.mark.parametrize(
    "invalid_site_position",
    [
        pytest.param("308", id="numeric-string"),
        pytest.param(308.0, id="float"),
        pytest.param(True, id="boolean"),
        pytest.param(0, id="zero"),
        pytest.param(-1, id="negative"),
        pytest.param(None, id="missing"),
        pytest.param([308], id="list"),
    ],
)
def test_direct_dataset_boundary_rejects_invalid_site_position_when_position_valid(
    invalid_site_position: object,
) -> None:
    payload = _valid_payload_with_explicit_positions(
        position_value=308,
        site_position_value=invalid_site_position,
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    message = str(exc_info.value)
    assert "dataset.site_metadata" in message
    assert "site_position" in message
    assert "position" in message


@pytest.mark.parametrize(
    "invalid_position",
    [
        pytest.param("308", id="numeric-string"),
        pytest.param(308.0, id="float"),
        pytest.param(True, id="boolean"),
        pytest.param(0, id="zero"),
        pytest.param(-1, id="negative"),
        pytest.param(None, id="missing"),
        pytest.param([308], id="list"),
    ],
)
def test_direct_dataset_boundary_rejects_invalid_position_when_site_position_valid(
    invalid_position: object,
) -> None:
    payload = _valid_payload_with_explicit_positions(
        position_value=invalid_position,
        site_position_value=308,
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    message = str(exc_info.value)
    assert "dataset.site_metadata" in message
    assert "position" in message


@pytest.mark.parametrize(
    ("position_value", "site_position_value"),
    [
        pytest.param(308, 309, id="position-valid-site-position-disagrees"),
        pytest.param(309, 308, id="position-disagrees-site-position-valid"),
    ],
)
def test_direct_dataset_boundary_rejects_disagreeing_position_and_site_position_metadata(
    position_value: object,
    site_position_value: object,
) -> None:
    payload = _valid_payload_with_explicit_positions(
        position_value=position_value,
        site_position_value=site_position_value,
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    message = str(exc_info.value)
    assert "dataset.site_metadata" in message
    assert "position" in message


def test_direct_dataset_boundary_rejects_invalid_position_and_site_position_metadata() -> (
    None
):
    payload = _valid_payload_with_explicit_positions(
        position_value="308",
        site_position_value=308.0,
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    message = str(exc_info.value)
    assert "dataset.site_metadata" in message
    assert "position" in message


@pytest.mark.parametrize("column_name", ["position", "site_position"])
@pytest.mark.parametrize(
    "invalid_position",
    [
        pytest.param(True, id="boolean"),
        pytest.param(0, id="zero"),
        pytest.param(-1, id="negative"),
        pytest.param(308.0, id="float"),
        pytest.param("308", id="numeric-string"),
        pytest.param(None, id="missing"),
        pytest.param([308], id="list"),
    ],
)
def test_direct_dataset_boundary_rejects_loose_explicit_position_metadata(
    column_name: str,
    invalid_position: object,
) -> None:
    payload = _valid_payload()
    site_metadata = payload["site_metadata"].copy(deep=True)
    site_metadata.loc[:, column_name] = pd.Series(
        [182, invalid_position],
        index=site_metadata.index,
        dtype="object",
    )
    payload["site_metadata"] = site_metadata

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    message = str(exc_info.value)
    assert "dataset.site_metadata" in message
    assert column_name in message
    assert "position" in message


def test_protein_scoped_site_key_encodes_numpy_integer_like_python_int() -> None:
    python_int_key = build_protein_scoped_site_key(
        organism="rat",
        protein_namespace="protein_id",
        protein_identifier="AKT1",
        residue="T",
        position=308,
        field_name="test.site_key",
        error_type=ValueError,
    )
    numpy_int_key = build_protein_scoped_site_key(
        organism="rat",
        protein_namespace="protein_id",
        protein_identifier="AKT1",
        residue="T",
        position=np.int64(308),
        field_name="test.site_key",
        error_type=ValueError,
    )

    assert encode_site_key(python_int_key) == encode_site_key(numpy_int_key)


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


def test_model_boundary_validator_rejects_duplicate_sample_metadata_columns() -> None:
    payload = _valid_payload()
    payload["sample_metadata"] = pd.DataFrame(
        [["a", "batch_1"], ["b", "batch_2"]],
        columns=["condition", "condition"],
        index=pd.Index(["sample_a", "sample_b"], name="sample_id"),
    )

    with pytest.raises(
        _MODEL_BOUNDARY_ERRORS,
        match="dataset.sample_metadata.columns must be unique",
    ):
        AnalysisReadyPhosphoDataset(**payload)
    with pytest.raises(
        _MODEL_BOUNDARY_ERRORS,
        match="dataset.sample_metadata.columns must be unique",
    ):
        _BOUNDARY_VALIDATOR.run(**payload)


def test_model_boundary_validator_accepts_distinct_sample_metadata_columns() -> None:
    payload = _valid_payload()
    payload["sample_metadata"] = pd.DataFrame(
        {
            "condition": ["a", "b"],
            "batch": ["batch_1", "batch_2"],
        },
        index=pd.Index(["sample_a", "sample_b"], name="sample_id"),
    )

    dataset = AnalysisReadyPhosphoDataset(**payload)

    assert dataset.sample_metadata is not None
    assert list(dataset.sample_metadata.columns) == ["condition", "batch"]


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


def test_analysis_ready_numeric_matrix_missing_values_use_fast_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = _valid_payload()
    phospho = payload["phospho"].copy(deep=True)
    phospho.iloc[0, 0] = np.nan
    payload["phospho"] = phospho

    def fail_if_scalar_fallback_runs(value: object) -> bool:
        raise AssertionError(f"unexpected scalar missing-value scan for {value!r}")

    monkeypatch.setattr(
        dataset_models,
        "_is_missing_value",
        fail_if_scalar_fallback_runs,
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    assert str(exc_info.value) == _NO_MISSING_PROCESSING_STATE_ERROR


def test_analysis_ready_object_matrix_missing_values_use_fallback(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = _valid_payload()
    phospho = payload["phospho"].copy(deep=True).astype(object)
    phospho.iloc[0, 0] = None
    payload["phospho"] = phospho
    observed_values: list[object] = []
    original_is_missing_value = dataset_models._is_missing_value

    def spy_is_missing_value(value: object) -> bool:
        observed_values.append(value)
        return original_is_missing_value(value)

    monkeypatch.setattr(dataset_models, "_is_missing_value", spy_is_missing_value)

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    assert str(exc_info.value) == _NO_MISSING_PROCESSING_STATE_ERROR
    assert observed_values
    assert any(value is None for value in observed_values)


@pytest.mark.parametrize(
    ("dtype", "missing_value"),
    [
        pytest.param("float64", np.nan, id="float-nan"),
        pytest.param("Float64", pd.NA, id="nullable-float-pd-na"),
        pytest.param("object", None, id="object-none"),
    ],
)
def test_analysis_ready_missing_value_detection_preserves_error_behavior(
    dtype: str,
    missing_value: object,
) -> None:
    payload = _valid_payload()
    phospho = payload["phospho"].copy(deep=True).astype(dtype)
    phospho.iloc[0, 0] = missing_value
    payload["phospho"] = phospho

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    assert str(exc_info.value) == _NO_MISSING_PROCESSING_STATE_ERROR


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
