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


def _valid_payload() -> dict[str, object]:
    site_index = pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id")
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


def test_model_boundary_validator_parity_for_misaligned_site_metadata() -> None:
    payload = _valid_payload()
    payload["site_metadata"] = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                _CENTRED_Y_SEQUENCE,
                _CENTRED_T_SEQUENCE,
            ],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T309;"], name="site_id"),
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
