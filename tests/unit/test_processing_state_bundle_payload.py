from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.datasets.processing_state import (
    ComparisonState,
    DatasetProcessingState,
    MissingDataState,
    NormalisationState,
    SiteMatrixState,
    TotalProteinCorrectionDiagnostics,
    TotalProteinCorrectionState,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
)
from phospy.io.bundles._shared.processing_state import (
    processing_state_from_payload,
    processing_state_to_payload,
)


class _CustomDiagnostic:
    pass


def _intensity_scale_state():
    return intensity_scale_state_from_payload(
        {
            "phospho": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "total": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "quantity": "phospho_total_log_ratio",
        }
    )


def _processing_state_with_diagnostics(diagnostics):
    return DatasetProcessingState(
        intensity_scale=_intensity_scale_state(),
        missing_data=MissingDataState(
            policy="forbid",
            min_observed_values=None,
            complete_matrix=True,
            imputed=False,
        ),
        normalisation=NormalisationState(policy="none"),
        total_protein_correction=TotalProteinCorrectionState(
            policy="subtract_log_total",
            applied=True,
            formula="log2_phospho - log2_total",
            requires_log_scale=True,
            input_scale="log2",
            output_scale="log2_ratio",
            quantitative_meaning="phospho_total_log_ratio",
            diagnostics=diagnostics,
        ),
        site_matrix=SiteMatrixState(
            policy="as_input",
            constructed=False,
            missing_data_policy="drop_any_missing",
            minimum_observed_values=None,
            duplicate_site_policy="max_mean_signal",
        ),
        comparisons=ComparisonState(
            policy="none",
            sample_group_column="comparison_group",
            pairs=None,
        ),
    )


def _processing_payload_with_diagnostics(diagnostics):
    return {
        "intensity_scale": {
            "phospho": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "total": {
                "kind": "log2",
                "transformed": True,
                "established_by": "bundle.fixture",
            },
            "quantity": "phospho_total_log_ratio",
        },
        "missing_data": {
            "policy": "forbid",
            "min_observed_values": None,
            "complete_matrix": True,
            "imputed": False,
        },
        "normalisation": {"policy": "none"},
        "total_protein_correction": {
            "policy": "subtract_log_total",
            "applied": True,
            "formula": "log2_phospho - log2_total",
            "requires_log_scale": True,
            "input_scale": "log2",
            "output_scale": "log2_ratio",
            "quantitative_meaning": "phospho_total_log_ratio",
            "diagnostics": diagnostics,
        },
        "site_matrix": {
            "policy": "as_input",
            "constructed": False,
            "missing_data_policy": "drop_any_missing",
            "minimum_observed_values": None,
            "duplicate_site_policy": "max_mean_signal",
        },
        "comparisons": {
            "policy": "none",
            "sample_group_column": "comparison_group",
            "pairs": None,
        },
    }


def test_processing_state_payload_round_trip_preserves_total_correction_fields() -> (
    None
):
    diagnostics = {
        "policy": "subtract_log_total",
        "requested_policy": "subtract_log_total",
        "resolved_policy": "subtract_log_total",
        "formula": "log2_phospho - log2_total",
        "requires_log_scale": True,
        "input_scale": "log2",
        "output_scale": "log2_ratio",
        "quantitative_meaning": "phospho_total_log_ratio",
        "output_quantity": "phospho_total_log_ratio",
        "matched_rows": 3,
        "hashes": {
            "input_phospho_hash": "abc123",
            "output_phospho_hash": "def456",
        },
        "notes": ["stable", "json"],
    }
    state = _processing_state_with_diagnostics(diagnostics)

    payload = processing_state_to_payload(state)
    diagnostics_payload = payload["total_protein_correction"]["diagnostics"]
    assert diagnostics_payload["diagnostics_schema_version"] == 1
    assert (
        payload["total_protein_correction"]["quantitative_meaning"]
        == "phospho_total_log_ratio"
    )
    restored = processing_state_from_payload(payload)
    correction = restored.total_protein_correction

    assert correction.policy == "subtract_log_total"
    assert correction.applied is True
    assert correction.formula == "log2_phospho - log2_total"
    assert correction.requires_log_scale is True
    assert correction.input_scale == "log2"
    assert correction.output_scale == "log2_ratio"
    assert correction.quantitative_meaning == "phospho_total_log_ratio"
    assert restored.intensity_scale.quantity.value == "phospho_total_log_ratio"
    assert isinstance(correction.diagnostics, TotalProteinCorrectionDiagnostics)
    assert correction.diagnostics is not None
    assert correction.diagnostics.to_payload() == diagnostics_payload


def test_processing_state_payload_loads_new_versioned_diagnostics() -> None:
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "requires_log_scale": True,
            "matched_rows": 2,
            "legacy_debug_note": "preserve-me",
        }
    )

    restored = processing_state_from_payload(payload)
    correction = restored.total_protein_correction

    assert correction.diagnostics is not None
    diagnostics_payload = correction.diagnostics.to_payload()
    assert diagnostics_payload["diagnostics_schema_version"] == 1
    assert diagnostics_payload["matched_rows"] == 2
    assert diagnostics_payload["legacy_debug_note"] == "preserve-me"


def test_processing_state_payload_loads_legacy_output_quantity_as_quantitative_meaning() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {"output_quantity": "phospho_total_log_ratio"}
    )
    payload["total_protein_correction"].pop("quantitative_meaning", None)

    restored = processing_state_from_payload(payload)

    assert (
        restored.total_protein_correction.quantitative_meaning
        == "phospho_total_log_ratio"
    )


def test_processing_state_payload_resaves_legacy_diagnostics_as_versioned_schema() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "output_quantity": "phospho_total_log_ratio",
            "legacy_numeric_hint": 7,
        }
    )
    loaded = processing_state_from_payload(payload)

    rewritten_payload = processing_state_to_payload(loaded)
    diagnostics = rewritten_payload["total_protein_correction"]["diagnostics"]

    assert diagnostics["diagnostics_schema_version"] == 1
    assert diagnostics["quantitative_meaning"] == "phospho_total_log_ratio"
    assert diagnostics["output_quantity"] == "phospho_total_log_ratio"
    assert diagnostics["legacy_numeric_hint"] == 7


def test_processing_state_payload_converts_tuple_diagnostics_to_json_arrays() -> None:
    state = _processing_state_with_diagnostics(
        {
            "row_ids": ("row_a", "row_b"),
            "nested": {"tokens": ("a", 1, True)},
        }
    )

    payload = processing_state_to_payload(state)
    diagnostics = payload["total_protein_correction"]["diagnostics"]

    assert diagnostics == {
        "diagnostics_schema_version": 1,
        "row_ids": ["row_a", "row_b"],
        "nested": {"tokens": ["a", 1, True]},
    }


@pytest.mark.parametrize(
    ("diagnostics", "message"),
    (
        pytest.param(
            {"bad": _CustomDiagnostic()},
            "unsupported value type",
            id="custom_object",
        ),
        pytest.param(
            {"bad": np.array([1, 2, 3])},
            "numpy.ndarray",
            id="numpy_array",
        ),
        pytest.param(
            {"bad": pd.Series([1, 2, 3])},
            "pandas.Series",
            id="pandas_series",
        ),
        pytest.param(
            {"bad": {1, 2, 3}},
            "builtins.set",
            id="set",
        ),
        pytest.param(
            {1: "value"},
            "must contain only string keys",
            id="non_string_key",
        ),
    ),
)
def test_processing_state_payload_rejects_unsupported_diagnostics(
    diagnostics,
    message: str,
) -> None:
    with pytest.raises(PhosPyInputError, match=message):
        state = _processing_state_with_diagnostics(diagnostics)
        processing_state_to_payload(state)


def test_processing_state_from_payload_rejects_non_finite_diagnostic_float() -> None:
    payload = _processing_payload_with_diagnostics({"nan_value": float("nan")})

    with pytest.raises(PhosPyInputError, match="finite float"):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_non_string_diagnostic_keys() -> None:
    payload = _processing_payload_with_diagnostics({1: "value"})

    with pytest.raises(PhosPyInputError, match="must contain only string keys"):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_unsupported_diagnostic_schema_version() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 99,
            "matched_rows": 2,
        }
    )

    with pytest.raises(
        PhosPyInputError,
        match="diagnostics_schema_version=99.*unsupported",
    ):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_malformed_versioned_diagnostics() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "matched_rows": "three",
        }
    )

    with pytest.raises(PhosPyInputError, match="matched_rows must be an int"):
        processing_state_from_payload(payload)
