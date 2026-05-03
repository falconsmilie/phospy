from __future__ import annotations

import pytest

from phospy.datasets.processing_state import (
    ComparisonState,
    DatasetProcessingState,
    MissingDataDiagnostics,
    MissingDataState,
    NormalisationState,
    SiteMatrixState,
    SiteSequenceResolutionState,
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


def _intensity_scale_state(*, quantity: str = "phospho_total_log_ratio"):
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
            "quantity": quantity,
        }
    )


def _processing_state_with_diagnostics(
    diagnostics,
    *,
    quantitative_meaning: str = "phospho_total_log_ratio",
    missing_data_diagnostics=None,
):
    return DatasetProcessingState(
        intensity_scale=_intensity_scale_state(quantity=quantitative_meaning),
        site_sequence_resolution=SiteSequenceResolutionState(
            configured=False,
            mode=None,
            flank_size=None,
            fasta_sha256=None,
            resolved_site_count=0,
            unresolved_site_count=0,
            unresolved_counts_by_reason={},
        ),
        missing_data=MissingDataState(
            policy="forbid",
            min_observed_values=None,
            complete_matrix=True,
            imputed=False,
            diagnostics=missing_data_diagnostics,
        ),
        normalisation=NormalisationState(policy="none"),
        total_protein_correction=TotalProteinCorrectionState(
            policy="subtract_log_total",
            applied=True,
            formula="log2_phospho - log2_total",
            requires_log_scale=True,
            input_scale="log2",
            output_scale="log2_ratio",
            quantitative_meaning=quantitative_meaning,
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


def _processing_payload_with_diagnostics(
    diagnostics,
    *,
    quantitative_meaning: str = "phospho_total_log_ratio",
    missing_data_diagnostics=None,
):
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
            "quantity": quantitative_meaning,
        },
        "missing_data": {
            "policy": "forbid",
            "min_observed_values": None,
            "complete_matrix": True,
            "imputed": False,
            "diagnostics": missing_data_diagnostics,
        },
        "normalisation": {"policy": "none"},
        "total_protein_correction": {
            "policy": "subtract_log_total",
            "applied": True,
            "formula": "log2_phospho - log2_total",
            "requires_log_scale": True,
            "input_scale": "log2",
            "output_scale": "log2_ratio",
            "quantitative_meaning": quantitative_meaning,
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
        "diagnostics_schema_version": 1,
        "policy": "subtract_log_total",
        "requested_policy": "subtract_log_total",
        "resolved_policy": "subtract_log_total",
        "formula": "log2_phospho - log2_total",
        "requires_log_scale": True,
        "input_scale": "log2",
        "output_scale": "log2_ratio",
        "quantitative_meaning": "phospho_total_log_ratio",
        "matched_rows": 3,
        "total_table_hash": "abc123",
        "input_phospho_hash": "def456",
        "output_phospho_hash": "ghi789",
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
            "quantitative_meaning": "phospho_total_log_ratio",
            "requires_log_scale": True,
            "matched_rows": 2,
        }
    )

    restored = processing_state_from_payload(payload)
    correction = restored.total_protein_correction

    assert correction.diagnostics is not None
    diagnostics_payload = correction.diagnostics.to_payload()
    assert diagnostics_payload["diagnostics_schema_version"] == 1
    assert diagnostics_payload["matched_rows"] == 2


def test_processing_state_payload_round_trip_preserves_mixed_total_correction_state() -> (
    None
):
    mixed_meaning = "mixed_phospho_total_log_ratio_and_phosphosite_log_abundance"
    diagnostics = {
        "diagnostics_schema_version": 1,
        "policy": "subtract_log_total",
        "requested_policy": "subtract_log_total",
        "resolved_policy": "subtract_log_total",
        "quantitative_meaning": mixed_meaning,
        "corrected_row_count": 2,
        "uncorrected_row_count": 1,
        "unmatched_policy": "allow_uncorrected",
        "corrected_phosphosite_row_ids": ["SITE_A", "SITE_B"],
        "corrected_phosphosite_to_total_protein_row_id": {
            "SITE_A": "TP_A",
            "SITE_B": "TP_B",
        },
        "unmatched_phosphosite_row_ids": ["SITE_C"],
        "uncorrected_phosphosite_row_reasons": {
            "SITE_C": "no_matching_total_protein_row_retained_by_unmatched_policy_allow_uncorrected"
        },
    }
    state = _processing_state_with_diagnostics(
        diagnostics,
        quantitative_meaning=mixed_meaning,
    )
    payload = processing_state_to_payload(state)
    restored = processing_state_from_payload(payload)
    assert restored.intensity_scale.quantity is not None
    assert restored.intensity_scale.quantity.value == mixed_meaning
    correction = restored.total_protein_correction
    assert correction.quantitative_meaning == mixed_meaning
    assert correction.diagnostics is not None
    assert (
        correction.diagnostics.to_payload()
        == (payload["total_protein_correction"]["diagnostics"])
    )


def test_processing_state_from_payload_rejects_unversioned_diagnostics() -> None:
    payload = _processing_payload_with_diagnostics(
        {
            "policy": "subtract_log_total",
            "matched_rows": 2,
        }
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction."
            "diagnostics.diagnostics_schema_version is required"
        ),
    ):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_unknown_versioned_diagnostics_fields() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "matched_rows": 2,
            "legacy_debug_note": "not-allowed",
        }
    )

    with pytest.raises(
        PhosPyInputError,
        match="contains unsupported field\\(s\\): legacy_debug_note",
    ):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_non_object_diagnostics() -> None:
    payload = _processing_payload_with_diagnostics("not-an-object")

    with pytest.raises(PhosPyInputError, match="must be an object"):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_non_string_diagnostic_keys() -> None:
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            1: "value",
        }
    )

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
            "quantitative_meaning": "phospho_total_log_ratio",
            "matched_rows": "three",
        }
    )

    with pytest.raises(PhosPyInputError, match="matched_rows must be an int"):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_missing_total_correction_diagnostics_quantitative_meaning_key() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "matched_rows": 2,
        }
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction."
            "diagnostics.quantitative_meaning is required"
        ),
    ):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_missing_total_correction_quantitative_meaning_key() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "matched_rows": 2,
        }
    )
    payload["total_protein_correction"].pop("quantitative_meaning", None)

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction."
            "quantitative_meaning is required"
        ),
    ):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_missing_total_correction_diagnostics_key() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "matched_rows": 2,
        }
    )
    payload["total_protein_correction"].pop("diagnostics", None)

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction."
            "diagnostics is required"
        ),
    ):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_applied_total_correction_with_null_diagnostics() -> (
    None
):
    payload = _processing_payload_with_diagnostics(None)

    with pytest.raises(
        PhosPyInputError,
        match=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics "
            "must be an object with"
        ),
    ):
        processing_state_from_payload(payload)


def test_processing_state_payload_round_trip_preserves_missing_data_diagnostics() -> (
    None
):
    missing_data_diagnostics = {
        "diagnostics_schema_version": 1,
        "missing_data_policy": "impute_row_median",
        "imputation_method_id": "row_median",
        "imputation_method_family": "deterministic_row_statistic",
        "input_missing_cell_count": 2,
        "output_missing_cell_count": 0,
        "imputed_cell_count": 2,
        "affected_row_count": 2,
        "affected_column_count": 2,
        "affected_row_ids": ["row_a", "row_b"],
        "affected_column_ids": ["sample_1", "sample_2"],
        "imputed_row_ids": ["row_a"],
        "imputed_column_ids": ["sample_2"],
        "dropped_row_ids": ["row_c"],
        "random_seed": None,
        "method_parameters": {"min_observed_values": 1},
        "matrix_scale_requirement": None,
        "stage_order": ["missing_data"],
        "missingness_mask_hash": "abc123",
        "left_censored_assumption": False,
        "rows_not_imputable": [],
        "per_column_distribution_parameters": {
            "sample_1": {
                "observed_count": 2,
                "missing_count": 1,
                "q": 0.01,
            }
        },
        "dropped_rows_above_max_missing_fraction": ["row_c"],
    }
    state = _processing_state_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        },
        missing_data_diagnostics=missing_data_diagnostics,
    )

    payload = processing_state_to_payload(state)
    restored = processing_state_from_payload(payload)

    assert isinstance(restored.missing_data.diagnostics, MissingDataDiagnostics)
    assert restored.missing_data.diagnostics is not None
    assert (
        restored.missing_data.diagnostics.to_payload()
        == payload["missing_data"]["diagnostics"]
    )


def test_processing_state_payload_without_missing_data_diagnostics_deserializes() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        },
    )
    payload["missing_data"].pop("diagnostics", None)

    restored = processing_state_from_payload(payload)

    assert restored.missing_data.diagnostics is None
