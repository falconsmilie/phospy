from __future__ import annotations

from dataclasses import replace

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
)
from phospy.io.bundles._shared.processing_state import (
    processing_state_from_payload,
    processing_state_to_payload,
)
from phospy.science.datasets.preprocessing.policy_models import (
    MissingDataPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.processing_state import (
    ComparisonState,
    DatasetProcessingState,
    MissingDataDiagnostics,
    MissingDataDiagnosticsV1,
    MissingDataDiagnosticsV2,
    MissingDataState,
    NormalisationState,
    RuvReadinessState,
    SiteMatrixState,
    SiteSequenceResolutionRowDiagnostic,
    SiteSequenceResolutionState,
    TotalProteinCorrectionDiagnostics,
    TotalProteinCorrectionState,
)
from phospy.science.transformations.models import QuantitativeMeaning


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
        },
        legacy_quantitative_meaning_policy="migrate_unverified",
    )


def _processing_state_with_diagnostics(
    diagnostics,
    *,
    quantitative_meaning: str = "phospho_total_log_ratio",
    missing_data_diagnostics=None,
    ruv_readiness: RuvReadinessState | None = None,
):
    missing_policy, missing_imputed = _missing_data_policy_and_imputed_flag(
        missing_data_diagnostics
    )
    return DatasetProcessingState(
        intensity_scale=_intensity_scale_state(quantity=quantitative_meaning),
        site_sequence_resolution=SiteSequenceResolutionState(
            configured=False,
            mode=None,
            flank_size=None,
            fasta_source_path=None,
            fasta_source_label=None,
            fasta_sha256=None,
            resolver_version=None,
            resolved_site_count=0,
            unresolved_site_count=0,
            unresolved_counts_by_reason={},
            filled_missing_count=0,
            replaced_existing_count=0,
            preserved_existing_count=0,
            existing_sequence_conflict_count=0,
            conflict_policy=None,
            row_diagnostics=(),
        ),
        missing_data=MissingDataState(
            policy=missing_policy,
            min_observed_values=None,
            complete_matrix=True,
            imputed=missing_imputed,
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
            duplicate_site_policy="error",
        ),
        comparisons=ComparisonState(
            policy="none",
            sample_group_column="comparison_group",
            pairs=None,
        ),
        ruv_readiness=(
            RuvReadinessState(
                enabled=False,
                ready=False,
                reasons=("not configured",),
                control_feature_column="is_control_feature",
                replicate_group_column="replicate_group",
                batch_column="batch",
                control_feature_count=0,
                replicate_group_count=0,
                batch_count=0,
                requires_complete_matrix=True,
                matrix_complete=True,
                imputation_method_id=None,
                missingness_mask_preserved=False,
            )
            if ruv_readiness is None
            else ruv_readiness
        ),
    )


def _processing_payload_with_diagnostics(
    diagnostics,
    *,
    quantitative_meaning: str = "phospho_total_log_ratio",
    missing_data_diagnostics=None,
):
    missing_policy, missing_imputed = _missing_data_policy_and_imputed_flag(
        missing_data_diagnostics
    )
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
            "policy": missing_policy,
            "min_observed_values": None,
            "complete_matrix": True,
            "imputed": missing_imputed,
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
            "duplicate_site_policy": "error",
        },
        "comparisons": {
            "policy": "none",
            "sample_group_column": "comparison_group",
            "pairs": None,
        },
    }


def _missing_data_policy_and_imputed_flag(missing_data_diagnostics) -> tuple[str, bool]:
    if missing_data_diagnostics is None:
        return "forbid", False
    payload = (
        missing_data_diagnostics.to_payload()
        if isinstance(missing_data_diagnostics, MissingDataDiagnostics)
        else missing_data_diagnostics
    )
    if not isinstance(payload, dict) or payload.get("imputed_cell_count", 0) <= 0:
        return "forbid", False
    policy = payload.get("missing_data_policy")
    return (policy if isinstance(policy, str) else "impute_row_median"), True


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
    assert correction.policy is TotalProteinCorrectionPolicy.SUBTRACT_LOG_TOTAL
    assert (
        correction.quantitative_meaning is QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO
    )
    assert restored.missing_data.policy is MissingDataPolicy.FORBID


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


def test_processing_state_from_payload_rejects_unknown_total_correction_quantitative_meaning() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "not_a_supported_meaning",
        },
    )
    payload["total_protein_correction"]["quantitative_meaning"] = (
        "not_a_supported_meaning"
    )

    with pytest.raises(PhosPyInputError, match="must be one of:"):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_unknown_missing_data_policy() -> None:
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        },
        missing_data_diagnostics={
            "diagnostics_schema_version": 1,
            "missing_data_policy": "unknown_policy",
            "input_missing_cell_count": 0,
            "output_missing_cell_count": 0,
            "imputed_cell_count": 0,
            "affected_row_count": 0,
            "affected_column_count": 0,
            "affected_row_ids": [],
            "affected_column_ids": [],
            "imputed_row_ids": [],
            "imputed_column_ids": [],
            "dropped_row_ids": [],
            "method_parameters": {},
            "stage_order": ["missing_data"],
            "missingness_mask_hash": "hash",
            "rows_not_imputable": [],
        },
    )

    with pytest.raises(PhosPyInputError, match="must be one of:"):
        processing_state_from_payload(payload)


def test_processing_state_from_payload_rejects_unknown_missing_data_diagnostics_fields() -> (
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
        missing_data_diagnostics={
            "diagnostics_schema_version": 1,
            "missing_data_policy": "forbid",
            "input_missing_cell_count": 0,
            "output_missing_cell_count": 0,
            "imputed_cell_count": 0,
            "affected_row_count": 0,
            "affected_column_count": 0,
            "affected_row_ids": [],
            "affected_column_ids": [],
            "imputed_row_ids": [],
            "imputed_column_ids": [],
            "dropped_row_ids": [],
            "method_parameters": {},
            "stage_order": ["missing_data"],
            "missingness_mask_hash": "hash",
            "rows_not_imputable": [],
            "legacy_debug_note": "unsupported",
        },
    )

    with pytest.raises(
        PhosPyInputError,
        match="contains unsupported field\\(s\\): legacy_debug_note",
    ):
        processing_state_from_payload(payload)


def test_missing_data_diagnostics_json_round_trip_stays_stable() -> None:
    payload = {
        "diagnostics_schema_version": 1,
        "missing_data_policy": "impute_row_median",
        "imputation_method_id": "row_median",
        "imputation_method_family": "deterministic_row_statistic",
        "input_missing_cell_count": 2,
        "output_missing_cell_count": 0,
        "imputed_cell_count": 2,
        "affected_row_count": 1,
        "affected_column_count": 1,
        "affected_row_ids": ["row_a"],
        "affected_column_ids": ["sample_1"],
        "imputed_row_ids": ["row_a"],
        "imputed_column_ids": ["sample_1"],
        "dropped_row_ids": [],
        "method_parameters": {
            "min_observed_values": 1,
            "input_scale": "linear",
            "imputation_operation_order": "no_intensity_transform",
        },
        "imputation_input_scale": "linear",
        "imputation_input_scale_source": "caller_selected",
        "imputation_operation_order": "no_intensity_transform",
        "stage_order": ["missing_data"],
        "missingness_mask_hash": "mask-hash",
        "imputation_mask_hash": "imputation-mask-hash",
        "rows_not_imputable": [],
    }

    diagnostics = MissingDataDiagnostics.from_payload(
        payload,
        field_name="dataset.metadata.processing_state.missing_data.diagnostics",
    )

    assert diagnostics.to_payload() == {
        **payload,
        "row_medians_used": {},
    }


def _group_aware_v2_diagnostics_payload() -> dict[str, object]:
    return {
        "diagnostics_schema_version": 2,
        "missing_data_policy": "impute_group_aware",
        "imputation_method_id": "group_aware_knn_minprob",
        "imputation_method_family": "group_aware_mixed_mechanism",
        "input_missing_cell_count": 6,
        "output_missing_cell_count": 0,
        "imputed_cell_count": 5,
        "affected_row_count": 3,
        "affected_column_count": 4,
        "affected_row_ids": ["mixed", "absence", "unsupported"],
        "affected_column_ids": ["a1", "a2", "b1", "b2"],
        "imputed_row_ids": ["mixed", "absence"],
        "imputed_column_ids": ["a1", "a2", "b1", "b2"],
        "dropped_row_ids": ["unsupported"],
        "imputed_row_count": 2,
        "imputed_column_count": 4,
        "dropped_row_count": 1,
        "random_seed": 42,
        "method_parameters": {
            "group_column": "condition",
            "min_partial_observed_fraction": 0.75,
            "min_reference_observed_fraction": 0.75,
            "k": 2,
            "distance": "nan_euclidean",
            "no_overlap_policy": "error",
            "no_overlap_policy_version": 1,
            "q": 0.01,
            "width": 0.3,
            "seed": 42,
            "knn_target_cell_count": 1,
            "minprob_target_cell_count": 4,
            "knn_target_mask_hash": "knn-target-hash",
            "minprob_target_mask_hash": "minprob-target-hash",
            "mechanism_input": "original_retained_matrix",
            "resolved_group_samples": {
                "A": ["a1", "a2"],
                "B": ["b1", "b2"],
            },
            "input_scale": "log2",
            "imputation_operation_order": "no_intensity_transform",
        },
        "matrix_scale_requirement": "log2",
        "imputation_input_scale": "log2",
        "imputation_input_scale_source": "method_required",
        "imputation_operation_order": "no_intensity_transform",
        "stage_order": ["missing_data"],
        "missingness_mask_hash": "missingness-hash",
        "imputation_mask_hash": "overall-imputation-hash",
        "left_censored_assumption": True,
        "rows_not_imputable": ["unsupported"],
        "row_medians_used": {},
        "neighbour_count": 2,
        "distance_metric": "nan_euclidean",
        "knn_no_overlap_policy": "error",
        "knn_no_overlap_policy_version": 1,
        "dropped_rows_above_max_missing_fraction": [],
        "per_column_distribution_parameters": {
            column: {
                "observed_count": 2,
                "missing_count": 1,
                "q": 0.01,
                "width": 0.3,
                "lower_q_quantile": 1.0,
                "lower_tail_mean": 1.0,
                "observed_sd": 0.5,
                "imputation_mean": 0.73,
                "imputation_sd": 0.15,
            }
            for column in ("a1", "a2", "b1", "b2")
        },
        "group_aware": {
            "group_column": "condition",
            "observed_group_sizes": {"A": 2, "B": 2},
            "min_partial_observed_fraction": 0.75,
            "min_reference_observed_fraction": 0.75,
            "retained_row_count": 2,
            "dropped_unsupported_row_count": 1,
            "knn_routed_cell_count": 1,
            "minprob_routed_cell_count": 4,
            "knn_imputed_cell_count": 1,
            "minprob_imputed_cell_count": 4,
            "knn_target_mask_hash": "knn-target-hash",
            "minprob_target_mask_hash": "minprob-target-hash",
            "knn_imputation_mask_hash": "knn-imputation-hash",
            "minprob_imputation_mask_hash": "minprob-imputation-hash",
            "unsupported_partial_group_count": 1,
            "unsupported_absence_group_count": 0,
            "unsupported_partial_row_ids": ["unsupported"],
            "unsupported_absence_row_ids": [],
            "routed_rows": [
                {
                    "row_id": "mixed",
                    "knn_imputed_columns": ["a1"],
                    "minprob_imputed_columns": ["b1", "b2"],
                    "knn_group_labels": ["A"],
                    "minprob_group_labels": ["B"],
                },
                {
                    "row_id": "absence",
                    "knn_imputed_columns": [],
                    "minprob_imputed_columns": ["a1", "a2"],
                    "knn_group_labels": [],
                    "minprob_group_labels": ["A"],
                },
            ],
            "rejected_rows": [
                {
                    "row_id": "unsupported",
                    "unsupported_partial_groups": ["A"],
                    "unsupported_absence_groups": [],
                    "observed_finite_count_by_group": {"A": 1, "B": 2},
                    "observed_fraction_by_group": {"A": 0.5, "B": 1.0},
                }
            ],
            "minprob_left_censored_assumption": True,
            "route_categories": [
                "partial_observation_knn",
                "asymmetric_absence_minprob",
                "unsupported_partial",
                "unsupported_absence",
            ],
            "mechanism_input": "original_retained_matrix",
        },
    }


def test_group_aware_v2_diagnostics_dispatch_and_bundle_round_trip() -> None:
    diagnostics_payload = _group_aware_v2_diagnostics_payload()
    diagnostics = MissingDataDiagnostics.from_payload(
        diagnostics_payload,
        field_name="dataset.metadata.processing_state.missing_data.diagnostics",
    )
    assert isinstance(diagnostics, MissingDataDiagnosticsV2)
    assert diagnostics.group_aware.knn_imputed_cell_count == 1
    assert diagnostics.group_aware.rejected_rows[0].unsupported_partial_groups == ("A",)

    state = _processing_state_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        },
        missing_data_diagnostics=diagnostics,
    )
    restored = processing_state_from_payload(processing_state_to_payload(state))
    assert isinstance(restored.missing_data.diagnostics, MissingDataDiagnosticsV2)
    assert restored.missing_data.diagnostics.to_payload() == diagnostics.to_payload()


def test_historical_group_aware_v1_diagnostics_round_trip_without_rewriting() -> None:
    payload = _group_aware_v2_diagnostics_payload()
    payload["diagnostics_schema_version"] = 1
    payload.pop("group_aware")

    diagnostics = MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")

    assert type(diagnostics) is MissingDataDiagnosticsV1
    restored = MissingDataDiagnostics.from_payload(
        diagnostics.to_payload(), field_name="diagnostics"
    )
    assert type(restored) is MissingDataDiagnosticsV1
    assert restored.to_payload() == diagnostics.to_payload()
    assert restored.to_payload()["diagnostics_schema_version"] == 1
    assert "group_aware" not in restored.to_payload()


def test_missing_data_diagnostics_schema_versions_reject_unknown_fields() -> None:
    v1_payload = _group_aware_v2_diagnostics_payload()
    v1_payload["diagnostics_schema_version"] = 1
    with pytest.raises(PhosPyInputError, match="group_aware"):
        MissingDataDiagnosticsV1.from_mapping(v1_payload, field_name="diagnostics")

    v2_payload = _group_aware_v2_diagnostics_payload()
    assert isinstance(v2_payload["group_aware"], dict)
    v2_payload["group_aware"]["unknown_route_fact"] = True
    with pytest.raises(PhosPyInputError, match="unknown_route_fact"):
        MissingDataDiagnostics.from_payload(v2_payload, field_name="diagnostics")


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("overlapping_cells", "both mechanisms"),
        ("unknown_column", "resolved_group_samples"),
        ("wrong_group", "must match routed columns"),
    ],
)
def test_group_aware_v2_rejects_invalid_routed_cell_attribution(
    mutation: str, message: str
) -> None:
    payload = _group_aware_v2_diagnostics_payload()
    group_aware = payload["group_aware"]
    assert isinstance(group_aware, dict)
    routed_rows = group_aware["routed_rows"]
    assert isinstance(routed_rows, list)
    first = routed_rows[0]
    assert isinstance(first, dict)
    if mutation == "overlapping_cells":
        first["minprob_imputed_columns"] = ["a1", "b1", "b2"]
    elif mutation == "unknown_column":
        first["knn_imputed_columns"] = ["unknown_sample"]
        payload["imputed_column_ids"] = ["unknown_sample", "a1", "a2", "b1", "b2"]
        payload["imputed_column_count"] = 5
    else:
        second = routed_rows[1]
        assert isinstance(second, dict)
        second["minprob_group_labels"] = ["B"]

    with pytest.raises(PhosPyInputError, match=message):
        MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("missing_evidence", "observed_finite_count_by_group"),
        ("inconsistent_fraction", "observed fractions must match counts"),
        ("partial_above_threshold", "below the threshold"),
        ("absence_with_reference", "lack a sufficiently observed reference"),
    ],
)
def test_group_aware_v2_rejects_invalid_unsupported_route_evidence(
    mutation: str, message: str
) -> None:
    payload = _group_aware_v2_diagnostics_payload()
    group_aware = payload["group_aware"]
    assert isinstance(group_aware, dict)
    rejected_rows = group_aware["rejected_rows"]
    assert isinstance(rejected_rows, list)
    rejected = rejected_rows[0]
    assert isinstance(rejected, dict)
    if mutation == "missing_evidence":
        rejected.pop("observed_finite_count_by_group")
    elif mutation == "inconsistent_fraction":
        rejected["observed_fraction_by_group"]["A"] = 0.25
    elif mutation == "partial_above_threshold":
        group_aware["min_partial_observed_fraction"] = 0.5
        payload["method_parameters"]["min_partial_observed_fraction"] = 0.5
    else:
        rejected["unsupported_partial_groups"] = []
        rejected["unsupported_absence_groups"] = ["A"]
        rejected["observed_finite_count_by_group"] = {"A": 0, "B": 2}
        rejected["observed_fraction_by_group"] = {"A": 0.0, "B": 1.0}
        group_aware["unsupported_partial_group_count"] = 0
        group_aware["unsupported_absence_group_count"] = 1
        group_aware["unsupported_partial_row_ids"] = []
        group_aware["unsupported_absence_row_ids"] = ["unsupported"]

    with pytest.raises(PhosPyInputError, match=message):
        MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")


def test_group_aware_v2_routing_records_are_recursively_immutable() -> None:
    payload = _group_aware_v2_diagnostics_payload()
    diagnostics = MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")
    assert isinstance(diagnostics, MissingDataDiagnosticsV2)
    expected = diagnostics.to_payload()

    payload["group_aware"]["routed_rows"][0]["knn_imputed_columns"].append("a2")
    payload["group_aware"]["rejected_rows"][0]["observed_fraction_by_group"]["A"] = 1.0
    assert diagnostics.to_payload() == expected

    serialized = diagnostics.to_payload()
    serialized["group_aware"]["routed_rows"][0]["knn_group_labels"].append("B")
    serialized["group_aware"]["rejected_rows"][0]["observed_finite_count_by_group"][
        "A"
    ] = 2
    assert diagnostics.to_payload() == expected


@pytest.mark.parametrize(
    ("section", "field_name", "invalid_value"),
    [
        ("method_parameters", "q", None),
        ("method_parameters", "width", 0.0),
        ("method_parameters", "k", 3),
        ("method_parameters", "seed", 43),
        ("method_parameters", "group_column", "batch"),
        ("method_parameters", "min_partial_observed_fraction", 0.6),
        ("method_parameters", "no_overlap_policy", "column_mean"),
        ("top", "neighbour_count", None),
        ("top", "distance_metric", "euclidean"),
        ("top", "knn_no_overlap_policy", None),
        ("top", "matrix_scale_requirement", None),
        ("top", "imputation_input_scale", None),
        ("top", "imputation_operation_order", None),
        ("top", "per_column_distribution_parameters", None),
        ("top", "dropped_row_ids", []),
        ("top", "rows_not_imputable", []),
        ("group_aware", "group_column", "batch"),
    ],
)
def test_group_aware_v2_rejects_incomplete_or_contradictory_provenance(
    section: str,
    field_name: str,
    invalid_value: object,
) -> None:
    payload = _group_aware_v2_diagnostics_payload()
    target = payload if section == "top" else payload[section]
    assert isinstance(target, dict)
    if invalid_value is None:
        target.pop(field_name)
    else:
        target[field_name] = invalid_value

    with pytest.raises(PhosPyInputError):
        MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")


@pytest.mark.parametrize(
    "rejected_group_field",
    ["unsupported_partial_groups", "unsupported_absence_groups"],
)
def test_group_aware_v2_rejects_unobserved_rejected_groups(
    rejected_group_field: str,
) -> None:
    payload = _group_aware_v2_diagnostics_payload()
    group_aware = payload["group_aware"]
    assert isinstance(group_aware, dict)
    rejected_rows = group_aware["rejected_rows"]
    assert isinstance(rejected_rows, list)
    rejected_row = rejected_rows[0]
    assert isinstance(rejected_row, dict)
    rejected_row[rejected_group_field] = ["NOT_AN_OBSERVED_GROUP"]

    with pytest.raises(PhosPyInputError, match="observed_group_sizes"):
        MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")


@pytest.mark.parametrize(
    ("partial_groups", "absence_groups", "message"),
    [
        (["A", "A"], [], "must contain unique groups"),
        ([], ["A", "A"], "must contain unique groups"),
        (["A"], ["A"], "both unsupported categories"),
    ],
)
def test_group_aware_v2_rejects_duplicate_or_overlapping_rejected_groups(
    partial_groups: list[str],
    absence_groups: list[str],
    message: str,
) -> None:
    payload = _group_aware_v2_diagnostics_payload()
    group_aware = payload["group_aware"]
    assert isinstance(group_aware, dict)
    rejected_rows = group_aware["rejected_rows"]
    assert isinstance(rejected_rows, list)
    rejected_row = rejected_rows[0]
    assert isinstance(rejected_row, dict)
    rejected_row["unsupported_partial_groups"] = partial_groups
    rejected_row["unsupported_absence_groups"] = absence_groups

    with pytest.raises(PhosPyInputError, match=message):
        MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")


def test_group_aware_v2_rejects_imputed_row_that_is_also_rejected() -> None:
    payload = _group_aware_v2_diagnostics_payload()
    payload["imputed_row_ids"] = ["mixed", "unsupported"]

    with pytest.raises(PhosPyInputError, match="disjoint from rejected/dropped rows"):
        MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")


def test_group_aware_v2_direct_construction_enforces_routing_identities() -> None:
    diagnostics = MissingDataDiagnostics.from_payload(
        _group_aware_v2_diagnostics_payload(), field_name="diagnostics"
    )
    assert isinstance(diagnostics, MissingDataDiagnosticsV2)
    group_aware = diagnostics.group_aware
    rejected_row = group_aware.rejected_rows[0]

    with pytest.raises(PhosPyInputError, match="must contain unique groups"):
        replace(rejected_row, unsupported_partial_groups=("A", "A"))

    with pytest.raises(PhosPyInputError, match="both unsupported categories"):
        replace(rejected_row, unsupported_absence_groups=("A",))

    unobserved_record = replace(
        rejected_row, unsupported_partial_groups=("NOT_AN_OBSERVED_GROUP",)
    )
    with pytest.raises(PhosPyInputError, match="observed_group_sizes"):
        replace(group_aware, rejected_rows=(unobserved_record,))

    with pytest.raises(PhosPyInputError, match="disjoint from rejected/dropped rows"):
        replace(diagnostics, imputed_row_ids=("mixed", "unsupported"))


def test_group_aware_v2_direct_construction_enforces_required_minprob_facts() -> None:
    diagnostics = MissingDataDiagnostics.from_payload(
        _group_aware_v2_diagnostics_payload(), field_name="diagnostics"
    )
    assert isinstance(diagnostics, MissingDataDiagnosticsV2)
    method_parameters = diagnostics.to_payload()["method_parameters"]
    assert isinstance(method_parameters, dict)
    method_parameters.pop("q")

    with pytest.raises(PhosPyInputError, match=r"method_parameters\.q"):
        replace(diagnostics, method_parameters=method_parameters)


def test_group_aware_v2_rejects_incomplete_per_column_distribution_provenance() -> None:
    payload = _group_aware_v2_diagnostics_payload()
    distributions = payload["per_column_distribution_parameters"]
    assert isinstance(distributions, dict)
    distributions["a1"] = {"q": 0.01, "width": 0.3}

    with pytest.raises(PhosPyInputError, match="observed_count"):
        MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")

    payload = _group_aware_v2_diagnostics_payload()
    distributions = payload["per_column_distribution_parameters"]
    assert isinstance(distributions, dict)
    a1_parameters = distributions["a1"]
    assert isinstance(a1_parameters, dict)
    a1_parameters["q"] = 0.02

    with pytest.raises(PhosPyInputError, match="must agree with method q"):
        MissingDataDiagnostics.from_payload(payload, field_name="diagnostics")


def test_processing_state_rejects_policy_that_contradicts_diagnostics() -> None:
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
            "requires_log_scale": True,
            "matched_rows": 2,
        },
        missing_data_diagnostics=_group_aware_v2_diagnostics_payload(),
    )
    missing_data = payload["missing_data"]
    assert isinstance(missing_data, dict)
    missing_data["policy"] = "impute_knn"

    with pytest.raises(
        PhosPyInputError,
        match=r"missing_data\.policy must match.*diagnostics\.missing_data_policy",
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
        "method_parameters": {
            "min_observed_values": 1,
            "input_scale": "linear",
            "imputation_operation_order": "no_intensity_transform",
        },
        "matrix_scale_requirement": None,
        "imputation_input_scale": "linear",
        "imputation_input_scale_source": "caller_selected",
        "imputation_operation_order": "no_intensity_transform",
        "stage_order": ["missing_data"],
        "missingness_mask_hash": "abc123",
        "imputation_mask_hash": "imputation-mask-hash",
        "left_censored_assumption": False,
        "rows_not_imputable": [],
        "row_medians_used": {"row_a": 1.25},
        "neighbour_count": 3,
        "distance_metric": "nan_euclidean",
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
    assert restored.missing_data.imputation_input_scale == "linear"
    assert restored.missing_data.imputation_operation_order == "no_intensity_transform"
    assert (
        restored.missing_data.diagnostics.to_payload()
        == payload["missing_data"]["diagnostics"]
    )


def test_processing_state_payload_round_trip_preserves_knn_fallback_diagnostics() -> (
    None
):
    missing_data_diagnostics = {
        "diagnostics_schema_version": 1,
        "missing_data_policy": "impute_knn",
        "imputation_method_id": "knn",
        "imputation_method_family": "nearest_neighbour",
        "input_missing_cell_count": 2,
        "output_missing_cell_count": 0,
        "imputed_cell_count": 2,
        "affected_row_count": 1,
        "affected_column_count": 2,
        "affected_row_ids": ["target"],
        "affected_column_ids": ["sample_b", "sample_c"],
        "imputed_row_ids": ["target"],
        "imputed_column_ids": ["sample_b", "sample_c"],
        "dropped_row_ids": [],
        "random_seed": None,
        "method_parameters": {
            "k": 1,
            "distance": "nan_euclidean",
            "max_missing_fraction_per_row": 1.0,
            "no_overlap_policy": "column_mean_with_caveat",
            "no_overlap_policy_version": 1,
            "input_scale": "linear",
            "imputation_operation_order": "no_intensity_transform",
        },
        "matrix_scale_requirement": None,
        "imputation_input_scale": "linear",
        "imputation_input_scale_source": "caller_selected",
        "imputation_operation_order": "no_intensity_transform",
        "stage_order": ["missing_data"],
        "missingness_mask_hash": "missingness-mask-hash",
        "imputation_mask_hash": "aggregate-mask-hash",
        "left_censored_assumption": False,
        "rows_not_imputable": [],
        "row_medians_used": {},
        "dropped_rows_above_max_missing_fraction": [],
        "neighbour_count": 1,
        "distance_metric": "nan_euclidean",
        "knn_no_overlap_policy": "column_mean_with_caveat",
        "knn_no_overlap_policy_version": 1,
        "knn_nearest_neighbour_imputed_cell_count": 1,
        "knn_nearest_neighbour_imputed_row_ids": ["target"],
        "knn_nearest_neighbour_imputed_column_ids": ["sample_b"],
        "knn_column_mean_fallback_imputed_cell_count": 1,
        "knn_column_mean_fallback_row_ids": ["target"],
        "knn_column_mean_fallback_column_ids": ["sample_c"],
        "knn_nearest_neighbour_imputation_mask_hash": "nearest-mask-hash",
        "knn_column_mean_fallback_imputation_mask_hash": "fallback-mask-hash",
        "knn_fully_column_mean_fallback_row_ids": [],
        "diagnostic_caveat_codes": ["knn_column_mean_fallback_used"],
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
    restored_payload = restored.missing_data.diagnostics.to_payload()
    assert restored_payload["knn_column_mean_fallback_imputed_cell_count"] == 1
    assert restored_payload["diagnostic_caveat_codes"] == [
        "knn_column_mean_fallback_used"
    ]


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


def test_processing_state_payload_missing_data_diagnostics_defaults_row_medians_used_for_legacy_payload() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        }
    )
    payload["missing_data"]["diagnostics"] = {
        "diagnostics_schema_version": 1,
        "missing_data_policy": "impute_row_median",
        "imputation_method_id": "row_median",
        "imputation_method_family": "deterministic_row_statistic",
        "input_missing_cell_count": 1,
        "output_missing_cell_count": 0,
        "imputed_cell_count": 1,
        "affected_row_count": 1,
        "affected_column_count": 1,
        "affected_row_ids": ["row_a"],
        "affected_column_ids": ["sample_1"],
        "imputed_row_ids": ["row_a"],
        "imputed_column_ids": ["sample_1"],
        "dropped_row_ids": [],
        "method_parameters": {"min_observed_values": 1},
        "matrix_scale_requirement": None,
        "stage_order": ["missing_data"],
        "missingness_mask_hash": "legacy-mask-hash",
        "imputation_mask_hash": "legacy-imputation-mask-hash",
        "left_censored_assumption": False,
        "rows_not_imputable": [],
        "dropped_rows_above_max_missing_fraction": [],
        "neighbour_count": None,
        "distance_metric": None,
    }
    payload["missing_data"]["policy"] = "impute_row_median"
    payload["missing_data"]["imputed"] = True

    restored = processing_state_from_payload(payload)

    assert restored.missing_data.diagnostics is not None
    assert restored.missing_data.diagnostics.to_payload()["row_medians_used"] == {}


def test_processing_state_payload_round_trip_preserves_ruv_readiness() -> None:
    state = _processing_state_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        },
        ruv_readiness=RuvReadinessState(
            enabled=True,
            ready=True,
            reasons=(),
            control_feature_column="is_control_feature",
            replicate_group_column="replicate_group",
            batch_column="batch",
            control_feature_count=3,
            replicate_group_count=2,
            batch_count=2,
            requires_complete_matrix=True,
            matrix_complete=True,
            imputation_method_id="row_median",
            missingness_mask_preserved=True,
        ),
    )

    payload = processing_state_to_payload(state)
    restored = processing_state_from_payload(payload)

    assert restored.ruv_readiness.enabled is True
    assert restored.ruv_readiness.ready is True
    assert restored.ruv_readiness.reasons == ()
    assert restored.ruv_readiness.imputation_method_id == "row_median"
    assert restored.ruv_readiness.missingness_mask_preserved is True


def test_processing_state_payload_without_ruv_readiness_uses_backward_compatible_default() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        }
    )

    restored = processing_state_from_payload(payload)

    assert restored.ruv_readiness.enabled is False
    assert restored.ruv_readiness.ready is False
    assert "not configured" in set(restored.ruv_readiness.reasons)


def test_processing_state_payload_round_trip_preserves_site_sequence_resolution_fields() -> (
    None
):
    state = _processing_state_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        }
    )
    state = DatasetProcessingState(
        intensity_scale=state.intensity_scale,
        site_sequence_resolution=SiteSequenceResolutionState(
            configured=True,
            mode="replace_existing",
            flank_size=7,
            fasta_source_path="C:/data/proteome.fasta",
            fasta_source_label="dataset.site_sequence_resolution",
            fasta_sha256="abcdef123456",
            resolver_version="phospy.science.sequences.resolver.v1",
            resolved_site_count=11,
            unresolved_site_count=3,
            unresolved_counts_by_reason={"missing_accession": 2, "site_not_found": 1},
            filled_missing_count=5,
            replaced_existing_count=2,
            preserved_existing_count=4,
            existing_sequence_conflict_count=2,
            conflict_policy="replace_existing",
            row_diagnostics=(
                SiteSequenceResolutionRowDiagnostic(
                    row_index=0,
                    row_id="MAPK14;S5;",
                    site_id="MAPK14;S5;",
                    status="resolved",
                    existing_site_sequence=None,
                    fasta_site_sequence="AASAA",
                    resolved_site_sequence="AASAA",
                    action="fill_missing",
                    reason="missing site_sequence resolved from FASTA",
                    conflict_policy="replace_existing",
                    resolver_version="phospy.science.sequences.resolver.v1",
                    fasta_source_path="C:/data/proteome.fasta",
                    fasta_sha256="abcdef123456",
                ),
                SiteSequenceResolutionRowDiagnostic(
                    row_index=1,
                    row_id="GSK3B;T6;",
                    site_id="GSK3B;T6;",
                    status="existing_sequence_conflict",
                    existing_site_sequence="XXXXX",
                    fasta_site_sequence="CCTCC",
                    resolved_site_sequence="CCTCC",
                    action="replace_existing",
                    reason="existing site_sequence conflicts with FASTA-derived sequence",
                    conflict_policy="replace_existing",
                    resolver_version="phospy.science.sequences.resolver.v1",
                    fasta_source_path="C:/data/proteome.fasta",
                    fasta_sha256="abcdef123456",
                ),
            ),
        ),
        missing_data=state.missing_data,
        normalisation=state.normalisation,
        total_protein_correction=state.total_protein_correction,
        site_matrix=state.site_matrix,
        comparisons=state.comparisons,
        ruv_readiness=state.ruv_readiness,
    )

    payload = processing_state_to_payload(state)
    restored = processing_state_from_payload(payload)
    resolution = restored.site_sequence_resolution
    assert resolution.configured is True
    assert resolution.mode == "replace_existing"
    assert resolution.flank_size == 7
    assert resolution.fasta_source_path == "C:/data/proteome.fasta"
    assert resolution.fasta_source_label == "dataset.site_sequence_resolution"
    assert resolution.fasta_sha256 == "abcdef123456"
    assert resolution.resolver_version == "phospy.science.sequences.resolver.v1"
    assert resolution.resolved_site_count == 11
    assert resolution.unresolved_site_count == 3
    assert resolution.unresolved_counts_by_reason == {
        "missing_accession": 2,
        "site_not_found": 1,
    }
    assert resolution.filled_missing_count == 5
    assert resolution.replaced_existing_count == 2
    assert resolution.preserved_existing_count == 4
    assert resolution.existing_sequence_conflict_count == 2
    assert resolution.conflict_policy == "replace_existing"
    assert len(resolution.row_diagnostics) == 2
    assert resolution.row_diagnostics[1].action == "replace_existing"
    assert resolution.row_diagnostics[1].fasta_site_sequence == "CCTCC"


def test_processing_state_payload_without_site_sequence_resolution_uses_backward_compatible_default() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        }
    )

    restored = processing_state_from_payload(payload)
    resolution = restored.site_sequence_resolution
    assert resolution.configured is False
    assert resolution.mode is None
    assert resolution.flank_size is None
    assert resolution.fasta_source_path is None
    assert resolution.fasta_source_label is None
    assert resolution.fasta_sha256 is None
    assert resolution.resolver_version is None
    assert resolution.resolved_site_count == 0
    assert resolution.unresolved_site_count == 0
    assert resolution.unresolved_counts_by_reason == {}
    assert resolution.filled_missing_count == 0
    assert resolution.replaced_existing_count == 0
    assert resolution.preserved_existing_count == 0
    assert resolution.existing_sequence_conflict_count == 0
    assert resolution.conflict_policy is None
    assert resolution.row_diagnostics == ()


def test_processing_state_payload_with_legacy_site_sequence_resolution_fields_deserializes() -> (
    None
):
    payload = _processing_payload_with_diagnostics(
        {
            "diagnostics_schema_version": 1,
            "policy": "subtract_log_total",
            "requested_policy": "subtract_log_total",
            "resolved_policy": "subtract_log_total",
            "quantitative_meaning": "phospho_total_log_ratio",
        }
    )
    payload["site_sequence_resolution"] = {
        "configured": True,
        "mode": "fill_missing_only",
        "flank_size": 5,
        "fasta_sha256": "legacy_digest",
        "resolved_site_count": 10,
        "unresolved_site_count": 1,
        "unresolved_counts_by_reason": {"missing_accession": 1},
        "conflict_policy": "preserve_existing",
        "row_diagnostics": [
            {
                "row_index": 0,
                "row_id": "MAPK14;S5;",
                "site_id": "MAPK14;S5;",
                "status": "existing_sequence_conflict",
                "existing_site_sequence": "XXXXX",
                "fasta_site_sequence": "AASAA",
                "resolved_site_sequence": "XXXXX",
                "action": "preserve_existing",
                "reason": "existing site_sequence conflicts with FASTA-derived sequence",
                "conflict_policy": "preserve_existing",
                "resolver_version": "phospy.science.sequences.resolver.v1",
                "fasta_source_path": "C:/data/proteome.fasta",
                "fasta_sha256": "legacy_digest",
            }
        ],
    }

    restored = processing_state_from_payload(payload)
    resolution = restored.site_sequence_resolution
    assert resolution.configured is True
    assert resolution.mode == "fill_missing_only"
    assert resolution.flank_size == 5
    assert resolution.fasta_sha256 == "legacy_digest"
    assert resolution.resolved_site_count == 10
    assert resolution.unresolved_site_count == 1
    assert resolution.unresolved_counts_by_reason == {"missing_accession": 1}
    assert resolution.fasta_source_path is None
    assert resolution.fasta_source_label is None
    assert resolution.resolver_version is None
    assert resolution.filled_missing_count == 0
    assert resolution.replaced_existing_count == 0
    assert resolution.preserved_existing_count == 0
    assert resolution.existing_sequence_conflict_count == 1
    assert resolution.conflict_policy == "preserve_existing"
    assert len(resolution.row_diagnostics) == 1
    assert resolution.row_diagnostics[0].action == "preserve_existing"
