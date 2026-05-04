from __future__ import annotations

import pytest

from phospy.datasets.processing_state import (
    ComparisonState,
    DatasetProcessingState,
    MissingDataDiagnostics,
    MissingDataState,
    NormalisationState,
    RuvReadinessState,
    SiteMatrixState,
    SiteSequenceResolutionRowDiagnostic,
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
    ruv_readiness: RuvReadinessState | None = None,
):
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
        "left_censored_assumption": False,
        "rows_not_imputable": [],
        "dropped_rows_above_max_missing_fraction": [],
        "neighbour_count": None,
        "distance_metric": None,
    }

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
            resolver_version="phospy.sequences.resolver.v1",
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
                    resolver_version="phospy.sequences.resolver.v1",
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
                    resolver_version="phospy.sequences.resolver.v1",
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
    assert resolution.resolver_version == "phospy.sequences.resolver.v1"
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
                "resolver_version": "phospy.sequences.resolver.v1",
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
    assert resolution.existing_sequence_conflict_count == 0
    assert resolution.conflict_policy == "preserve_existing"
    assert len(resolution.row_diagnostics) == 1
    assert resolution.row_diagnostics[0].action == "preserve_existing"
