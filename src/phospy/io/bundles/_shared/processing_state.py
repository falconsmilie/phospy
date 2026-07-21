"""Dataset-processing-state payload serialization for bundle manifests."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_from_payload,
    intensity_scale_state_to_payload,
)
from phospy.io.bundles._shared.primitives import (
    require_bool,
    require_int,
    require_mapping,
    require_str,
)
from phospy.science.datasets.preprocessing.policy_models import (
    MissingDataPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.datasets.processing_state import (
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
    default_ruv_readiness_state,
)
from phospy.science.transformations.models import QuantitativeMeaning


@dataclass(frozen=True, slots=True)
class _ProcessingStatePayloads:
    missing_data: Mapping[str, object]
    normalisation: Mapping[str, object]
    total_protein_correction: Mapping[str, object]
    site_matrix: Mapping[str, object]
    comparisons: Mapping[str, object]
    site_sequence_resolution: Mapping[str, object]
    intensity_scale: Mapping[str, object]
    ruv_readiness_raw: object


def processing_state_to_payload(state: DatasetProcessingState) -> dict[str, object]:
    """Serialize dataset processing state to manifest payload."""

    correction_diagnostics = _normalize_optional_total_correction_diagnostics(
        state.total_protein_correction.diagnostics,
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics"
        ),
    )
    missing_data_diagnostics = _normalize_optional_missing_data_diagnostics(
        state.missing_data.diagnostics,
        field_name="dataset.metadata.processing_state.missing_data.diagnostics",
    )
    return {
        "intensity_scale": intensity_scale_state_to_payload(state.intensity_scale),
        "site_sequence_resolution": {
            "configured": bool(state.site_sequence_resolution.configured),
            "mode": state.site_sequence_resolution.mode,
            "flank_size": state.site_sequence_resolution.flank_size,
            "fasta_source_path": state.site_sequence_resolution.fasta_source_path,
            "fasta_source_label": state.site_sequence_resolution.fasta_source_label,
            "fasta_sha256": state.site_sequence_resolution.fasta_sha256,
            "resolver_version": state.site_sequence_resolution.resolver_version,
            "resolved_site_count": int(
                state.site_sequence_resolution.resolved_site_count
            ),
            "unresolved_site_count": int(
                state.site_sequence_resolution.unresolved_site_count
            ),
            "unresolved_counts_by_reason": {
                str(key): int(value)
                for key, value in state.site_sequence_resolution.unresolved_counts_by_reason.items()
            },
            "filled_missing_count": int(
                state.site_sequence_resolution.filled_missing_count
            ),
            "replaced_existing_count": int(
                state.site_sequence_resolution.replaced_existing_count
            ),
            "preserved_existing_count": int(
                state.site_sequence_resolution.preserved_existing_count
            ),
            "existing_sequence_conflict_count": int(
                state.site_sequence_resolution.existing_sequence_conflict_count
            ),
            "conflict_policy": state.site_sequence_resolution.conflict_policy,
            "row_diagnostics": [
                _site_sequence_row_diagnostic_to_payload(item)
                for item in state.site_sequence_resolution.row_diagnostics
            ],
        },
        "missing_data": {
            "policy": state.missing_data.policy.value,
            "min_observed_values": state.missing_data.min_observed_values,
            "complete_matrix": state.missing_data.complete_matrix,
            "imputed": state.missing_data.imputed,
            "has_missing_values": state.missing_data.has_missing_values,
            "missing_value_count": state.missing_data.missing_value_count,
            "diagnostics": (
                None
                if missing_data_diagnostics is None
                else missing_data_diagnostics.to_payload()
            ),
        },
        "normalisation": {"policy": str(state.normalisation.policy)},
        "total_protein_correction": {
            "policy": state.total_protein_correction.policy.value,
            "applied": state.total_protein_correction.applied,
            "formula": state.total_protein_correction.formula,
            "requires_log_scale": state.total_protein_correction.requires_log_scale,
            "input_scale": state.total_protein_correction.input_scale,
            "output_scale": state.total_protein_correction.output_scale,
            "quantitative_meaning": (
                None
                if state.total_protein_correction.quantitative_meaning is None
                else state.total_protein_correction.quantitative_meaning.value
            ),
            "diagnostics": (
                None
                if correction_diagnostics is None
                else correction_diagnostics.to_payload()
            ),
        },
        "site_matrix": {
            "policy": state.site_matrix.policy,
            "constructed": state.site_matrix.constructed,
            "missing_data_policy": state.site_matrix.missing_data_policy,
            "minimum_observed_values": state.site_matrix.minimum_observed_values,
            "duplicate_site_policy": state.site_matrix.duplicate_site_policy,
        },
        "comparisons": {
            "policy": state.comparisons.policy,
            "sample_group_column": state.comparisons.sample_group_column,
            "pairs": (
                None
                if state.comparisons.pairs is None
                else [list(pair) for pair in state.comparisons.pairs]
            ),
        },
        "ruv_readiness": {
            "enabled": state.ruv_readiness.enabled,
            "ready": state.ruv_readiness.ready,
            "reasons": list(state.ruv_readiness.reasons),
            "control_feature_column": state.ruv_readiness.control_feature_column,
            "replicate_group_column": state.ruv_readiness.replicate_group_column,
            "batch_column": state.ruv_readiness.batch_column,
            "control_feature_count": state.ruv_readiness.control_feature_count,
            "replicate_group_count": state.ruv_readiness.replicate_group_count,
            "batch_count": state.ruv_readiness.batch_count,
            "requires_complete_matrix": state.ruv_readiness.requires_complete_matrix,
            "matrix_complete": state.ruv_readiness.matrix_complete,
            "imputation_method_id": state.ruv_readiness.imputation_method_id,
            "missingness_mask_preserved": state.ruv_readiness.missingness_mask_preserved,
        },
    }


def processing_state_from_payload(
    payload: Mapping[str, object],
) -> DatasetProcessingState:
    """Deserialize dataset processing state from manifest payload."""

    payloads = _require_processing_state_payloads(payload)
    minimum_observed_values = _require_optional_int(
        payloads.missing_data.get("min_observed_values"),
        field_name="dataset.metadata.processing_state.missing_data.min_observed_values",
    )
    missing_data_diagnostics = _parse_optional_missing_data_diagnostics(
        payloads.missing_data.get("diagnostics"),
        field_name="dataset.metadata.processing_state.missing_data.diagnostics",
    )
    site_matrix_minimum_observed_values = _require_optional_int(
        payloads.site_matrix.get("minimum_observed_values"),
        field_name="dataset.metadata.processing_state.site_matrix.minimum_observed_values",
    )
    correction_applied = require_bool(
        payloads.total_protein_correction.get("applied"),
        field_name="dataset.metadata.processing_state.total_protein_correction.applied",
    )
    _require_payload_key(
        payloads.total_protein_correction,
        key="requires_log_scale",
        field_name="dataset.metadata.processing_state.total_protein_correction",
    )
    _require_payload_key(
        payloads.total_protein_correction,
        key="quantitative_meaning",
        field_name="dataset.metadata.processing_state.total_protein_correction",
    )
    _require_payload_key(
        payloads.total_protein_correction,
        key="diagnostics",
        field_name="dataset.metadata.processing_state.total_protein_correction",
    )
    requires_log_scale = _require_optional_bool(
        payloads.total_protein_correction.get("requires_log_scale"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction."
            "requires_log_scale"
        ),
    )
    intensity_scale_state = intensity_scale_state_from_payload(
        payloads.intensity_scale,
        legacy_quantitative_meaning_policy="migrate_unverified",
    )
    correction_diagnostics = _parse_total_correction_diagnostics(
        payloads.total_protein_correction.get("diagnostics"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics"
        ),
    )
    correction_quantitative_meaning = _require_total_correction_quantitative_meaning(
        correction_payload=payloads.total_protein_correction,
        correction_diagnostics=correction_diagnostics,
    )
    return DatasetProcessingState(
        intensity_scale=intensity_scale_state,
        site_sequence_resolution=_parse_site_sequence_resolution_state(
            payloads.site_sequence_resolution
        ),
        missing_data=_parse_missing_data_state(
            payloads.missing_data,
            minimum_observed_values=minimum_observed_values,
            diagnostics=missing_data_diagnostics,
        ),
        normalisation=_parse_normalisation_state(payloads.normalisation),
        total_protein_correction=_parse_total_protein_correction_state(
            payloads.total_protein_correction,
            applied=correction_applied,
            requires_log_scale=requires_log_scale,
            quantitative_meaning=correction_quantitative_meaning,
            diagnostics=correction_diagnostics,
        ),
        site_matrix=_parse_site_matrix_state(
            payloads.site_matrix,
            minimum_observed_values=site_matrix_minimum_observed_values,
        ),
        comparisons=_parse_comparison_state(payloads.comparisons),
        ruv_readiness=_parse_ruv_readiness_state(payloads.ruv_readiness_raw),
    )


def _require_processing_state_payloads(
    payload: Mapping[str, object],
) -> _ProcessingStatePayloads:
    missing_data_payload = require_mapping(
        payload.get("missing_data"),
        field_name="dataset.metadata.processing_state.missing_data",
    )
    normalisation_payload = require_mapping(
        payload.get("normalisation"),
        field_name="dataset.metadata.processing_state.normalisation",
    )
    correction_payload = require_mapping(
        payload.get("total_protein_correction"),
        field_name="dataset.metadata.processing_state.total_protein_correction",
    )
    site_matrix_payload = require_mapping(
        payload.get("site_matrix"),
        field_name="dataset.metadata.processing_state.site_matrix",
    )
    comparisons_payload = require_mapping(
        payload.get("comparisons"),
        field_name="dataset.metadata.processing_state.comparisons",
    )
    site_sequence_resolution_raw = payload.get("site_sequence_resolution")
    site_sequence_resolution_payload = (
        {}
        if site_sequence_resolution_raw is None
        else require_mapping(
            site_sequence_resolution_raw,
            field_name="dataset.metadata.processing_state.site_sequence_resolution",
        )
    )
    intensity_scale_payload = require_mapping(
        payload.get("intensity_scale"),
        field_name="dataset.metadata.processing_state.intensity_scale",
    )
    return _ProcessingStatePayloads(
        missing_data=missing_data_payload,
        normalisation=normalisation_payload,
        total_protein_correction=correction_payload,
        site_matrix=site_matrix_payload,
        comparisons=comparisons_payload,
        site_sequence_resolution=site_sequence_resolution_payload,
        intensity_scale=intensity_scale_payload,
        ruv_readiness_raw=payload.get("ruv_readiness"),
    )


def _parse_site_sequence_resolution_state(
    payload: Mapping[str, object],
) -> SiteSequenceResolutionState:
    row_diagnostics = _parse_optional_site_sequence_row_diagnostics(
        payload.get("row_diagnostics"),
        field_name=(
            "dataset.metadata.processing_state.site_sequence_resolution.row_diagnostics"
        ),
    )
    raw_existing_conflict_count = payload.get("existing_sequence_conflict_count")
    existing_sequence_conflict_count = (
        _count_site_sequence_conflict_rows(row_diagnostics)
        if raw_existing_conflict_count is None
        else require_int(
            raw_existing_conflict_count,
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "existing_sequence_conflict_count"
            ),
        )
    )
    return SiteSequenceResolutionState(
        configured=require_bool(
            payload.get("configured", False),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution.configured"
            ),
        ),
        mode=_require_optional_str(
            payload.get("mode"),
            field_name="dataset.metadata.processing_state.site_sequence_resolution.mode",
        ),
        flank_size=_require_optional_int(
            payload.get("flank_size"),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution.flank_size"
            ),
        ),
        fasta_source_path=_require_optional_str(
            payload.get("fasta_source_path"),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "fasta_source_path"
            ),
        ),
        fasta_source_label=_require_optional_str(
            payload.get("fasta_source_label"),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "fasta_source_label"
            ),
        ),
        fasta_sha256=_require_optional_str(
            payload.get("fasta_sha256"),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "fasta_sha256"
            ),
        ),
        resolver_version=_require_optional_str(
            payload.get("resolver_version"),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "resolver_version"
            ),
        ),
        resolved_site_count=require_int(
            payload.get("resolved_site_count", 0),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "resolved_site_count"
            ),
        ),
        unresolved_site_count=require_int(
            payload.get("unresolved_site_count", 0),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "unresolved_site_count"
            ),
        ),
        unresolved_counts_by_reason=_parse_optional_string_int_mapping(
            payload.get("unresolved_counts_by_reason"),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "unresolved_counts_by_reason"
            ),
        ),
        filled_missing_count=require_int(
            payload.get("filled_missing_count", 0),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "filled_missing_count"
            ),
        ),
        replaced_existing_count=require_int(
            payload.get("replaced_existing_count", 0),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "replaced_existing_count"
            ),
        ),
        preserved_existing_count=require_int(
            payload.get("preserved_existing_count", 0),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "preserved_existing_count"
            ),
        ),
        existing_sequence_conflict_count=existing_sequence_conflict_count,
        conflict_policy=_require_optional_str(
            payload.get("conflict_policy"),
            field_name=(
                "dataset.metadata.processing_state.site_sequence_resolution."
                "conflict_policy"
            ),
        ),
        row_diagnostics=row_diagnostics,
    )


def _parse_missing_data_state(
    payload: Mapping[str, object],
    *,
    minimum_observed_values: int | None,
    diagnostics: MissingDataDiagnostics | None,
) -> MissingDataState:
    has_missing_values = _require_optional_bool(
        payload.get("has_missing_values"),
        field_name="dataset.metadata.processing_state.missing_data.has_missing_values",
    )
    missing_value_count = _require_optional_non_negative_int(
        payload.get("missing_value_count"),
        field_name="dataset.metadata.processing_state.missing_data.missing_value_count",
    )
    return MissingDataState(
        policy=MissingDataPolicy.parse(
            require_str(
                payload.get("policy"),
                field_name="dataset.metadata.processing_state.missing_data.policy",
            ),
            field_name="dataset.metadata.processing_state.missing_data.policy",
        ),
        min_observed_values=minimum_observed_values,
        complete_matrix=require_bool(
            payload.get("complete_matrix"),
            field_name="dataset.metadata.processing_state.missing_data.complete_matrix",
        ),
        imputed=require_bool(
            payload.get("imputed"),
            field_name="dataset.metadata.processing_state.missing_data.imputed",
        ),
        diagnostics=diagnostics,
        has_missing_values=has_missing_values,
        missing_value_count=missing_value_count,
    )


def _parse_normalisation_state(
    payload: Mapping[str, object],
) -> NormalisationState:
    return NormalisationState(
        policy=require_str(
            payload.get("policy"),
            field_name="dataset.metadata.processing_state.normalisation.policy",
        )
    )


def _parse_total_protein_correction_state(
    payload: Mapping[str, object],
    *,
    applied: bool,
    requires_log_scale: bool | None,
    quantitative_meaning: QuantitativeMeaning,
    diagnostics: TotalProteinCorrectionDiagnostics,
) -> TotalProteinCorrectionState:
    return TotalProteinCorrectionState(
        policy=TotalProteinCorrectionPolicy.parse(
            require_str(
                payload.get("policy"),
                field_name=(
                    "dataset.metadata.processing_state.total_protein_correction.policy"
                ),
            ),
            field_name=(
                "dataset.metadata.processing_state.total_protein_correction.policy"
            ),
        ),
        applied=applied,
        formula=_require_optional_str(
            payload.get("formula"),
            field_name="dataset.metadata.processing_state.total_protein_correction.formula",
        ),
        requires_log_scale=requires_log_scale,
        input_scale=_require_optional_str(
            payload.get("input_scale"),
            field_name=(
                "dataset.metadata.processing_state.total_protein_correction.input_scale"
            ),
        ),
        output_scale=_require_optional_str(
            payload.get("output_scale"),
            field_name=(
                "dataset.metadata.processing_state.total_protein_correction."
                "output_scale"
            ),
        ),
        quantitative_meaning=quantitative_meaning,
        diagnostics=diagnostics,
    )


def _parse_site_matrix_state(
    payload: Mapping[str, object],
    *,
    minimum_observed_values: int | None,
) -> SiteMatrixState:
    return SiteMatrixState(
        policy=require_str(
            payload.get("policy"),
            field_name="dataset.metadata.processing_state.site_matrix.policy",
        ),
        constructed=require_bool(
            payload.get("constructed"),
            field_name="dataset.metadata.processing_state.site_matrix.constructed",
        ),
        missing_data_policy=require_str(
            payload.get("missing_data_policy"),
            field_name=(
                "dataset.metadata.processing_state.site_matrix.missing_data_policy"
            ),
        ),
        minimum_observed_values=minimum_observed_values,
        duplicate_site_policy=require_str(
            payload.get("duplicate_site_policy"),
            field_name=(
                "dataset.metadata.processing_state.site_matrix.duplicate_site_policy"
            ),
        ),
    )


def _parse_comparison_state(
    payload: Mapping[str, object],
) -> ComparisonState:
    return ComparisonState(
        policy=require_str(
            payload.get("policy"),
            field_name="dataset.metadata.processing_state.comparisons.policy",
        ),
        sample_group_column=require_str(
            payload.get("sample_group_column"),
            field_name=(
                "dataset.metadata.processing_state.comparisons.sample_group_column"
            ),
        ),
        pairs=_parse_optional_pairs(
            payload.get("pairs"),
            field_name="dataset.metadata.processing_state.comparisons.pairs",
        ),
    )


def _require_optional_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return require_int(value, field_name=field_name)


def _require_optional_bool(value: object, *, field_name: str) -> bool | None:
    if value is None:
        return None
    return require_bool(value, field_name=field_name)


def _require_optional_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return require_str(value, field_name=field_name)


def _require_payload_key(
    payload: Mapping[str, object],
    *,
    key: str,
    field_name: str,
) -> None:
    if key not in payload:
        raise PhosPyInputError(f"{field_name}.{key} is required")


def _parse_total_correction_diagnostics(
    value: object,
    *,
    field_name: str,
) -> TotalProteinCorrectionDiagnostics:
    if value is None:
        raise PhosPyInputError(
            f"{field_name} must be an object with "
            f"{field_name}.diagnostics_schema_version"
        )
    return TotalProteinCorrectionDiagnostics.from_payload(value, field_name=field_name)


def _normalize_optional_total_correction_diagnostics(
    value: object,
    *,
    field_name: str,
) -> TotalProteinCorrectionDiagnostics | None:
    if value is None:
        return None
    if isinstance(value, TotalProteinCorrectionDiagnostics):
        return value
    return TotalProteinCorrectionDiagnostics.from_payload(value, field_name=field_name)


def _parse_optional_missing_data_diagnostics(
    value: object,
    *,
    field_name: str,
) -> MissingDataDiagnostics | None:
    if value is None:
        return None
    return _coerce_missing_data_diagnostics(value, field_name=field_name)


def _normalize_optional_missing_data_diagnostics(
    value: object,
    *,
    field_name: str,
) -> MissingDataDiagnostics | None:
    """Normalise versioned missing-data diagnostics payloads for bundle I/O."""

    if value is None:
        return None
    return _coerce_missing_data_diagnostics(value, field_name=field_name)


def _coerce_missing_data_diagnostics(
    value: object,
    *,
    field_name: str,
) -> MissingDataDiagnostics:
    payload = value.to_payload() if isinstance(value, MissingDataDiagnostics) else value
    return MissingDataDiagnostics.from_payload(payload, field_name=field_name)


def _parse_optional_pairs(
    value: object,
    *,
    field_name: str,
) -> tuple[tuple[str, str], ...] | None:
    if value is None:
        return None
    if not isinstance(value, list):
        raise PhosPyInputError(f"{field_name} must be an array of [left, right] pairs")
    parsed: list[tuple[str, str]] = []
    for item in value:
        if (
            not isinstance(item, list)
            or len(item) != 2
            or not isinstance(item[0], str)
            or not isinstance(item[1], str)
        ):
            raise PhosPyInputError(
                f"{field_name} must contain only [left_group, right_group] string pairs"
            )
        left = item[0].strip()
        right = item[1].strip()
        if not left or not right:
            raise PhosPyInputError(
                f"{field_name} entries must contain non-empty group names"
            )
        parsed.append((left, right))
    return tuple(parsed)


def _parse_ruv_readiness_state(value: object) -> RuvReadinessState:
    if value is None:
        return default_ruv_readiness_state()
    payload = require_mapping(
        value,
        field_name="dataset.metadata.processing_state.ruv_readiness",
    )
    reasons_raw = payload.get("reasons")
    if reasons_raw is None:
        reasons = ()
    else:
        reasons = _parse_required_string_tuple(
            reasons_raw,
            field_name="dataset.metadata.processing_state.ruv_readiness.reasons",
        )
    control_feature_count = _require_non_negative_int(
        payload.get("control_feature_count"),
        field_name=(
            "dataset.metadata.processing_state.ruv_readiness.control_feature_count"
        ),
    )
    replicate_group_count = _require_non_negative_int(
        payload.get("replicate_group_count"),
        field_name=(
            "dataset.metadata.processing_state.ruv_readiness.replicate_group_count"
        ),
    )
    batch_count = _require_optional_non_negative_int(
        payload.get("batch_count"),
        field_name="dataset.metadata.processing_state.ruv_readiness.batch_count",
    )
    return RuvReadinessState(
        enabled=require_bool(
            payload.get("enabled"),
            field_name="dataset.metadata.processing_state.ruv_readiness.enabled",
        ),
        ready=require_bool(
            payload.get("ready"),
            field_name="dataset.metadata.processing_state.ruv_readiness.ready",
        ),
        reasons=reasons,
        control_feature_column=require_str(
            payload.get("control_feature_column"),
            field_name=(
                "dataset.metadata.processing_state.ruv_readiness.control_feature_column"
            ),
        ),
        replicate_group_column=require_str(
            payload.get("replicate_group_column"),
            field_name=(
                "dataset.metadata.processing_state.ruv_readiness.replicate_group_column"
            ),
        ),
        batch_column=_require_optional_str(
            payload.get("batch_column"),
            field_name="dataset.metadata.processing_state.ruv_readiness.batch_column",
        ),
        control_feature_count=control_feature_count,
        replicate_group_count=replicate_group_count,
        batch_count=batch_count,
        requires_complete_matrix=require_bool(
            payload.get("requires_complete_matrix"),
            field_name=(
                "dataset.metadata.processing_state.ruv_readiness."
                "requires_complete_matrix"
            ),
        ),
        matrix_complete=require_bool(
            payload.get("matrix_complete"),
            field_name="dataset.metadata.processing_state.ruv_readiness.matrix_complete",
        ),
        imputation_method_id=_require_optional_str(
            payload.get("imputation_method_id"),
            field_name=(
                "dataset.metadata.processing_state.ruv_readiness.imputation_method_id"
            ),
        ),
        missingness_mask_preserved=require_bool(
            payload.get("missingness_mask_preserved"),
            field_name=(
                "dataset.metadata.processing_state.ruv_readiness."
                "missingness_mask_preserved"
            ),
        ),
    )


def _require_total_correction_quantitative_meaning(
    *,
    correction_payload: Mapping[str, object],
    correction_diagnostics: TotalProteinCorrectionDiagnostics,
) -> QuantitativeMeaning:
    direct = require_str(
        correction_payload.get("quantitative_meaning"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction."
            "quantitative_meaning"
        ),
    )
    if "quantitative_meaning" not in correction_diagnostics:
        raise PhosPyInputError(
            "dataset.metadata.processing_state.total_protein_correction."
            "diagnostics.quantitative_meaning is required"
        )
    from_diagnostics = require_str(
        correction_diagnostics.get("quantitative_meaning"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction."
            "diagnostics.quantitative_meaning"
        ),
    )
    if from_diagnostics != direct:
        raise PhosPyInputError(
            "dataset.metadata.processing_state.total_protein_correction."
            "quantitative_meaning must match "
            "dataset.metadata.processing_state.total_protein_correction.diagnostics."
            "quantitative_meaning"
        )
    try:
        return QuantitativeMeaning(direct)
    except ValueError as exc:
        supported = ", ".join(member.value for member in QuantitativeMeaning)
        raise PhosPyInputError(
            "dataset.metadata.processing_state.total_protein_correction."
            "quantitative_meaning must be one of: "
            f"{supported}"
        ) from exc


def _parse_optional_string_int_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, int]:
    if value is None:
        return {}
    mapping = require_mapping(value, field_name=field_name)
    parsed: dict[str, int] = {}
    for raw_key, raw_value in mapping.items():
        key = require_str(raw_key, field_name=f"{field_name}.<key>")
        parsed[key] = require_int(raw_value, field_name=f"{field_name}.{key}")
    return parsed


def _site_sequence_row_diagnostic_to_payload(
    item: SiteSequenceResolutionRowDiagnostic,
) -> dict[str, object]:
    return {
        "row_index": int(item.row_index),
        "row_id": item.row_id,
        "site_id": item.site_id,
        "status": item.status,
        "existing_site_sequence": item.existing_site_sequence,
        "fasta_site_sequence": item.fasta_site_sequence,
        "resolved_site_sequence": item.resolved_site_sequence,
        "action": item.action,
        "reason": item.reason,
        "conflict_policy": item.conflict_policy,
        "resolver_version": item.resolver_version,
        "fasta_source_path": item.fasta_source_path,
        "fasta_sha256": item.fasta_sha256,
    }


def _parse_optional_site_sequence_row_diagnostics(
    value: object,
    *,
    field_name: str,
) -> tuple[SiteSequenceResolutionRowDiagnostic, ...]:
    if value is None:
        return ()
    if not isinstance(value, list):
        raise PhosPyInputError(f"{field_name} must be an array of diagnostic objects")
    parsed: list[SiteSequenceResolutionRowDiagnostic] = []
    for position, item in enumerate(value):
        entry_field = f"{field_name}[{position}]"
        payload = require_mapping(item, field_name=entry_field)
        parsed.append(
            SiteSequenceResolutionRowDiagnostic(
                row_index=require_int(
                    payload.get("row_index"),
                    field_name=f"{entry_field}.row_index",
                ),
                row_id=require_str(
                    payload.get("row_id"),
                    field_name=f"{entry_field}.row_id",
                ),
                site_id=_require_optional_str(
                    payload.get("site_id"),
                    field_name=f"{entry_field}.site_id",
                ),
                status=require_str(
                    payload.get("status"),
                    field_name=f"{entry_field}.status",
                ),
                existing_site_sequence=_require_optional_str(
                    payload.get("existing_site_sequence"),
                    field_name=f"{entry_field}.existing_site_sequence",
                ),
                fasta_site_sequence=_require_optional_str(
                    payload.get("fasta_site_sequence"),
                    field_name=f"{entry_field}.fasta_site_sequence",
                ),
                resolved_site_sequence=_require_optional_str(
                    payload.get("resolved_site_sequence"),
                    field_name=f"{entry_field}.resolved_site_sequence",
                ),
                action=require_str(
                    payload.get("action"),
                    field_name=f"{entry_field}.action",
                ),
                reason=_require_optional_str(
                    payload.get("reason"),
                    field_name=f"{entry_field}.reason",
                ),
                conflict_policy=_require_optional_str(
                    payload.get("conflict_policy"),
                    field_name=f"{entry_field}.conflict_policy",
                ),
                resolver_version=_require_optional_str(
                    payload.get("resolver_version"),
                    field_name=f"{entry_field}.resolver_version",
                ),
                fasta_source_path=_require_optional_str(
                    payload.get("fasta_source_path"),
                    field_name=f"{entry_field}.fasta_source_path",
                ),
                fasta_sha256=_require_optional_str(
                    payload.get("fasta_sha256"),
                    field_name=f"{entry_field}.fasta_sha256",
                ),
            )
        )
    return tuple(parsed)


def _count_site_sequence_conflict_rows(
    row_diagnostics: tuple[SiteSequenceResolutionRowDiagnostic, ...],
) -> int:
    return sum(
        1
        for item in row_diagnostics
        if "conflict" in item.status.lower()
        or (item.reason is not None and "conflict" in item.reason.lower())
    )


def _parse_required_string_tuple(
    value: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        raise PhosPyInputError(f"{field_name} must be an array of strings")
    parsed: list[str] = []
    for position, item in enumerate(value):
        parsed.append(require_str(item, field_name=f"{field_name}[{position}]"))
    return tuple(parsed)


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    parsed = require_int(value, field_name=field_name)
    if parsed < 0:
        raise PhosPyInputError(f"{field_name} must be >= 0")
    return parsed


def _require_optional_non_negative_int(value: object, *, field_name: str) -> int | None:
    if value is None:
        return None
    return _require_non_negative_int(value, field_name=field_name)
