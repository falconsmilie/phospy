"""Dataset-processing-state payload serialization for bundle manifests."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.datasets.processing_state import (
    ComparisonState,
    DatasetProcessingState,
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
    intensity_scale_state_to_payload,
)
from phospy.io.bundles._shared.primitives import (
    require_bool,
    require_int,
    require_mapping,
    require_str,
)


def processing_state_to_payload(state: DatasetProcessingState) -> dict[str, object]:
    """Serialize dataset processing state to manifest payload."""

    correction_diagnostics = _normalize_optional_total_correction_diagnostics(
        state.total_protein_correction.diagnostics,
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics"
        ),
    )
    return {
        "intensity_scale": intensity_scale_state_to_payload(state.intensity_scale),
        "site_sequence_resolution": {
            "configured": bool(state.site_sequence_resolution.configured),
            "mode": state.site_sequence_resolution.mode,
            "flank_size": state.site_sequence_resolution.flank_size,
            "fasta_sha256": state.site_sequence_resolution.fasta_sha256,
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
        },
        "missing_data": {
            "policy": state.missing_data.policy,
            "min_observed_values": state.missing_data.min_observed_values,
            "complete_matrix": state.missing_data.complete_matrix,
            "imputed": state.missing_data.imputed,
        },
        "normalisation": {"policy": state.normalisation.policy},
        "total_protein_correction": {
            "policy": state.total_protein_correction.policy,
            "applied": state.total_protein_correction.applied,
            "formula": state.total_protein_correction.formula,
            "requires_log_scale": state.total_protein_correction.requires_log_scale,
            "input_scale": state.total_protein_correction.input_scale,
            "output_scale": state.total_protein_correction.output_scale,
            "quantitative_meaning": (
                state.total_protein_correction.quantitative_meaning
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
    }


def processing_state_from_payload(
    payload: Mapping[str, object],
) -> DatasetProcessingState:
    """Deserialize dataset processing state from manifest payload."""

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
    minimum_observed_values = _require_optional_int(
        missing_data_payload.get("min_observed_values"),
        field_name="dataset.metadata.processing_state.missing_data.min_observed_values",
    )
    site_matrix_minimum_observed_values = _require_optional_int(
        site_matrix_payload.get("minimum_observed_values"),
        field_name="dataset.metadata.processing_state.site_matrix.minimum_observed_values",
    )
    correction_applied = require_bool(
        correction_payload.get("applied"),
        field_name="dataset.metadata.processing_state.total_protein_correction.applied",
    )
    _require_payload_key(
        correction_payload,
        key="requires_log_scale",
        field_name="dataset.metadata.processing_state.total_protein_correction",
    )
    _require_payload_key(
        correction_payload,
        key="quantitative_meaning",
        field_name="dataset.metadata.processing_state.total_protein_correction",
    )
    _require_payload_key(
        correction_payload,
        key="diagnostics",
        field_name="dataset.metadata.processing_state.total_protein_correction",
    )
    requires_log_scale = _require_optional_bool(
        correction_payload.get("requires_log_scale"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction."
            "requires_log_scale"
        ),
    )
    intensity_scale_state = intensity_scale_state_from_payload(intensity_scale_payload)
    correction_diagnostics = _parse_total_correction_diagnostics(
        correction_payload.get("diagnostics"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics"
        ),
    )
    correction_quantitative_meaning = _require_total_correction_quantitative_meaning(
        correction_payload=correction_payload,
        correction_diagnostics=correction_diagnostics,
    )
    return DatasetProcessingState(
        intensity_scale=intensity_scale_state,
        site_sequence_resolution=SiteSequenceResolutionState(
            configured=require_bool(
                site_sequence_resolution_payload.get("configured", False),
                field_name=(
                    "dataset.metadata.processing_state.site_sequence_resolution."
                    "configured"
                ),
            ),
            mode=_require_optional_str(
                site_sequence_resolution_payload.get("mode"),
                field_name=(
                    "dataset.metadata.processing_state.site_sequence_resolution.mode"
                ),
            ),
            flank_size=_require_optional_int(
                site_sequence_resolution_payload.get("flank_size"),
                field_name=(
                    "dataset.metadata.processing_state.site_sequence_resolution."
                    "flank_size"
                ),
            ),
            fasta_sha256=_require_optional_str(
                site_sequence_resolution_payload.get("fasta_sha256"),
                field_name=(
                    "dataset.metadata.processing_state.site_sequence_resolution."
                    "fasta_sha256"
                ),
            ),
            resolved_site_count=require_int(
                site_sequence_resolution_payload.get("resolved_site_count", 0),
                field_name=(
                    "dataset.metadata.processing_state.site_sequence_resolution."
                    "resolved_site_count"
                ),
            ),
            unresolved_site_count=require_int(
                site_sequence_resolution_payload.get("unresolved_site_count", 0),
                field_name=(
                    "dataset.metadata.processing_state.site_sequence_resolution."
                    "unresolved_site_count"
                ),
            ),
            unresolved_counts_by_reason=_parse_optional_string_int_mapping(
                site_sequence_resolution_payload.get("unresolved_counts_by_reason"),
                field_name=(
                    "dataset.metadata.processing_state.site_sequence_resolution."
                    "unresolved_counts_by_reason"
                ),
            ),
        ),
        missing_data=MissingDataState(
            policy=require_str(
                missing_data_payload.get("policy"),
                field_name="dataset.metadata.processing_state.missing_data.policy",
            ),
            min_observed_values=minimum_observed_values,
            complete_matrix=require_bool(
                missing_data_payload.get("complete_matrix"),
                field_name=(
                    "dataset.metadata.processing_state.missing_data.complete_matrix"
                ),
            ),
            imputed=require_bool(
                missing_data_payload.get("imputed"),
                field_name="dataset.metadata.processing_state.missing_data.imputed",
            ),
        ),
        normalisation=NormalisationState(
            policy=require_str(
                normalisation_payload.get("policy"),
                field_name="dataset.metadata.processing_state.normalisation.policy",
            )
        ),
        total_protein_correction=TotalProteinCorrectionState(
            policy=require_str(
                correction_payload.get("policy"),
                field_name=(
                    "dataset.metadata.processing_state.total_protein_correction.policy"
                ),
            ),
            applied=correction_applied,
            formula=_require_optional_str(
                correction_payload.get("formula"),
                field_name=(
                    "dataset.metadata.processing_state.total_protein_correction.formula"
                ),
            ),
            requires_log_scale=requires_log_scale,
            input_scale=_require_optional_str(
                correction_payload.get("input_scale"),
                field_name=(
                    "dataset.metadata.processing_state.total_protein_correction."
                    "input_scale"
                ),
            ),
            output_scale=_require_optional_str(
                correction_payload.get("output_scale"),
                field_name=(
                    "dataset.metadata.processing_state.total_protein_correction."
                    "output_scale"
                ),
            ),
            quantitative_meaning=correction_quantitative_meaning,
            diagnostics=correction_diagnostics,
        ),
        site_matrix=SiteMatrixState(
            policy=require_str(
                site_matrix_payload.get("policy"),
                field_name="dataset.metadata.processing_state.site_matrix.policy",
            ),
            constructed=require_bool(
                site_matrix_payload.get("constructed"),
                field_name=(
                    "dataset.metadata.processing_state.site_matrix.constructed"
                ),
            ),
            missing_data_policy=require_str(
                site_matrix_payload.get("missing_data_policy"),
                field_name=(
                    "dataset.metadata.processing_state.site_matrix.missing_data_policy"
                ),
            ),
            minimum_observed_values=site_matrix_minimum_observed_values,
            duplicate_site_policy=require_str(
                site_matrix_payload.get("duplicate_site_policy"),
                field_name=(
                    "dataset.metadata.processing_state."
                    "site_matrix.duplicate_site_policy"
                ),
            ),
        ),
        comparisons=ComparisonState(
            policy=require_str(
                comparisons_payload.get("policy"),
                field_name="dataset.metadata.processing_state.comparisons.policy",
            ),
            sample_group_column=require_str(
                comparisons_payload.get("sample_group_column"),
                field_name=(
                    "dataset.metadata.processing_state.comparisons.sample_group_column"
                ),
            ),
            pairs=_parse_optional_pairs(
                comparisons_payload.get("pairs"),
                field_name="dataset.metadata.processing_state.comparisons.pairs",
            ),
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


def _require_total_correction_quantitative_meaning(
    *,
    correction_payload: Mapping[str, object],
    correction_diagnostics: TotalProteinCorrectionDiagnostics,
) -> str:
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
    return direct


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
