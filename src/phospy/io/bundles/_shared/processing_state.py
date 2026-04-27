"""Dataset-processing-state payload serialization for bundle manifests."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.datasets.processing_state import (
    ComparisonState,
    DatasetProcessingState,
    MissingDataState,
    NormalisationState,
    SiteMatrixState,
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
    validate_json_safe_mapping,
)
from phospy.transformations.models import QuantitativeMeaning


def processing_state_to_payload(state: DatasetProcessingState) -> dict[str, object]:
    """Serialize dataset processing state to manifest payload."""

    correction_diagnostics = _require_optional_json_safe_mapping(
        state.total_protein_correction.diagnostics,
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics"
        ),
    )
    return {
        "intensity_scale": intensity_scale_state_to_payload(state.intensity_scale),
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
            "diagnostics": correction_diagnostics,
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
    if "requires_log_scale" in correction_payload:
        requires_log_scale = _require_optional_bool(
            correction_payload.get("requires_log_scale"),
            field_name=(
                "dataset.metadata.processing_state.total_protein_correction."
                "requires_log_scale"
            ),
        )
    else:
        # Legacy minimal payloads encoded only policy/applied.
        requires_log_scale = False if not correction_applied else None
    inferred_quantity = _infer_legacy_quantitative_meaning(
        intensity_scale_payload=intensity_scale_payload,
        correction_payload=correction_payload,
        correction_applied=correction_applied,
    )
    intensity_scale_state = intensity_scale_state_from_payload(
        intensity_scale_payload,
        fallback_quantity=inferred_quantity,
    )
    correction_diagnostics = _require_optional_json_safe_mapping(
        correction_payload.get("diagnostics"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction.diagnostics"
        ),
    )
    correction_quantitative_meaning = _resolve_total_correction_quantitative_meaning(
        correction_payload=correction_payload,
        correction_diagnostics=correction_diagnostics,
        fallback=None,
    )
    return DatasetProcessingState(
        intensity_scale=intensity_scale_state,
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


def _require_optional_json_safe_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, object] | None:
    if value is None:
        return None
    return validate_json_safe_mapping(value, field_name=field_name)


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


def _resolve_total_correction_quantitative_meaning(
    *,
    correction_payload: Mapping[str, object],
    correction_diagnostics: Mapping[str, object] | None,
    fallback: str | None,
) -> str | None:
    direct = _require_optional_str(
        correction_payload.get("quantitative_meaning"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction."
            "quantitative_meaning"
        ),
    )
    if direct is not None:
        return direct
    if correction_diagnostics is None:
        return fallback
    from_diagnostics = _require_optional_str(
        correction_diagnostics.get("quantitative_meaning"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction."
            "diagnostics.quantitative_meaning"
        ),
    )
    if from_diagnostics is not None:
        return from_diagnostics
    legacy_output_quantity = _require_optional_str(
        correction_diagnostics.get("output_quantity"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction."
            "diagnostics.output_quantity"
        ),
    )
    if legacy_output_quantity is not None:
        return legacy_output_quantity
    return fallback


def _infer_legacy_quantitative_meaning(
    *,
    intensity_scale_payload: Mapping[str, object],
    correction_payload: Mapping[str, object],
    correction_applied: bool,
) -> QuantitativeMeaning:
    if "quantity" in intensity_scale_payload:
        return QuantitativeMeaning.UNKNOWN
    correction_quantitative_meaning = _resolve_total_correction_quantitative_meaning(
        correction_payload=correction_payload,
        correction_diagnostics=_require_optional_json_safe_mapping(
            correction_payload.get("diagnostics"),
            field_name=(
                "dataset.metadata.processing_state.total_protein_correction.diagnostics"
            ),
        ),
        fallback=None,
    )
    if correction_quantitative_meaning is not None:
        try:
            return QuantitativeMeaning(correction_quantitative_meaning)
        except ValueError:
            return QuantitativeMeaning.UNKNOWN
    correction_policy = _require_optional_str(
        correction_payload.get("policy"),
        field_name="dataset.metadata.processing_state.total_protein_correction.policy",
    )
    correction_output_scale = _require_optional_str(
        correction_payload.get("output_scale"),
        field_name=(
            "dataset.metadata.processing_state.total_protein_correction.output_scale"
        ),
    )
    if correction_applied and (
        correction_policy == "subtract_log_total"
        or correction_output_scale == "log2_ratio"
    ):
        return QuantitativeMeaning.PHOSPHO_TOTAL_LOG_RATIO
    if correction_applied:
        return QuantitativeMeaning.UNKNOWN
    phospho_payload = require_mapping(
        intensity_scale_payload.get("phospho"),
        field_name="dataset.metadata.processing_state.intensity_scale.phospho",
    )
    kind = _require_optional_str(
        phospho_payload.get("kind"),
        field_name="dataset.metadata.processing_state.intensity_scale.phospho.kind",
    )
    if kind == "linear":
        return QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE
    if kind == "log2":
        return QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE
    return QuantitativeMeaning.UNKNOWN
