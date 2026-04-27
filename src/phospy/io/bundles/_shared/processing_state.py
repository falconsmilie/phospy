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
)


def processing_state_to_payload(state: DatasetProcessingState) -> dict[str, object]:
    """Serialize dataset processing state to manifest payload."""

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
            "diagnostics": (
                None
                if state.total_protein_correction.diagnostics is None
                else dict(state.total_protein_correction.diagnostics)
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
    return DatasetProcessingState(
        intensity_scale=intensity_scale_state_from_payload(
            require_mapping(
                payload.get("intensity_scale"),
                field_name="dataset.metadata.processing_state.intensity_scale",
            )
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
            diagnostics=_require_optional_mapping(
                correction_payload.get("diagnostics"),
                field_name=(
                    "dataset.metadata.processing_state.total_protein_correction."
                    "diagnostics"
                ),
            ),
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


def _require_optional_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, object] | None:
    if value is None:
        return None
    payload = require_mapping(value, field_name=field_name)
    return {str(key): item for key, item in payload.items()}


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
