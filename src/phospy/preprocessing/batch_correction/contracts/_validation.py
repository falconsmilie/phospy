"""Scalar construction validation for batch-correction contracts."""

from __future__ import annotations

from collections.abc import Collection, Sequence
from enum import Enum

from phospy.errors.input import PhosPyInputError
from phospy.policies.config_values import (
    coerce_policy_enum,
    require_non_empty_string,
    require_optional_int_at_least,
    require_supported_literal,
)


def validate_batch_correction_config(
    *,
    method: object,
    batch_column: object,
    condition_column: object,
    preserve_condition_effects: object,
    supported_methods: Collection[str],
) -> None:
    """Validate public batch-correction intent config fields."""

    require_supported_literal(
        method,
        field_name="dataset build request preprocessing_config.batch_correction.method",
        supported_values=supported_methods,
        error_type=PhosPyInputError,
    )
    require_non_empty_string(
        batch_column,
        field_name=(
            "dataset build request preprocessing_config.batch_correction.batch_column"
        ),
        error_type=PhosPyInputError,
    )
    require_non_empty_string(
        condition_column,
        field_name=(
            "dataset build request preprocessing_config.batch_correction."
            "condition_column"
        ),
        error_type=PhosPyInputError,
    )
    if preserve_condition_effects is not True:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.batch_correction."
            "preserve_condition_effects must be True because "
            "linear_residualize_batch preserves condition effects by design"
        )


def validate_internal_batch_correction_request(
    *,
    method: object,
    batch_column: object,
    condition_columns: object,
    replicate_column: object | None,
    control_site_source: object,
    control_site_mode: object,
    missing_value_policy: object,
    imputation_policy: object,
    n_unwanted_factors: object | None,
    stage_order: object,
    diagnostics_enabled: object,
    method_type: type[Enum],
    control_site_source_type: type[Enum],
    control_site_mode_type: type[Enum],
    missing_value_policy_type: type[Enum],
    imputation_policy_type: type[Enum],
    stage_order_type: type[Enum],
    missing_value_policy_reject: str,
    missing_value_policy_allow_temporary_imputation: str,
    missing_value_policy_use_existing_imputation_provenance: str,
    imputation_policy_none: str,
) -> dict[str, object]:
    """Validate scalar construction invariants for an internal request."""

    prefix = "internal batch-correction request."
    resolved_method = _coerce_internal_enum(
        method_type, method, field_name=f"{prefix}method"
    )
    resolved_control_site_source = _coerce_internal_enum(
        control_site_source_type,
        control_site_source,
        field_name=f"{prefix}control_site_source",
    )
    resolved_control_site_mode = _coerce_internal_enum(
        control_site_mode_type,
        control_site_mode,
        field_name=f"{prefix}control_site_mode",
    )
    resolved_missing_value_policy = _coerce_internal_enum(
        missing_value_policy_type,
        missing_value_policy,
        field_name=f"{prefix}missing_value_policy",
    )
    resolved_imputation_policy = _coerce_internal_enum(
        imputation_policy_type,
        imputation_policy,
        field_name=f"{prefix}imputation_policy",
    )
    resolved_stage_order = _coerce_internal_enum(
        stage_order_type, stage_order, field_name=f"{prefix}stage_order"
    )
    resolved_batch_column = require_non_empty_string(
        batch_column,
        field_name=f"{prefix}batch_column",
        error_type=PhosPyInputError,
    )
    resolved_condition_columns = _require_non_empty_string_sequence(
        condition_columns,
        field_name=f"{prefix}condition_columns",
    )
    resolved_replicate_column = None
    if replicate_column is not None:
        resolved_replicate_column = require_non_empty_string(
            replicate_column,
            field_name=f"{prefix}replicate_column",
            error_type=PhosPyInputError,
            when_provided=True,
        )
    if (
        str(resolved_method.value) == "ruv_iii_style"
        and resolved_replicate_column is None
    ):
        raise PhosPyInputError(
            f"{prefix}replicate_column is required when method='ruv_iii_style'"
        )
    resolved_n_unwanted_factors = require_optional_int_at_least(
        n_unwanted_factors,
        field_name=f"{prefix}n_unwanted_factors",
        minimum=1,
        error_type=PhosPyInputError,
    )
    if (
        str(resolved_method.value) == "ruv_iii_style"
        and resolved_n_unwanted_factors is None
    ):
        raise PhosPyInputError(
            f"{prefix}n_unwanted_factors is required when method='ruv_iii_style'"
        )
    if not isinstance(diagnostics_enabled, bool):
        raise PhosPyInputError(f"{prefix}diagnostics_enabled must be a bool")

    _validate_internal_missing_imputation_pair(
        missing_value_policy=str(resolved_missing_value_policy.value),
        imputation_policy=str(resolved_imputation_policy.value),
        field_prefix=prefix,
        missing_value_policy_reject=missing_value_policy_reject,
        missing_value_policy_allow_temporary_imputation=(
            missing_value_policy_allow_temporary_imputation
        ),
        missing_value_policy_use_existing_imputation_provenance=(
            missing_value_policy_use_existing_imputation_provenance
        ),
        imputation_policy_none=imputation_policy_none,
    )
    return {
        "method": resolved_method,
        "batch_column": resolved_batch_column,
        "condition_columns": resolved_condition_columns,
        "replicate_column": resolved_replicate_column,
        "control_site_source": resolved_control_site_source,
        "control_site_mode": resolved_control_site_mode,
        "missing_value_policy": resolved_missing_value_policy,
        "imputation_policy": resolved_imputation_policy,
        "n_unwanted_factors": resolved_n_unwanted_factors,
        "stage_order": resolved_stage_order,
        "diagnostics_enabled": diagnostics_enabled,
    }


def _coerce_internal_enum(
    enum_type: type[Enum], value: object, *, field_name: str
) -> Enum:
    return coerce_policy_enum(
        enum_type,
        value,
        field_name=field_name,
        error_type=PhosPyInputError,
    )


def _require_non_empty_string_sequence(
    value: object, *, field_name: str
) -> tuple[str, ...]:
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise PhosPyInputError(
            f"{field_name} must be a non-empty sequence of non-empty strings"
        )
    resolved = tuple(value)
    if not resolved:
        raise PhosPyInputError(
            f"{field_name} must be a non-empty sequence of non-empty strings"
        )
    if any(not isinstance(item, str) or not item.strip() for item in resolved):
        raise PhosPyInputError(f"{field_name} must contain only non-empty strings")
    if len(set(resolved)) != len(resolved):
        raise PhosPyInputError(f"{field_name} must not contain duplicates")
    return resolved


def _validate_internal_missing_imputation_pair(
    *,
    missing_value_policy: str,
    imputation_policy: str,
    field_prefix: str,
    missing_value_policy_reject: str,
    missing_value_policy_allow_temporary_imputation: str,
    missing_value_policy_use_existing_imputation_provenance: str,
    imputation_policy_none: str,
) -> None:
    if (
        missing_value_policy == missing_value_policy_reject
        and imputation_policy != imputation_policy_none
    ):
        raise PhosPyInputError(
            f"{field_prefix}imputation_policy must be '{imputation_policy_none}' "
            f"when missing_value_policy='{missing_value_policy_reject}'"
        )
    if (
        missing_value_policy == missing_value_policy_allow_temporary_imputation
        and imputation_policy == imputation_policy_none
    ):
        raise PhosPyInputError(
            f"{field_prefix}imputation_policy must request a temporary "
            "imputation policy when missing_value_policy="
            f"'{missing_value_policy_allow_temporary_imputation}'"
        )
    if (
        missing_value_policy == missing_value_policy_use_existing_imputation_provenance
        and imputation_policy != imputation_policy_none
    ):
        raise PhosPyInputError(
            f"{field_prefix}imputation_policy must be '{imputation_policy_none}' "
            "when relying on existing imputation provenance"
        )


__all__ = [
    "validate_batch_correction_config",
    "validate_internal_batch_correction_request",
]
