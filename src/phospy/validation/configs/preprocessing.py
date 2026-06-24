"""Preprocessing public-config validation helpers."""

from __future__ import annotations

import math
from collections.abc import Collection, Sequence
from enum import Enum

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.validation.common.config_values import (
    coerce_policy_enum,
    require_instance,
    require_non_empty_string,
    require_supported_literal,
)
from phospy.validation.common.numbers import require_optional_int_at_least
from phospy.validation.common.paths import require_local_filesystem_path

_INCOMPATIBLE_SITE_MATRIX_MISSING_DATA_POLICIES = frozenset(
    {"retain_missing", "require_min_observed_values"}
)


def validate_preprocessing_section_type(
    value: object,
    *,
    field_name: str,
    expected_type: type[object],
) -> None:
    """Validate one nested preprocessing section object type."""

    require_instance(
        value,
        expected_type=expected_type,
        field_name=field_name,
        error_type=PhosPyInputError,
    )


def validate_intensity_transform_config(
    *,
    policy: object,
    pseudocount: object,
    supported_policies: Collection[str],
) -> None:
    """Validate public intensity-transform config fields."""

    require_supported_literal(
        policy,
        field_name=(
            "dataset build request preprocessing_config.intensity_transform.policy"
        ),
        supported_values=supported_policies,
        error_type=PhosPyInputError,
    )

    if isinstance(pseudocount, bool) or not isinstance(pseudocount, (int, float)):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.intensity_transform."
            "pseudocount must be a float or int"
        )
    if not math.isfinite(float(pseudocount)):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.intensity_transform."
            "pseudocount must be finite"
        )
    if pseudocount < 0:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.intensity_transform."
            "pseudocount must be greater than or equal to 0"
        )


def validate_normalisation_config(
    *,
    policy: object,
    supported_policies: Collection[str],
) -> None:
    """Validate public normalisation config fields."""

    require_supported_literal(
        policy,
        field_name="dataset build request preprocessing_config.normalisation.policy",
        supported_values=supported_policies,
        error_type=PhosPyInputError,
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
    """Validate internal future SPS/RUV-style correction request fields.

    This guard covers scalar construction invariants only. It does not inspect
    data matrices, select controls, estimate factors, or decide whether a
    dataset is scientifically eligible for correction.
    """

    prefix = "internal batch-correction request."
    resolved_method = _coerce_internal_enum(
        method_type,
        method,
        field_name=f"{prefix}method",
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
        stage_order_type,
        stage_order,
        field_name=f"{prefix}stage_order",
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


def validate_group_coverage_filter_config(
    *,
    enabled: object,
    group_column: object | None,
    min_finite_observations_per_group: object | None,
    min_finite_fraction_per_group: object | None,
    min_groups_passing_threshold: object,
) -> None:
    """Validate public group-aware coverage filter config fields."""

    prefix = "dataset build request preprocessing_config.group_coverage_filter."
    if not isinstance(enabled, bool):
        raise PhosPyInputError(f"{prefix}enabled must be a bool")

    if enabled:
        require_non_empty_string(
            group_column,
            field_name=f"{prefix}group_column",
            error_type=PhosPyInputError,
        )
    elif group_column is not None:
        require_non_empty_string(
            group_column,
            field_name=f"{prefix}group_column",
            error_type=PhosPyInputError,
            when_provided=True,
        )

    if isinstance(min_groups_passing_threshold, bool) or not isinstance(
        min_groups_passing_threshold, int
    ):
        raise PhosPyInputError(f"{prefix}min_groups_passing_threshold must be an int")
    if min_groups_passing_threshold < 1:
        raise PhosPyInputError(
            f"{prefix}min_groups_passing_threshold must be greater than or equal to 1"
        )

    has_count_threshold = min_finite_observations_per_group is not None
    has_fraction_threshold = min_finite_fraction_per_group is not None

    if has_count_threshold and has_fraction_threshold:
        raise PhosPyInputError(
            f"{prefix}min_finite_observations_per_group and "
            "min_finite_fraction_per_group are mutually exclusive"
        )
    if enabled and not has_count_threshold and not has_fraction_threshold:
        raise PhosPyInputError(
            f"{prefix}must set exactly one finite-observation threshold when enabled"
        )

    if has_count_threshold:
        if isinstance(min_finite_observations_per_group, bool) or not isinstance(
            min_finite_observations_per_group, int
        ):
            raise PhosPyInputError(
                f"{prefix}min_finite_observations_per_group must be an int"
            )
        if min_finite_observations_per_group < 1:
            raise PhosPyInputError(
                f"{prefix}min_finite_observations_per_group must be greater than "
                "or equal to 1"
            )

    if has_fraction_threshold:
        if isinstance(min_finite_fraction_per_group, bool) or not isinstance(
            min_finite_fraction_per_group, (int, float)
        ):
            raise PhosPyInputError(
                f"{prefix}min_finite_fraction_per_group must be a float"
            )
        fraction = float(min_finite_fraction_per_group)
        if not math.isfinite(fraction) or not (0.0 < fraction <= 1.0):
            raise PhosPyInputError(
                f"{prefix}min_finite_fraction_per_group must satisfy 0 < value <= 1"
            )


def validate_protein_aware_preparation_config(
    *,
    policy: object,
    protein_mapping_policy: object,
    supported_policies: Collection[str],
    supported_mapping_policies: Collection[str],
) -> None:
    """Validate public protein-aware preparation intent config fields."""

    require_supported_literal(
        policy,
        field_name=(
            "dataset build request preprocessing_config.protein_aware_preparation."
            "policy"
        ),
        supported_values=supported_policies,
        error_type=PhosPyInputError,
    )
    require_supported_literal(
        protein_mapping_policy,
        field_name=(
            "dataset build request preprocessing_config.protein_aware_preparation."
            "protein_mapping_policy"
        ),
        supported_values=supported_mapping_policies,
        error_type=PhosPyInputError,
    )


def validate_protein_aware_sample_alignment_config(
    *,
    protein_mapping_policy: object,
    allow_reordered_samples: object,
    supported_mapping_policies: Collection[str],
) -> None:
    """Validate protein-aware alignment diagnostics configuration fields."""

    require_supported_literal(
        protein_mapping_policy,
        field_name=(
            "protein-aware sample alignment diagnostics protein_mapping_policy"
        ),
        supported_values=supported_mapping_policies,
        error_type=PhosPyInputError,
    )
    if not isinstance(allow_reordered_samples, bool):
        raise PhosPyInputError(
            "protein-aware sample alignment diagnostics allow_reordered_samples "
            "must be a bool"
        )


def validate_localisation_config(
    *,
    mode: object,
    min_confidence: object,
    confidence_column: object,
    waiver_reason: object | None,
    supported_modes: Collection[str],
    mode_allow_missing_with_waiver: str,
) -> None:
    """Validate public localisation-confidence config fields."""

    resolved_mode = require_supported_literal(
        mode,
        field_name="dataset build request preprocessing_config.localisation.mode",
        supported_values=supported_modes,
        error_type=PhosPyInputError,
    )
    if isinstance(min_confidence, bool) or not isinstance(min_confidence, (int, float)):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.localisation.min_confidence "
            "must be a float between 0.0 and 1.0"
        )
    min_confidence_value = float(min_confidence)
    if not math.isfinite(min_confidence_value) or not (
        0.0 <= min_confidence_value <= 1.0
    ):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.localisation.min_confidence "
            "must be between 0.0 and 1.0"
        )

    require_non_empty_string(
        confidence_column,
        field_name=(
            "dataset build request preprocessing_config.localisation.confidence_column"
        ),
        error_type=PhosPyInputError,
    )

    if resolved_mode != mode_allow_missing_with_waiver:
        if waiver_reason is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.localisation."
                "waiver_reason must be None unless "
                "localisation.mode='allow_missing_with_waiver'"
            )
        return

    require_non_empty_string(
        waiver_reason,
        field_name=(
            "dataset build request preprocessing_config.localisation.waiver_reason"
        ),
        error_type=PhosPyInputError,
        when_provided=True,
    )


def validate_missing_data_config(
    *,
    policy: object,
    min_observed_values: object | None,
    q: object | None,
    width: object | None,
    seed: object | None,
    k: object | None,
    distance: object | None,
    max_missing_fraction_per_row: object | None,
    supported_policies: Collection[str],
    policy_forbid: str,
    policy_impute_row_median: str,
    policy_impute_minprob: str,
    policy_impute_knn: str,
) -> None:
    """Validate public missing-data config fields."""

    resolved_policy = require_supported_literal(
        policy,
        field_name="dataset build request preprocessing_config.missing_data.policy",
        supported_values=supported_policies,
        error_type=PhosPyInputError,
    )
    if resolved_policy == policy_forbid:
        if min_observed_values is not None:
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.min_observed_values must be "
                "None when missing_data.policy='forbid'"
            )
        if (
            q is not None
            or width is not None
            or seed is not None
            or k is not None
            or distance is not None
            or max_missing_fraction_per_row is not None
        ):
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.q, .width, .seed, .k, "
                ".distance, and "
                ".max_missing_fraction_per_row must be None when "
                "missing_data.policy='forbid'"
            )
        return

    if resolved_policy == policy_impute_row_median:
        if isinstance(min_observed_values, bool) or not isinstance(
            min_observed_values, int
        ):
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.min_observed_values must be an "
                "int when missing_data.policy='impute_row_median'"
            )
        if min_observed_values < 1:
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.min_observed_values "
                "must be greater than or equal to 1 when "
                "missing_data.policy='impute_row_median'"
            )
        if (
            q is not None
            or width is not None
            or seed is not None
            or k is not None
            or distance is not None
            or max_missing_fraction_per_row is not None
        ):
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.q, .width, .seed, .k, "
                ".distance, and "
                ".max_missing_fraction_per_row must be None when "
                "missing_data.policy='impute_row_median'"
            )
        return

    if resolved_policy == policy_impute_minprob:
        if min_observed_values is not None:
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.min_observed_values must be "
                "None when missing_data.policy='impute_minprob'"
            )
        if isinstance(q, bool) or not isinstance(q, (int, float)):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.q must "
                "be a float when missing_data.policy='impute_minprob'"
            )
        q_value = float(q)
        if not math.isfinite(q_value) or not (0.0 < q_value < 0.5):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.q must "
                "satisfy 0 < q < 0.5 when missing_data.policy='impute_minprob'"
            )
        if isinstance(width, bool) or not isinstance(width, (int, float)):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.width "
                "must be a float when missing_data.policy='impute_minprob'"
            )
        width_value = float(width)
        if not math.isfinite(width_value) or not (0.0 < width_value <= 1.0):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.width "
                "must satisfy 0 < width <= 1.0 when "
                "missing_data.policy='impute_minprob'"
            )
        if isinstance(seed, bool) or not isinstance(seed, int):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.seed "
                "must be an int when missing_data.policy='impute_minprob'"
            )
        if seed < 0:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.seed "
                "must be greater than or equal to 0 when "
                "missing_data.policy='impute_minprob'"
            )
        if isinstance(max_missing_fraction_per_row, bool) or not isinstance(
            max_missing_fraction_per_row, (int, float)
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data."
                "max_missing_fraction_per_row must be a float when "
                "missing_data.policy='impute_minprob'"
            )
        max_missing_fraction_value = float(max_missing_fraction_per_row)
        if not math.isfinite(max_missing_fraction_value) or not (
            0.0 < max_missing_fraction_value <= 1.0
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data."
                "max_missing_fraction_per_row must satisfy 0 < value <= 1 when "
                "missing_data.policy='impute_minprob'"
            )
        if k is not None or distance is not None:
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.k and .distance "
                "must be None when missing_data.policy='impute_minprob'"
            )
        return

    if resolved_policy == policy_impute_knn:
        if min_observed_values is not None:
            raise PhosPyInputError(
                "dataset build request "
                "preprocessing_config.missing_data.min_observed_values must be "
                "None when missing_data.policy='impute_knn'"
            )
        if q is not None or width is not None or seed is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.q, "
                ".width, and .seed must be None when "
                "missing_data.policy='impute_knn'"
            )
        if isinstance(k, bool) or not isinstance(k, int):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.k "
                "must be an int when missing_data.policy='impute_knn'"
            )
        if k < 1:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.k "
                "must be greater than or equal to 1 when "
                "missing_data.policy='impute_knn'"
            )
        if not isinstance(distance, str) or distance.strip() != "nan_euclidean":
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.distance "
                "must be 'nan_euclidean' when missing_data.policy='impute_knn'"
            )
        if isinstance(max_missing_fraction_per_row, bool) or not isinstance(
            max_missing_fraction_per_row, (int, float)
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data."
                "max_missing_fraction_per_row must be a float when "
                "missing_data.policy='impute_knn'"
            )
        max_missing_fraction_value = float(max_missing_fraction_per_row)
        if not math.isfinite(max_missing_fraction_value) or not (
            0.0 < max_missing_fraction_value <= 1.0
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data."
                "max_missing_fraction_per_row must satisfy 0 < value <= 1 when "
                "missing_data.policy='impute_knn'"
            )
        return

    raise PhosPyInputError(
        "dataset build request preprocessing_config contains an unsupported "
        "missing_data.policy"
    )


def validate_site_matrix_config(
    *,
    policy: object,
    duplicate_site_policy: object,
    missing_data_policy: object,
    minimum_observed_values: object | None,
    supported_policies: Collection[str],
    supported_duplicate_policies: Collection[str],
    supported_missing_data_policies: Collection[str],
    policy_as_input: str,
    duplicate_policy_default: str,
    supported_missing_data_policy: str,
) -> None:
    """Validate public site-matrix config fields."""

    resolved_policy = require_supported_literal(
        policy,
        field_name="dataset build request preprocessing_config.site_matrix.policy",
        supported_values=supported_policies,
        error_type=PhosPyInputError,
    )
    resolved_duplicate_policy = require_supported_literal(
        duplicate_site_policy,
        field_name=(
            "dataset build request preprocessing_config.site_matrix."
            "duplicate_site_policy"
        ),
        supported_values=supported_duplicate_policies,
        error_type=PhosPyInputError,
    )
    resolved_missing_data_policy = str(missing_data_policy)
    if resolved_missing_data_policy not in supported_missing_data_policies:
        if (
            resolved_missing_data_policy
            in _INCOMPATIBLE_SITE_MATRIX_MISSING_DATA_POLICIES
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                f"missing_data_policy='{resolved_missing_data_policy}' is not "
                "supported for strict AnalysisReadyPhosphoDataset construction in "
                "the public complete-case builder lane. Use "
                "site_matrix.missing_data_policy="
                f"'{supported_missing_data_policy}'."
            )
        require_supported_literal(
            missing_data_policy,
            field_name=(
                "dataset build request preprocessing_config.site_matrix."
                "missing_data_policy"
            ),
            supported_values=supported_missing_data_policies,
            error_type=PhosPyInputError,
        )
    if minimum_observed_values is not None:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.site_matrix."
            "minimum_observed_values is not supported for strict "
            "AnalysisReadyPhosphoDataset construction and must be None"
        )

    if (
        resolved_policy == policy_as_input
        and resolved_duplicate_policy != duplicate_policy_default
    ):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.site_matrix."
            "duplicate_site_policy is only valid when "
            "site_matrix.policy='build_from_metadata'"
        )


def validate_comparison_building_config(
    *,
    policy: object,
    sample_group_column: object,
    pairs: object | None,
    supported_policies: Collection[str],
    policy_none: str,
    policy_sample_metadata_pairs: str,
) -> None:
    """Validate public comparison-building config fields."""

    resolved_policy = require_supported_literal(
        policy,
        field_name="dataset build request preprocessing_config.comparisons.policy",
        supported_values=supported_policies,
        error_type=PhosPyInputError,
    )
    require_non_empty_string(
        sample_group_column,
        field_name=(
            "dataset build request preprocessing_config.comparisons.sample_group_column"
        ),
        error_type=PhosPyInputError,
    )
    if resolved_policy == policy_none:
        if pairs is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs "
                "must be None when comparisons.policy='none'"
            )
        return
    if resolved_policy != policy_sample_metadata_pairs:
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "comparisons.policy"
        )
    if pairs is None:
        return
    if not isinstance(pairs, (tuple, list)):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.comparisons.pairs must be "
            "a sequence of (left_group, right_group) pairs when provided"
        )
    resolved_pairs = tuple(pairs)
    if not resolved_pairs:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.comparisons.pairs must "
            "contain at least one pair when provided"
        )
    _validate_comparison_pairs(resolved_pairs)


def validate_site_sequence_resolution_config(
    *,
    mode: object,
    conflict_policy: object | None,
    flank_size: object,
    accession_column: object,
    site_column: object,
    fasta_path: object | None,
    supported_modes: Collection[str],
    supported_conflict_policies: Collection[str],
) -> None:
    """Validate public site-sequence-resolution config fields."""

    require_supported_literal(
        mode,
        field_name=(
            "dataset build request preprocessing_config.site_sequence_resolution.mode"
        ),
        supported_values=supported_modes,
        error_type=PhosPyInputError,
    )
    if conflict_policy is not None:
        require_supported_literal(
            conflict_policy,
            field_name=(
                "dataset build request preprocessing_config.site_sequence_resolution."
                "conflict_policy"
            ),
            supported_values=supported_conflict_policies,
            error_type=PhosPyInputError,
        )

    require_non_empty_string(
        accession_column,
        field_name=(
            "dataset build request preprocessing_config.site_sequence_resolution."
            "accession_column"
        ),
        error_type=PhosPyInputError,
    )
    require_non_empty_string(
        site_column,
        field_name=(
            "dataset build request preprocessing_config.site_sequence_resolution."
            "site_column"
        ),
        error_type=PhosPyInputError,
    )
    if fasta_path is None:
        return
    require_local_filesystem_path(
        fasta_path,
        field_name=(
            "dataset build request preprocessing_config.site_sequence_resolution."
            "fasta_path"
        ),
        error_type=PhosPyInputError,
        when_provided=True,
    )
    if isinstance(flank_size, bool) or not isinstance(flank_size, int):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.site_sequence_resolution."
            "flank_size must be an int when fasta_path is provided"
        )
    if flank_size < 1:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.site_sequence_resolution."
            "flank_size must be greater than or equal to 1 when fasta_path is "
            "provided"
        )


def validate_total_protein_correction_identity_config(
    *,
    mode: object,
    phosphosite_key: object,
    total_protein_key: object,
    mapping_table: object | None,
    mapping_phosphosite_key: object | None,
    mapping_total_protein_key: object | None,
    matching_policy: object,
    duplicate_policy: object,
    unmatched_policy: object,
    supported_modes: Collection[str],
    supported_matching_policies: Collection[str],
    supported_duplicate_policies: Collection[str],
    supported_unmatched_policies: Collection[str],
    mode_direct: str,
    mode_mapping_table: str,
    matching_policy_gene_symbol_normalised: str,
) -> None:
    """Validate public total-protein identity config fields."""

    resolved_mode = require_supported_literal(
        mode,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.mode"
        ),
        supported_values=supported_modes,
        error_type=PhosPyInputError,
    )
    require_non_empty_string(
        phosphosite_key,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.phosphosite_key"
        ),
        error_type=PhosPyInputError,
    )
    require_non_empty_string(
        total_protein_key,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.total_protein_key"
        ),
        error_type=PhosPyInputError,
    )
    resolved_matching_policy = require_supported_literal(
        matching_policy,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.matching_policy"
        ),
        supported_values=supported_matching_policies,
        error_type=PhosPyInputError,
    )
    require_supported_literal(
        duplicate_policy,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.duplicate_policy"
        ),
        supported_values=supported_duplicate_policies,
        error_type=PhosPyInputError,
    )
    require_supported_literal(
        unmatched_policy,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.unmatched_policy"
        ),
        supported_values=supported_unmatched_policies,
        error_type=PhosPyInputError,
    )

    uses_gene_symbol_keys = "gene_symbol" in {
        str(phosphosite_key).strip().lower(),
        str(total_protein_key).strip().lower(),
        (
            ""
            if mapping_phosphosite_key is None
            else str(mapping_phosphosite_key).strip().lower()
        ),
        (
            ""
            if mapping_total_protein_key is None
            else str(mapping_total_protein_key).strip().lower()
        ),
    }
    if (
        resolved_matching_policy == matching_policy_gene_symbol_normalised
        and not uses_gene_symbol_keys
    ):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.matching_policy='gene_symbol_normalised' requires at "
            "least one gene_symbol identity key "
            "(phosphosite_key/total_protein_key/mapping keys)"
        )

    if resolved_mode == mode_direct:
        if mapping_table is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.mapping_table must be None when identity.mode='direct'"
            )
        if mapping_phosphosite_key is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.mapping_phosphosite_key must be None when "
                "identity.mode='direct'"
            )
        if mapping_total_protein_key is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                "identity.mapping_total_protein_key must be None when "
                "identity.mode='direct'"
            )
        return

    if resolved_mode != mode_mapping_table:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity contains an unsupported mode"
        )
    if mapping_table is None:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.mapping_table is required when identity.mode='mapping_table'"
        )
    if not isinstance(mapping_table, pd.DataFrame):
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.mapping_table must be a pandas DataFrame"
        )
    require_non_empty_string(
        mapping_phosphosite_key,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.mapping_phosphosite_key"
        ),
        error_type=PhosPyInputError,
        when_provided=True,
    )
    require_non_empty_string(
        mapping_total_protein_key,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.mapping_total_protein_key"
        ),
        error_type=PhosPyInputError,
        when_provided=True,
    )


def validate_total_protein_correction_config(
    *,
    policy: object,
    identity: object,
    supported_policies: Collection[str],
    identity_type: type[object],
) -> None:
    """Validate public total-protein-correction config fields."""

    require_supported_literal(
        policy,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction.policy"
        ),
        supported_values=supported_policies,
        error_type=PhosPyInputError,
    )
    require_instance(
        identity,
        expected_type=identity_type,
        field_name=(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity"
        ),
        error_type=PhosPyInputError,
    )


def _validate_comparison_pairs(resolved_pairs: Sequence[object]) -> None:
    seen_pairs: set[tuple[str, str]] = set()
    for pair in resolved_pairs:
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs "
                "must contain only (left_group, right_group) tuples"
            )
        left_group, right_group = pair
        if not isinstance(left_group, str) or not left_group.strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs "
                "must contain non-empty left_group strings"
            )
        if not isinstance(right_group, str) or not right_group.strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs "
                "must contain non-empty right_group strings"
            )
        left = left_group.strip()
        right = right_group.strip()
        if left == right:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs "
                "cannot contain self-comparison pairs"
            )
        canonical_pair = (left, right) if left <= right else (right, left)
        if canonical_pair in seen_pairs:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs "
                "contains duplicate pairs regardless of direction"
            )
        seen_pairs.add(canonical_pair)


def _coerce_internal_enum(
    enum_type: type[Enum],
    value: object,
    *,
    field_name: str,
) -> Enum:
    return coerce_policy_enum(
        enum_type,
        value,
        field_name=field_name,
        error_type=PhosPyInputError,
    )


def _require_non_empty_string_sequence(
    value: object,
    *,
    field_name: str,
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
