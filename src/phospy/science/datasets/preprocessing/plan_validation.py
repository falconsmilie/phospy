"""Validation for directly constructed preprocessing plans."""

from __future__ import annotations

import math
from typing import Any

from phospy.errors.input import PhosPyInputError
from phospy.science.configs.preprocessing import (
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICIES,
    DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_COLUMN_MEAN_WITH_CAVEAT,
)
from phospy.science.datasets.preprocessing.imputation_scale_policy import (
    reject_incompatible_imputation_stage_order,
    resolve_imputation_scale_policy,
)
from phospy.science.datasets.preprocessing.plan_rules import (
    PreprocessingBatchCorrectionPlanRuleFamily,
    PreprocessingCorePlanPolicyRuleFamily,
    PreprocessingGroupCoveragePlanRuleFamily,
    PreprocessingLocalisationPlanRuleFamily,
    PreprocessingRuvReadinessPlanRuleFamily,
    PreprocessingSiteMatrixComparisonPlanRuleFamily,
    PreprocessingSiteSequencePlanRuleFamily,
    PreprocessingTotalProteinCorrectionPlanRuleFamily,
)
from phospy.science.datasets.preprocessing.plan_stage_order import (
    PreprocessingStageOrderValidator,
    normalize_stage_order_resolution,
)

_MutablePreprocessingPlan = Any


def validate_resolved_preprocessing_plan(plan: _MutablePreprocessingPlan) -> None:
    """Normalize and validate directly supplied resolved plan fields."""

    _validate_core_transform_fields(plan)
    _validate_imputation_fields(plan)
    _validate_localisation_fields(plan)
    _validate_site_sequence_fields(plan)
    _validate_total_protein_fields(plan)
    _validate_site_matrix_comparison_fields(plan)
    _validate_ruv_readiness_fields(plan)
    _validate_batch_correction_fields(plan)
    _validate_stage_order_fields(plan)
    _validate_group_coverage_fields(plan)


def _validate_core_transform_fields(plan: _MutablePreprocessingPlan) -> None:
    core = PreprocessingCorePlanPolicyRuleFamily().run(
        intensity_transform_policy=plan.intensity_transform_policy,
        normalisation_policy=plan.normalisation_policy,
        missing_data_policy=plan.missing_data_policy,
    )
    _set(plan, "intensity_transform_policy", core.intensity_transform_policy)
    _set(plan, "normalisation_policy", core.normalisation_policy)
    _set(plan, "missing_data_policy", core.missing_data_policy)
    pseudocount = float(plan.intensity_transform_pseudocount)
    if not math.isfinite(pseudocount) or pseudocount < 0:
        raise PhosPyInputError(
            "dataset preprocessing plan intensity_transform_pseudocount "
            "(internal model) must be finite and greater than or equal to 0"
        )
    _set(plan, "intensity_transform_pseudocount", pseudocount)


def _validate_imputation_fields(plan: _MutablePreprocessingPlan) -> None:
    imputation_scale = resolve_imputation_scale_policy(
        missing_data_policy=plan.missing_data_policy,
        requested_input_scale=plan.missing_data_input_scale,
        intensity_transform_policy=plan.intensity_transform_policy,
    )
    _reject_inconsistent_resolved_imputation_field(
        field_name="missing_data_input_scale_source",
        supplied=plan.missing_data_input_scale_source,
        expected=imputation_scale.input_scale_source,
    )
    _reject_inconsistent_resolved_imputation_field(
        field_name="missing_data_imputation_operation_order",
        supplied=plan.missing_data_imputation_operation_order,
        expected=imputation_scale.operation_order,
    )
    _set(plan, "missing_data_input_scale", imputation_scale.input_scale)
    _set(
        plan,
        "missing_data_input_scale_source",
        imputation_scale.input_scale_source,
    )
    _set(
        plan,
        "missing_data_imputation_operation_order",
        imputation_scale.operation_order,
    )
    reject_incompatible_imputation_stage_order(
        stage_order=plan.stage_order,
        resolved_policy=imputation_scale,
    )
    _validate_missing_data_no_overlap_policy(plan)


def _validate_missing_data_no_overlap_policy(plan: _MutablePreprocessingPlan) -> None:
    policy_value = plan.missing_data_policy.value
    no_overlap_policy = plan.missing_data_no_overlap_policy
    if policy_value != "impute_knn":
        if no_overlap_policy is not None:
            raise PhosPyInputError(
                "dataset preprocessing plan missing_data_no_overlap_policy "
                "(internal model) must be None unless "
                "missing_data_policy='impute_knn'"
            )
        return
    if no_overlap_policy is None:
        _set(
            plan,
            "missing_data_no_overlap_policy",
            DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICY_COLUMN_MEAN_WITH_CAVEAT,
        )
        return
    normalized = str(no_overlap_policy).strip()
    if normalized not in DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICIES:
        supported = ", ".join(sorted(DATASET_MISSING_DATA_KNN_NO_OVERLAP_POLICIES))
        raise PhosPyInputError(
            "dataset preprocessing plan missing_data_no_overlap_policy "
            f"(internal model) must be one of: {supported}"
        )
    _set(plan, "missing_data_no_overlap_policy", normalized)


def _reject_inconsistent_resolved_imputation_field(
    *,
    field_name: str,
    supplied: object | None,
    expected: object | None,
) -> None:
    if supplied is None or supplied == expected:
        return
    raise PhosPyInputError(
        f"dataset preprocessing plan {field_name} (internal model) is "
        "inconsistent with the resolved missing-data input-scale policy"
    )


def _validate_localisation_fields(plan: _MutablePreprocessingPlan) -> None:
    localisation = PreprocessingLocalisationPlanRuleFamily().run(
        localisation_mode=plan.localisation_mode,
        localisation_min_confidence=plan.localisation_min_confidence,
        localisation_confidence_column=plan.localisation_confidence_column,
        localisation_waiver_reason=plan.localisation_waiver_reason,
    )
    _set(plan, "localisation_mode", localisation.localisation_mode)
    _set(
        plan,
        "localisation_min_confidence",
        localisation.localisation_min_confidence,
    )
    _set(
        plan,
        "localisation_confidence_column",
        localisation.localisation_confidence_column,
    )
    _set(
        plan,
        "localisation_waiver_reason",
        localisation.localisation_waiver_reason,
    )


def _validate_site_sequence_fields(plan: _MutablePreprocessingPlan) -> None:
    site_sequence = PreprocessingSiteSequencePlanRuleFamily().run(
        site_sequence_resolution_enabled=plan.site_sequence_resolution_enabled,
        site_sequence_resolution_fasta_path=plan.site_sequence_resolution_fasta_path,
        site_sequence_resolution_mode=plan.site_sequence_resolution_mode,
        site_sequence_resolution_conflict_policy=(
            plan.site_sequence_resolution_conflict_policy
        ),
        site_sequence_resolution_flank_size=plan.site_sequence_resolution_flank_size,
        site_sequence_resolution_accession_column=(
            plan.site_sequence_resolution_accession_column
        ),
        site_sequence_resolution_site_column=plan.site_sequence_resolution_site_column,
    )
    _set(
        plan,
        "site_sequence_resolution_enabled",
        site_sequence.site_sequence_resolution_enabled,
    )
    _set(
        plan,
        "site_sequence_resolution_fasta_path",
        site_sequence.site_sequence_resolution_fasta_path,
    )
    _set(
        plan,
        "site_sequence_resolution_mode",
        site_sequence.site_sequence_resolution_mode,
    )
    _set(
        plan,
        "site_sequence_resolution_conflict_policy",
        site_sequence.site_sequence_resolution_conflict_policy,
    )
    _set(
        plan,
        "site_sequence_resolution_flank_size",
        site_sequence.site_sequence_resolution_flank_size,
    )
    _set(
        plan,
        "site_sequence_resolution_accession_column",
        site_sequence.site_sequence_resolution_accession_column,
    )
    _set(
        plan,
        "site_sequence_resolution_site_column",
        site_sequence.site_sequence_resolution_site_column,
    )


def _validate_total_protein_fields(plan: _MutablePreprocessingPlan) -> None:
    total_protein = PreprocessingTotalProteinCorrectionPlanRuleFamily().run(
        total_protein_correction_policy=plan.total_protein_correction_policy,
        total_protein_correction_identity_policy=(
            plan.total_protein_correction_identity_policy
        ),
        protein_aware_preparation_policy=plan.protein_aware_preparation_policy,
        protein_aware_preparation_mapping_policy=(
            plan.protein_aware_preparation_mapping_policy
        ),
    )
    _set(
        plan,
        "total_protein_correction_policy",
        total_protein.total_protein_correction_policy,
    )
    _set(
        plan,
        "total_protein_correction_identity_policy",
        total_protein.total_protein_correction_identity_policy,
    )
    _set(
        plan,
        "protein_aware_preparation_policy",
        total_protein.protein_aware_preparation_policy,
    )
    _set(
        plan,
        "protein_aware_preparation_mapping_policy",
        total_protein.protein_aware_preparation_mapping_policy,
    )


def _validate_site_matrix_comparison_fields(
    plan: _MutablePreprocessingPlan,
) -> None:
    site_matrix_comparison = PreprocessingSiteMatrixComparisonPlanRuleFamily().run(
        site_matrix_policy=plan.site_matrix_policy,
        site_matrix_duplicate_site_policy=plan.site_matrix_duplicate_site_policy,
        site_matrix_missing_data_policy=plan.site_matrix_missing_data_policy,
        site_matrix_minimum_observed_values=plan.site_matrix_minimum_observed_values,
        comparison_building_policy=plan.comparison_building_policy,
        comparison_sample_group_column=plan.comparison_sample_group_column,
        comparison_pairs=plan.comparison_pairs,
    )
    _set(plan, "site_matrix_policy", site_matrix_comparison.site_matrix_policy)
    _set(
        plan,
        "site_matrix_duplicate_site_policy",
        site_matrix_comparison.site_matrix_duplicate_site_policy,
    )
    _set(
        plan,
        "site_matrix_missing_data_policy",
        site_matrix_comparison.site_matrix_missing_data_policy,
    )
    _set(
        plan,
        "site_matrix_minimum_observed_values",
        site_matrix_comparison.site_matrix_minimum_observed_values,
    )
    _set(
        plan,
        "comparison_building_policy",
        site_matrix_comparison.comparison_building_policy,
    )
    _set(
        plan,
        "comparison_sample_group_column",
        site_matrix_comparison.comparison_sample_group_column,
    )
    _set(plan, "comparison_pairs", site_matrix_comparison.comparison_pairs)


def _validate_ruv_readiness_fields(plan: _MutablePreprocessingPlan) -> None:
    ruv_readiness = PreprocessingRuvReadinessPlanRuleFamily().run(
        ruv_readiness_enabled=plan.ruv_readiness_enabled,
        ruv_readiness_control_feature_column=(
            plan.ruv_readiness_control_feature_column
        ),
        ruv_readiness_replicate_group_column=(
            plan.ruv_readiness_replicate_group_column
        ),
        ruv_readiness_batch_column=plan.ruv_readiness_batch_column,
    )
    _set(plan, "ruv_readiness_enabled", ruv_readiness.ruv_readiness_enabled)
    _set(
        plan,
        "ruv_readiness_control_feature_column",
        ruv_readiness.ruv_readiness_control_feature_column,
    )
    _set(
        plan,
        "ruv_readiness_replicate_group_column",
        ruv_readiness.ruv_readiness_replicate_group_column,
    )
    _set(
        plan,
        "ruv_readiness_batch_column",
        ruv_readiness.ruv_readiness_batch_column,
    )


def _validate_batch_correction_fields(plan: _MutablePreprocessingPlan) -> None:
    batch_correction = PreprocessingBatchCorrectionPlanRuleFamily().run(
        batch_correction_method=plan.batch_correction_method,
        batch_correction_batch_column=plan.batch_correction_batch_column,
        batch_correction_condition_column=plan.batch_correction_condition_column,
        batch_correction_condition_columns=plan.batch_correction_condition_columns,
        batch_correction_replicate_column=plan.batch_correction_replicate_column,
        batch_correction_control_site_set=plan.batch_correction_control_site_set,
        batch_correction_missingness_policy=(plan.batch_correction_missingness_policy),
        batch_correction_internal_request=plan.batch_correction_internal_request,
        batch_correction_preserve_condition_effects=(
            plan.batch_correction_preserve_condition_effects
        ),
    )
    _set(
        plan,
        "batch_correction_method",
        batch_correction.batch_correction_method,
    )
    _set(
        plan,
        "batch_correction_batch_column",
        batch_correction.batch_correction_batch_column,
    )
    _set(
        plan,
        "batch_correction_condition_column",
        batch_correction.batch_correction_condition_column,
    )
    _set(
        plan,
        "batch_correction_condition_columns",
        batch_correction.batch_correction_condition_columns,
    )
    _set(
        plan,
        "batch_correction_replicate_column",
        batch_correction.batch_correction_replicate_column,
    )
    _set(
        plan,
        "batch_correction_control_site_set",
        batch_correction.batch_correction_control_site_set,
    )
    _set(
        plan,
        "batch_correction_missingness_policy",
        batch_correction.batch_correction_missingness_policy,
    )
    _set(
        plan,
        "batch_correction_internal_request",
        batch_correction.batch_correction_internal_request,
    )
    _set(
        plan,
        "batch_correction_preserve_condition_effects",
        batch_correction.batch_correction_preserve_condition_effects,
    )


def _validate_stage_order_fields(plan: _MutablePreprocessingPlan) -> None:
    stage_order = tuple(str(stage).strip() for stage in plan.stage_order)
    _set(plan, "stage_order", stage_order)
    PreprocessingStageOrderValidator().run(
        stage_order=stage_order,
        batch_correction_requested=(
            plan.batch_correction_method != DATASET_BATCH_CORRECTION_METHOD_NONE
        ),
    )
    _set(
        plan,
        "stage_order_resolution",
        normalize_stage_order_resolution(
            stage_order=stage_order,
            stage_order_resolution=plan.stage_order_resolution,
        ),
    )


def _validate_group_coverage_fields(plan: _MutablePreprocessingPlan) -> None:
    group_coverage = PreprocessingGroupCoveragePlanRuleFamily().run(
        enabled=plan.group_coverage_filter_enabled,
        group_column=plan.group_coverage_filter_group_column,
        min_finite_observations_per_group=(
            plan.group_coverage_filter_min_finite_observations_per_group
        ),
        min_finite_fraction_per_group=(
            plan.group_coverage_filter_min_finite_fraction_per_group
        ),
        min_groups_passing_threshold=(
            plan.group_coverage_filter_min_groups_passing_threshold
        ),
        stage_order=plan.stage_order,
    )
    _set(
        plan,
        "group_coverage_filter_enabled",
        group_coverage.group_coverage_filter_enabled,
    )
    _set(
        plan,
        "group_coverage_filter_group_column",
        group_coverage.group_coverage_filter_group_column,
    )
    _set(
        plan,
        "group_coverage_filter_min_finite_observations_per_group",
        group_coverage.group_coverage_filter_min_finite_observations_per_group,
    )
    _set(
        plan,
        "group_coverage_filter_min_finite_fraction_per_group",
        group_coverage.group_coverage_filter_min_finite_fraction_per_group,
    )
    _set(
        plan,
        "group_coverage_filter_min_groups_passing_threshold",
        group_coverage.group_coverage_filter_min_groups_passing_threshold,
    )


def _set(plan: _MutablePreprocessingPlan, field_name: str, value: object) -> None:
    object.__setattr__(plan, field_name, value)


__all__ = ["validate_resolved_preprocessing_plan"]
