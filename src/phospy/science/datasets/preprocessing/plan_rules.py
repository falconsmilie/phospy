"""Internal rule families for execution-ready preprocessing plans."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.science.configs.preprocessing import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED,
    SPS_RUV_BATCH_CORRECTION_METHODS,
    CorrectionMissingnessPolicy,
    DatasetProteinAwarePreparationMappingPolicy,
    DatasetProteinAwarePreparationPolicy,
    InternalBatchCorrectionRequest,
)
from phospy.science.configs.preprocessing._validation import (
    reject_unsupported_ruv_iii_style_method,
    validate_group_coverage_filter_config,
)
from phospy.science.datasets.preprocessing.plan_constants import (
    DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    UNSUPPORTED_BATCH_CORRECTION_METHOD_RUV_III_STYLE,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    IntensityTransformPolicy,
    LocalisationEligibilityMode,
    MissingDataPolicy,
    NormalisationPolicy,
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
    TotalProteinCorrectionPolicy,
)


@dataclass(frozen=True, slots=True)
class ResolvedCorePreprocessingPlanPolicies:
    """Core quantitative preprocessing policies resolved for internal execution."""

    intensity_transform_policy: IntensityTransformPolicy
    normalisation_policy: NormalisationPolicy
    missing_data_policy: MissingDataPolicy


@dataclass(frozen=True, slots=True)
class ResolvedLocalisationPlanFields:
    """Localisation evidence policy resolved for internal execution."""

    localisation_mode: LocalisationEligibilityMode
    localisation_min_confidence: float
    localisation_confidence_column: str
    localisation_waiver_reason: str | None


@dataclass(frozen=True, slots=True)
class ResolvedDownstreamPreprocessingPlanPolicies:
    """Site, comparison, total-protein, and protein-aware plan policies."""

    site_sequence_resolution_mode: SiteSequenceResolutionMode
    site_sequence_resolution_conflict_policy: SiteSequenceConflictPolicy
    site_matrix_policy: SiteMatrixPolicy
    comparison_building_policy: ComparisonBuildingPolicy
    site_matrix_duplicate_site_policy: SiteMatrixDuplicateSitePolicy
    site_matrix_missing_data_policy: SiteMatrixMissingDataPolicy
    total_protein_correction_policy: TotalProteinCorrectionPolicy
    protein_aware_preparation_policy: DatasetProteinAwarePreparationPolicy
    protein_aware_preparation_mapping_policy: (
        DatasetProteinAwarePreparationMappingPolicy
    )


@dataclass(frozen=True, slots=True)
class ResolvedBatchCorrectionPlanFields:
    """Batch-correction plan fields resolved for internal execution."""

    batch_correction_method: str
    batch_correction_batch_column: str
    batch_correction_condition_column: str
    batch_correction_condition_columns: tuple[str, ...]
    batch_correction_replicate_column: str | None


class PreprocessingCorePlanPolicyRuleFamily:
    """Resolve core quantitative policy fields on an internal plan."""

    def run(
        self,
        *,
        intensity_transform_policy: IntensityTransformPolicy,
        normalisation_policy: NormalisationPolicy,
        missing_data_policy: MissingDataPolicy,
    ) -> ResolvedCorePreprocessingPlanPolicies:
        return ResolvedCorePreprocessingPlanPolicies(
            intensity_transform_policy=IntensityTransformPolicy.parse(
                intensity_transform_policy,
                field_name=(
                    "dataset preprocessing plan intensity_transform_policy "
                    "(internal model)"
                ),
            ),
            normalisation_policy=NormalisationPolicy.parse(
                normalisation_policy,
                field_name=(
                    "dataset preprocessing plan normalisation_policy (internal model)"
                ),
            ),
            missing_data_policy=MissingDataPolicy.parse(
                missing_data_policy,
                field_name=(
                    "dataset preprocessing plan missing_data_policy (internal model)"
                ),
            ),
        )


class PreprocessingLocalisationPlanRuleFamily:
    """Resolve and validate localisation evidence policy fields."""

    def run(
        self,
        *,
        localisation_mode: LocalisationEligibilityMode,
        localisation_min_confidence: float,
        localisation_confidence_column: str,
        localisation_waiver_reason: str | None,
    ) -> ResolvedLocalisationPlanFields:
        resolved_mode = LocalisationEligibilityMode.parse(
            localisation_mode,
            field_name="dataset preprocessing plan localisation_mode (internal model)",
        )
        min_confidence = float(localisation_min_confidence)
        if not math.isfinite(min_confidence) or not (0.0 <= min_confidence <= 1.0):
            raise PhosPyInputError(
                "dataset preprocessing plan localisation_min_confidence "
                "(internal model) must be between 0.0 and 1.0"
            )
        confidence_column = str(localisation_confidence_column).strip()
        if confidence_column == "":
            raise PhosPyInputError(
                "dataset preprocessing plan localisation_confidence_column "
                "(internal model) must be a non-empty string"
            )
        if (
            resolved_mode is not LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER
            and localisation_waiver_reason is not None
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan localisation_waiver_reason "
                "(internal model) must be None unless localisation_mode="
                "'allow_missing_with_waiver'"
            )
        if resolved_mode is LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER:
            waiver_reason = (
                "" if localisation_waiver_reason is None else localisation_waiver_reason
            ).strip()
            if waiver_reason == "":
                raise PhosPyInputError(
                    "dataset preprocessing plan localisation_waiver_reason "
                    "(internal model) must be provided when localisation_mode="
                    "'allow_missing_with_waiver'"
                )
            return ResolvedLocalisationPlanFields(
                localisation_mode=resolved_mode,
                localisation_min_confidence=min_confidence,
                localisation_confidence_column=confidence_column,
                localisation_waiver_reason=waiver_reason,
            )
        return ResolvedLocalisationPlanFields(
            localisation_mode=resolved_mode,
            localisation_min_confidence=min_confidence,
            localisation_confidence_column=confidence_column,
            localisation_waiver_reason=localisation_waiver_reason,
        )


class PreprocessingDownstreamPlanPolicyRuleFamily:
    """Resolve downstream site, comparison, and total-protein policy fields."""

    def run(
        self,
        *,
        site_sequence_resolution_mode: SiteSequenceResolutionMode,
        site_sequence_resolution_conflict_policy: SiteSequenceConflictPolicy,
        site_matrix_policy: SiteMatrixPolicy,
        comparison_building_policy: ComparisonBuildingPolicy,
        site_matrix_duplicate_site_policy: SiteMatrixDuplicateSitePolicy,
        site_matrix_missing_data_policy: SiteMatrixMissingDataPolicy,
        total_protein_correction_policy: TotalProteinCorrectionPolicy,
        protein_aware_preparation_policy: DatasetProteinAwarePreparationPolicy = (
            DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED
        ),
        protein_aware_preparation_mapping_policy: (
            DatasetProteinAwarePreparationMappingPolicy
        ) = DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
    ) -> ResolvedDownstreamPreprocessingPlanPolicies:
        return ResolvedDownstreamPreprocessingPlanPolicies(
            site_sequence_resolution_mode=SiteSequenceResolutionMode.parse(
                site_sequence_resolution_mode,
                field_name=(
                    "dataset preprocessing plan site_sequence_resolution_mode "
                    "(internal model)"
                ),
            ),
            site_sequence_resolution_conflict_policy=SiteSequenceConflictPolicy.parse(
                site_sequence_resolution_conflict_policy,
                field_name=(
                    "dataset preprocessing plan "
                    "site_sequence_resolution_conflict_policy (internal model)"
                ),
            ),
            site_matrix_policy=SiteMatrixPolicy.parse(
                site_matrix_policy,
                field_name=(
                    "dataset preprocessing plan site_matrix_policy (internal model)"
                ),
            ),
            comparison_building_policy=ComparisonBuildingPolicy.parse(
                comparison_building_policy,
                field_name=(
                    "dataset preprocessing plan comparison_building_policy "
                    "(internal model)"
                ),
            ),
            site_matrix_duplicate_site_policy=SiteMatrixDuplicateSitePolicy.parse(
                site_matrix_duplicate_site_policy,
                field_name=(
                    "dataset preprocessing plan site_matrix_duplicate_site_policy "
                    "(internal model)"
                ),
            ),
            site_matrix_missing_data_policy=SiteMatrixMissingDataPolicy.parse(
                site_matrix_missing_data_policy,
                field_name=(
                    "dataset preprocessing plan site_matrix_missing_data_policy "
                    "(internal model)"
                ),
            ),
            total_protein_correction_policy=TotalProteinCorrectionPolicy.parse(
                total_protein_correction_policy,
                field_name=(
                    "dataset preprocessing plan total_protein_correction_policy "
                    "(internal model)"
                ),
            ),
            protein_aware_preparation_policy=cast(
                DatasetProteinAwarePreparationPolicy,
                str(protein_aware_preparation_policy).strip(),
            ),
            protein_aware_preparation_mapping_policy=cast(
                DatasetProteinAwarePreparationMappingPolicy,
                str(protein_aware_preparation_mapping_policy).strip(),
            ),
        )


class PreprocessingBatchCorrectionPlanRuleFamily:
    """Resolve and validate internal batch-correction plan fields."""

    def run(
        self,
        *,
        batch_correction_method: str,
        batch_correction_batch_column: str,
        batch_correction_condition_column: str,
        batch_correction_condition_columns: tuple[str, ...],
        batch_correction_replicate_column: str | None,
        batch_correction_control_site_set: object | None,
        batch_correction_missingness_policy: CorrectionMissingnessPolicy | None,
        batch_correction_internal_request: InternalBatchCorrectionRequest | None,
        batch_correction_preserve_condition_effects: bool,
    ) -> ResolvedBatchCorrectionPlanFields:
        method = str(batch_correction_method).strip()
        if not method:
            method = DATASET_BATCH_CORRECTION_METHOD_NONE
        if method == UNSUPPORTED_BATCH_CORRECTION_METHOD_RUV_III_STYLE:
            reject_unsupported_ruv_iii_style_method(
                method,
                field_name=(
                    "dataset preprocessing plan batch_correction_method "
                    "(internal model)"
                ),
            )
        if (
            method
            not in {
                DATASET_BATCH_CORRECTION_METHOD_NONE,
                DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
            }
            and method not in SPS_RUV_BATCH_CORRECTION_METHODS
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan batch_correction_method "
                "(internal model) must be one of: none, "
                "linear_residualize_batch, sps_ruv_style"
            )
        batch_column = str(batch_correction_batch_column).strip()
        condition_column = str(batch_correction_condition_column).strip()
        if batch_column == "":
            raise PhosPyInputError(
                "dataset preprocessing plan batch_correction_batch_column "
                "(internal model) must be a non-empty string"
            )
        if condition_column == "":
            raise PhosPyInputError(
                "dataset preprocessing plan batch_correction_condition_column "
                "(internal model) must be a non-empty string"
            )
        condition_columns = tuple(
            str(column).strip() for column in batch_correction_condition_columns
        )
        if not condition_columns or any(column == "" for column in condition_columns):
            raise PhosPyInputError(
                "dataset preprocessing plan batch_correction_condition_columns "
                "(internal model) must be a non-empty tuple of non-empty strings"
            )
        if len(set(condition_columns)) != len(condition_columns):
            raise PhosPyInputError(
                "dataset preprocessing plan batch_correction_condition_columns "
                "(internal model) must not contain duplicates"
            )
        replicate_column = None
        if batch_correction_replicate_column is not None:
            replicate_column = str(batch_correction_replicate_column).strip()
            if replicate_column == "":
                raise PhosPyInputError(
                    "dataset preprocessing plan batch_correction_replicate_column "
                    "(internal model) must be a non-empty string when provided"
                )
        if method in SPS_RUV_BATCH_CORRECTION_METHODS:
            if batch_correction_internal_request is None:
                raise PhosPyInputError(
                    "dataset preprocessing plan SPS/RUV-style batch correction "
                    "requires batch_correction_internal_request"
                )
            if batch_correction_control_site_set is None:
                raise PhosPyInputError(
                    "dataset preprocessing plan SPS/RUV-style batch correction "
                    "requires batch_correction_control_site_set"
                )
            if batch_correction_missingness_policy is None:
                raise PhosPyInputError(
                    "dataset preprocessing plan SPS/RUV-style batch correction "
                    "requires batch_correction_missingness_policy"
                )
        if (
            method != DATASET_BATCH_CORRECTION_METHOD_NONE
            and batch_correction_preserve_condition_effects is not True
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan "
                "batch_correction_preserve_condition_effects (internal model) "
                "must be True for linear_residualize_batch"
            )
        return ResolvedBatchCorrectionPlanFields(
            batch_correction_method=method,
            batch_correction_batch_column=batch_column,
            batch_correction_condition_column=condition_column,
            batch_correction_condition_columns=condition_columns,
            batch_correction_replicate_column=replicate_column,
        )


class PreprocessingGroupCoveragePlanRuleFamily:
    """Validate group-aware coverage filtering plan coherence."""

    def run(
        self,
        *,
        enabled: bool,
        group_column: str | None,
        min_finite_observations_per_group: int | None,
        min_finite_fraction_per_group: float | None,
        min_groups_passing_threshold: int,
        stage_order: tuple[str, ...],
    ) -> None:
        validate_group_coverage_filter_config(
            enabled=enabled,
            group_column=group_column,
            min_finite_observations_per_group=min_finite_observations_per_group,
            min_finite_fraction_per_group=min_finite_fraction_per_group,
            min_groups_passing_threshold=min_groups_passing_threshold,
        )
        if (
            enabled
            and DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER not in stage_order
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan requests group-aware coverage "
                "filtering but stage_order does not include "
                "'group_coverage_filter'. Build plans from "
                "DatasetPreprocessingConfig or include the group_coverage_filter "
                "stage explicitly."
            )


__all__ = [
    "PreprocessingBatchCorrectionPlanRuleFamily",
    "PreprocessingCorePlanPolicyRuleFamily",
    "PreprocessingDownstreamPlanPolicyRuleFamily",
    "PreprocessingGroupCoveragePlanRuleFamily",
    "PreprocessingLocalisationPlanRuleFamily",
    "ResolvedBatchCorrectionPlanFields",
    "ResolvedCorePreprocessingPlanPolicies",
    "ResolvedDownstreamPreprocessingPlanPolicies",
    "ResolvedLocalisationPlanFields",
]
