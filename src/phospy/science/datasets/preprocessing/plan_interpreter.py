"""Preprocessing plan interpretation coordinator."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError
from phospy.science.configs.dataset import DatasetPreprocessingConfig
from phospy.science.configs.preprocessing import (
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
)
from phospy.science.configs.preprocessing._validation import (
    validate_protein_aware_preparation_config,
)
from phospy.science.datasets.preprocessing.imputation_scale_policy import (
    reject_unestablished_log2_imputation_input_scale,
    resolve_imputation_scale_policy,
)
from phospy.science.datasets.preprocessing.plan import PreprocessingPlan
from phospy.science.datasets.preprocessing.plan_assembly import (
    PreprocessingPlanAssembler,
)
from phospy.science.datasets.preprocessing.plan_config_resolution import (
    PreprocessingConfigPolicyResolver,
    resolve_site_sequence_resolution_conflict_policy,
)
from phospy.science.datasets.preprocessing.plan_resolved import (
    ResolvedCoreTransformPlanFields,
    ResolvedImputationScalePlanFields,
    ResolvedPreprocessingPlanFields,
    ResolvedStageOrderPlanFields,
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
    PreprocessingStageOrderPlanner,
)
from phospy.science.datasets.preprocessing.total_protein_identity import (
    TotalProteinCorrectionIdentityResolver,
)
from phospy.science.transformations.models import IntensityScaleKind


class PreprocessingPlanInterpreter:
    """Convert public preprocessing config into an execution-ready plan."""

    def __init__(
        self,
        *,
        plan_type: type[PreprocessingPlan] = PreprocessingPlan,
    ) -> None:
        self._plan_type = plan_type

    def run(
        self,
        config: DatasetPreprocessingConfig,
        *,
        declared_input_scale_kind: IntensityScaleKind | None = None,
    ) -> PreprocessingPlan:
        if not isinstance(config, DatasetPreprocessingConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config must be a "
                "DatasetPreprocessingConfig"
            )

        validate_protein_aware_preparation_config(
            policy=config.protein_aware_preparation.policy,
            protein_mapping_policy=(
                config.protein_aware_preparation.protein_mapping_policy
            ),
            supported_policies=DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
            supported_mapping_policies=(
                DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES
            ),
        )

        policies = PreprocessingConfigPolicyResolver().run(config)
        core_policies = PreprocessingCorePlanPolicyRuleFamily().run(
            intensity_transform_policy=policies.intensity_transform_policy,
            normalisation_policy=policies.normalisation_policy,
            missing_data_policy=policies.missing_data_policy,
        )
        imputation_scale = resolve_imputation_scale_policy(
            missing_data_policy=core_policies.missing_data_policy,
            requested_input_scale=config.missing_data.input_scale,
            intensity_transform_policy=core_policies.intensity_transform_policy,
        )
        reject_unestablished_log2_imputation_input_scale(
            resolved_policy=imputation_scale,
            intensity_transform_policy=core_policies.intensity_transform_policy,
            declared_input_scale_kind=declared_input_scale_kind,
        )
        localisation = PreprocessingLocalisationPlanRuleFamily().run(
            localisation_mode=policies.localisation_mode,
            localisation_min_confidence=float(config.localisation.min_confidence),
            localisation_confidence_column=str(
                config.localisation.confidence_column
            ).strip(),
            localisation_waiver_reason=(
                None
                if config.localisation.waiver_reason is None
                else str(config.localisation.waiver_reason).strip()
            ),
        )
        site_sequence = PreprocessingSiteSequencePlanRuleFamily().run(
            site_sequence_resolution_enabled=(
                policies.site_sequence_resolution_enabled
            ),
            site_sequence_resolution_fasta_path=(
                config.site_sequence_resolution.fasta_path
            ),
            site_sequence_resolution_mode=policies.site_sequence_resolution_mode,
            site_sequence_resolution_conflict_policy=(
                resolve_site_sequence_resolution_conflict_policy(
                    mode=policies.site_sequence_resolution_mode,
                    conflict_policy=config.site_sequence_resolution.conflict_policy,
                )
            ),
            site_sequence_resolution_flank_size=int(
                config.site_sequence_resolution.flank_size
            ),
            site_sequence_resolution_accession_column=(
                config.site_sequence_resolution.accession_column
            ),
            site_sequence_resolution_site_column=(
                config.site_sequence_resolution.site_column
            ),
        )
        group_coverage = PreprocessingGroupCoveragePlanRuleFamily().run(
            enabled=bool(config.group_coverage_filter.enabled),
            group_column=config.group_coverage_filter.group_column,
            min_finite_observations_per_group=(
                config.group_coverage_filter.min_finite_observations_per_group
            ),
            min_finite_fraction_per_group=(
                None
                if config.group_coverage_filter.min_finite_fraction_per_group is None
                else float(config.group_coverage_filter.min_finite_fraction_per_group)
            ),
            min_groups_passing_threshold=(
                config.group_coverage_filter.min_groups_passing_threshold
            ),
        )
        total_protein = PreprocessingTotalProteinCorrectionPlanRuleFamily().run(
            total_protein_correction_policy=policies.total_correction_policy,
            total_protein_correction_identity_policy=(
                TotalProteinCorrectionIdentityResolver().run(
                    config.total_protein_correction.identity
                )
            ),
            protein_aware_preparation_policy=config.protein_aware_preparation.policy,
            protein_aware_preparation_mapping_policy=(
                config.protein_aware_preparation.protein_mapping_policy
            ),
        )
        site_matrix_comparisons = PreprocessingSiteMatrixComparisonPlanRuleFamily().run(
            site_matrix_policy=policies.site_matrix_policy,
            site_matrix_duplicate_site_policy=(
                policies.site_matrix_duplicate_site_policy
            ),
            site_matrix_missing_data_policy=(policies.site_matrix_missing_data_policy),
            site_matrix_minimum_observed_values=(
                config.site_matrix.minimum_observed_values
            ),
            comparison_building_policy=policies.comparison_building_policy,
            comparison_sample_group_column=(config.comparisons.sample_group_column),
            comparison_pairs=(
                None
                if config.comparisons.pairs is None
                else tuple(config.comparisons.pairs)
            ),
        )
        ruv_readiness = PreprocessingRuvReadinessPlanRuleFamily().run(
            ruv_readiness_enabled=bool(config.ruv_readiness.enabled),
            ruv_readiness_control_feature_column=(
                config.ruv_readiness.control_feature_column
            ),
            ruv_readiness_replicate_group_column=(
                config.ruv_readiness.replicate_group_column
            ),
            ruv_readiness_batch_column=config.ruv_readiness.batch_column,
        )
        batch_correction = PreprocessingBatchCorrectionPlanRuleFamily().run(
            batch_correction_method=policies.batch_correction_method,
            batch_correction_batch_column=config.batch_correction.batch_column,
            batch_correction_condition_column=(
                policies.batch_correction_condition_column
            ),
            batch_correction_condition_columns=(
                policies.batch_correction_condition_columns
            ),
            batch_correction_replicate_column=(
                policies.batch_correction_replicate_column
            ),
            batch_correction_control_site_set=(
                policies.batch_correction_control_site_set
            ),
            batch_correction_missingness_policy=(
                policies.batch_correction_missingness_policy
            ),
            batch_correction_internal_request=(
                policies.batch_correction_internal_request
            ),
            batch_correction_preserve_condition_effects=(
                policies.batch_correction_preserve_condition_effects
            ),
        )
        stage_plan = PreprocessingStageOrderPlanner().run(
            site_sequence_resolution_enabled=(
                site_sequence.site_sequence_resolution_enabled
            ),
            intensity_transform_policy=core_policies.intensity_transform_policy,
            normalisation_policy=core_policies.normalisation_policy,
            site_matrix_policy=site_matrix_comparisons.site_matrix_policy,
            comparison_building_policy=(
                site_matrix_comparisons.comparison_building_policy
            ),
            localisation_mode=localisation.localisation_mode,
            missing_data_policy=core_policies.missing_data_policy,
            missing_data_input_scale=imputation_scale.input_scale,
            batch_correction_method=batch_correction.batch_correction_method,
            total_correction_policy=(total_protein.total_protein_correction_policy),
            group_coverage_filter_enabled=(
                group_coverage.group_coverage_filter_enabled
            ),
        )
        resolved = ResolvedPreprocessingPlanFields(
            core=ResolvedCoreTransformPlanFields(
                intensity_transform_policy=core_policies.intensity_transform_policy,
                intensity_transform_pseudocount=float(
                    config.intensity_transform.pseudocount
                ),
                normalisation_policy=core_policies.normalisation_policy,
                missing_data_policy=core_policies.missing_data_policy,
                missing_data_min_observed_values=(
                    config.missing_data.min_observed_values
                ),
                missing_data_q=(
                    None
                    if config.missing_data.q is None
                    else float(config.missing_data.q)
                ),
                missing_data_width=(
                    None
                    if config.missing_data.width is None
                    else float(config.missing_data.width)
                ),
                missing_data_seed=config.missing_data.seed,
                missing_data_k=config.missing_data.k,
                missing_data_distance=config.missing_data.distance,
                missing_data_max_missing_fraction_per_row=(
                    None
                    if config.missing_data.max_missing_fraction_per_row is None
                    else float(config.missing_data.max_missing_fraction_per_row)
                ),
            ),
            imputation=ResolvedImputationScalePlanFields(
                missing_data_input_scale=imputation_scale.input_scale,
                missing_data_input_scale_source=imputation_scale.input_scale_source,
                missing_data_imputation_operation_order=(
                    imputation_scale.operation_order
                ),
            ),
            localisation=localisation,
            site_sequence=site_sequence,
            group_coverage=group_coverage,
            total_protein=total_protein,
            site_matrix_comparisons=site_matrix_comparisons,
            ruv_readiness=ruv_readiness,
            batch_correction=batch_correction,
            stage_order=ResolvedStageOrderPlanFields(
                stage_order=stage_plan.stage_order,
                stage_order_resolution=stage_plan.stage_order_resolution,
            ),
        )
        return PreprocessingPlanAssembler(plan_type=self._plan_type).run(resolved)


__all__ = ["PreprocessingPlanInterpreter"]
