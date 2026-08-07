"""Explicit assembly of typed resolved preprocessing sections into a plan."""

from __future__ import annotations

from phospy.science.datasets.preprocessing.plan import PreprocessingPlan
from phospy.science.datasets.preprocessing.plan_resolved import (
    ResolvedPreprocessingPlanFields,
)


class PreprocessingPlanAssembler:
    """Build a validated :class:`PreprocessingPlan` from typed resolved sections."""

    def __init__(
        self,
        *,
        plan_type: type[PreprocessingPlan] = PreprocessingPlan,
    ) -> None:
        self._plan_type = plan_type

    def run(
        self,
        resolved: ResolvedPreprocessingPlanFields,
    ) -> PreprocessingPlan:
        core = resolved.core
        imputation = resolved.imputation
        localisation = resolved.localisation
        site_sequence = resolved.site_sequence
        group_coverage = resolved.group_coverage
        total_protein = resolved.total_protein
        site_matrix_comparisons = resolved.site_matrix_comparisons
        ruv_readiness = resolved.ruv_readiness
        batch_correction = resolved.batch_correction
        stage_order = resolved.stage_order

        return self._plan_type(
            intensity_transform_policy=core.intensity_transform_policy,
            intensity_transform_pseudocount=core.intensity_transform_pseudocount,
            normalisation_policy=core.normalisation_policy,
            missing_data_policy=core.missing_data_policy,
            missing_data_min_observed_values=core.missing_data_min_observed_values,
            missing_data_q=core.missing_data_q,
            missing_data_width=core.missing_data_width,
            missing_data_seed=core.missing_data_seed,
            missing_data_k=core.missing_data_k,
            missing_data_distance=core.missing_data_distance,
            missing_data_max_missing_fraction_per_row=(
                core.missing_data_max_missing_fraction_per_row
            ),
            missing_data_no_overlap_policy=core.missing_data_no_overlap_policy,
            missing_data_input_scale=imputation.missing_data_input_scale,
            missing_data_input_scale_source=(
                imputation.missing_data_input_scale_source
            ),
            missing_data_imputation_operation_order=(
                imputation.missing_data_imputation_operation_order
            ),
            localisation_mode=localisation.localisation_mode,
            localisation_min_confidence=localisation.localisation_min_confidence,
            localisation_confidence_column=(
                localisation.localisation_confidence_column
            ),
            localisation_waiver_reason=localisation.localisation_waiver_reason,
            site_sequence_resolution_enabled=(
                site_sequence.site_sequence_resolution_enabled
            ),
            site_sequence_resolution_fasta_path=(
                site_sequence.site_sequence_resolution_fasta_path
            ),
            site_sequence_resolution_mode=site_sequence.site_sequence_resolution_mode,
            site_sequence_resolution_conflict_policy=(
                site_sequence.site_sequence_resolution_conflict_policy
            ),
            site_sequence_resolution_flank_size=(
                site_sequence.site_sequence_resolution_flank_size
            ),
            site_sequence_resolution_accession_column=(
                site_sequence.site_sequence_resolution_accession_column
            ),
            site_sequence_resolution_site_column=(
                site_sequence.site_sequence_resolution_site_column
            ),
            group_coverage_filter_enabled=(
                group_coverage.group_coverage_filter_enabled
            ),
            group_coverage_filter_group_column=(
                group_coverage.group_coverage_filter_group_column
            ),
            group_coverage_filter_min_finite_observations_per_group=(
                group_coverage.group_coverage_filter_min_finite_observations_per_group
            ),
            group_coverage_filter_min_finite_fraction_per_group=(
                group_coverage.group_coverage_filter_min_finite_fraction_per_group
            ),
            group_coverage_filter_min_groups_passing_threshold=(
                group_coverage.group_coverage_filter_min_groups_passing_threshold
            ),
            total_protein_correction_policy=(
                total_protein.total_protein_correction_policy
            ),
            total_protein_correction_identity_policy=(
                total_protein.total_protein_correction_identity_policy
            ),
            protein_aware_preparation_policy=(
                total_protein.protein_aware_preparation_policy
            ),
            protein_aware_preparation_mapping_policy=(
                total_protein.protein_aware_preparation_mapping_policy
            ),
            site_matrix_policy=site_matrix_comparisons.site_matrix_policy,
            site_matrix_duplicate_site_policy=(
                site_matrix_comparisons.site_matrix_duplicate_site_policy
            ),
            site_matrix_missing_data_policy=(
                site_matrix_comparisons.site_matrix_missing_data_policy
            ),
            site_matrix_minimum_observed_values=(
                site_matrix_comparisons.site_matrix_minimum_observed_values
            ),
            comparison_building_policy=(
                site_matrix_comparisons.comparison_building_policy
            ),
            comparison_sample_group_column=(
                site_matrix_comparisons.comparison_sample_group_column
            ),
            comparison_pairs=site_matrix_comparisons.comparison_pairs,
            ruv_readiness_enabled=ruv_readiness.ruv_readiness_enabled,
            ruv_readiness_control_feature_column=(
                ruv_readiness.ruv_readiness_control_feature_column
            ),
            ruv_readiness_replicate_group_column=(
                ruv_readiness.ruv_readiness_replicate_group_column
            ),
            ruv_readiness_batch_column=ruv_readiness.ruv_readiness_batch_column,
            batch_correction_method=batch_correction.batch_correction_method,
            batch_correction_batch_column=(
                batch_correction.batch_correction_batch_column
            ),
            batch_correction_condition_column=(
                batch_correction.batch_correction_condition_column
            ),
            batch_correction_condition_columns=(
                batch_correction.batch_correction_condition_columns
            ),
            batch_correction_replicate_column=(
                batch_correction.batch_correction_replicate_column
            ),
            batch_correction_control_site_set=(
                batch_correction.batch_correction_control_site_set
            ),
            batch_correction_missingness_policy=(
                batch_correction.batch_correction_missingness_policy
            ),
            batch_correction_internal_request=(
                batch_correction.batch_correction_internal_request
            ),
            batch_correction_preserve_condition_effects=(
                batch_correction.batch_correction_preserve_condition_effects
            ),
            stage_order=stage_order.stage_order,
            stage_order_resolution=stage_order.stage_order_resolution,
        )


__all__ = ["PreprocessingPlanAssembler"]
