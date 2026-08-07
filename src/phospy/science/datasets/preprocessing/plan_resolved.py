"""Typed resolved preprocessing-plan sections."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.science.configs.preprocessing import (
    CorrectionMissingnessPolicy,
    DatasetComparisonPair,
    DatasetProteinAwarePreparationMappingPolicy,
    DatasetProteinAwarePreparationPolicy,
    InternalBatchCorrectionRequest,
)
from phospy.science.datasets.preprocessing.plan_stage_order import (
    PreprocessingStageOrderResolution,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    ImputationInputScale,
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
from phospy.science.datasets.preprocessing.total_protein_identity import (
    TotalProteinCorrectionIdentityPolicy,
)


@dataclass(frozen=True, slots=True)
class ResolvedCoreTransformPlanFields:
    """Core quantitative transform, normalisation, and missing-data fields."""

    intensity_transform_policy: IntensityTransformPolicy
    intensity_transform_pseudocount: float
    normalisation_policy: NormalisationPolicy
    missing_data_policy: MissingDataPolicy
    missing_data_min_observed_values: int | None
    missing_data_q: float | None
    missing_data_width: float | None
    missing_data_seed: int | None
    missing_data_k: int | None
    missing_data_distance: str | None
    missing_data_max_missing_fraction_per_row: float | None
    missing_data_no_overlap_policy: str | None = None


@dataclass(frozen=True, slots=True)
class ResolvedImputationScalePlanFields:
    """Resolved imputation input scale and operation-order fields."""

    missing_data_input_scale: ImputationInputScale | None
    missing_data_input_scale_source: str | None
    missing_data_imputation_operation_order: str | None


@dataclass(frozen=True, slots=True)
class ResolvedLocalisationPlanFields:
    """Localisation evidence policy resolved for internal execution."""

    localisation_mode: LocalisationEligibilityMode
    localisation_min_confidence: float
    localisation_confidence_column: str
    localisation_waiver_reason: str | None


@dataclass(frozen=True, slots=True)
class ResolvedSiteSequencePlanFields:
    """Resolved site-sequence configuration consumed by preprocessing."""

    site_sequence_resolution_enabled: bool
    site_sequence_resolution_fasta_path: str | None
    site_sequence_resolution_mode: SiteSequenceResolutionMode
    site_sequence_resolution_conflict_policy: SiteSequenceConflictPolicy
    site_sequence_resolution_flank_size: int
    site_sequence_resolution_accession_column: str
    site_sequence_resolution_site_column: str


@dataclass(frozen=True, slots=True)
class ResolvedGroupCoveragePlanFields:
    """Resolved group-aware coverage filtering fields."""

    group_coverage_filter_enabled: bool
    group_coverage_filter_group_column: str | None
    group_coverage_filter_min_finite_observations_per_group: int | None
    group_coverage_filter_min_finite_fraction_per_group: float | None
    group_coverage_filter_min_groups_passing_threshold: int


@dataclass(frozen=True, slots=True)
class ResolvedTotalProteinCorrectionPlanFields:
    """Resolved total-protein correction and protein-aware preparation fields."""

    total_protein_correction_policy: TotalProteinCorrectionPolicy
    total_protein_correction_identity_policy: TotalProteinCorrectionIdentityPolicy
    protein_aware_preparation_policy: DatasetProteinAwarePreparationPolicy
    protein_aware_preparation_mapping_policy: (
        DatasetProteinAwarePreparationMappingPolicy
    )


@dataclass(frozen=True, slots=True)
class ResolvedSiteMatrixComparisonPlanFields:
    """Resolved site-matrix and comparison-building fields."""

    site_matrix_policy: SiteMatrixPolicy
    site_matrix_duplicate_site_policy: SiteMatrixDuplicateSitePolicy
    site_matrix_missing_data_policy: SiteMatrixMissingDataPolicy
    site_matrix_minimum_observed_values: int | None
    comparison_building_policy: ComparisonBuildingPolicy
    comparison_sample_group_column: str
    comparison_pairs: tuple[DatasetComparisonPair, ...] | None


@dataclass(frozen=True, slots=True)
class ResolvedRuvReadinessPlanFields:
    """Resolved report-only RUV-readiness fields."""

    ruv_readiness_enabled: bool
    ruv_readiness_control_feature_column: str
    ruv_readiness_replicate_group_column: str
    ruv_readiness_batch_column: str | None


@dataclass(frozen=True, slots=True)
class ResolvedBatchCorrectionPlanFields:
    """Batch-correction plan fields resolved for internal execution."""

    batch_correction_method: str
    batch_correction_batch_column: str
    batch_correction_condition_column: str
    batch_correction_condition_columns: tuple[str, ...]
    batch_correction_replicate_column: str | None
    batch_correction_control_site_set: object | None
    batch_correction_missingness_policy: CorrectionMissingnessPolicy | None
    batch_correction_internal_request: InternalBatchCorrectionRequest | None
    batch_correction_preserve_condition_effects: bool


@dataclass(frozen=True, slots=True)
class ResolvedStageOrderPlanFields:
    """Resolved stage order and rationale records."""

    stage_order: tuple[str, ...]
    stage_order_resolution: tuple[PreprocessingStageOrderResolution, ...]


@dataclass(frozen=True, slots=True)
class ResolvedPreprocessingPlanFields:
    """Typed aggregate connecting rule resolution to plan construction."""

    core: ResolvedCoreTransformPlanFields
    imputation: ResolvedImputationScalePlanFields
    localisation: ResolvedLocalisationPlanFields
    site_sequence: ResolvedSiteSequencePlanFields
    group_coverage: ResolvedGroupCoveragePlanFields
    total_protein: ResolvedTotalProteinCorrectionPlanFields
    site_matrix_comparisons: ResolvedSiteMatrixComparisonPlanFields
    ruv_readiness: ResolvedRuvReadinessPlanFields
    batch_correction: ResolvedBatchCorrectionPlanFields
    stage_order: ResolvedStageOrderPlanFields


__all__ = [
    "ResolvedBatchCorrectionPlanFields",
    "ResolvedCoreTransformPlanFields",
    "ResolvedGroupCoveragePlanFields",
    "ResolvedImputationScalePlanFields",
    "ResolvedLocalisationPlanFields",
    "ResolvedPreprocessingPlanFields",
    "ResolvedRuvReadinessPlanFields",
    "ResolvedSiteMatrixComparisonPlanFields",
    "ResolvedSiteSequencePlanFields",
    "ResolvedStageOrderPlanFields",
    "ResolvedTotalProteinCorrectionPlanFields",
]
