"""Execution-ready preprocessing plan model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

from phospy.preparation_records.policies import (
    MissingDataPolicy,
    NormalisationPolicy,
    SiteSequenceConflictPolicy,
    SiteSequenceResolutionMode,
    TotalProteinCorrectionPolicy,
)
from phospy.preprocessing.batch_correction.contracts.configuration import (
    DATASET_BATCH_CORRECTION_METHOD_NONE,
)
from phospy.preprocessing.batch_correction.contracts.missingness import (
    CorrectionMissingnessPolicy,
)
from phospy.preprocessing.batch_correction.contracts.request import (
    InternalBatchCorrectionRequest,
)
from phospy.preprocessing.contracts.configuration import (
    DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED,
    DatasetComparisonPair,
    DatasetProteinAwarePreparationMappingPolicy,
    DatasetProteinAwarePreparationPolicy,
)
from phospy.preprocessing.contracts.policies import (
    ComparisonBuildingPolicy,
    ImputationInputScale,
    IntensityTransformPolicy,
    LocalisationEligibilityMode,
    SiteMatrixDuplicateSitePolicy,
    SiteMatrixMissingDataPolicy,
    SiteMatrixPolicy,
)
from phospy.preprocessing.contracts.stage_keys import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT,
    DATASET_PREPROCESSING_STAGE_PROTEIN_AWARE_PREPARATION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PREPROCESSING_STAGE_ORDER_RATIONALE_BATCH_CORRECTION,
    PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
    PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER,
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM,
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
)
from phospy.preprocessing.planning.ordering import (
    PreprocessingStageOrderResolution,
    PreprocessingStageOrderValidator,
    reject_external_corrected_output_after_downstream_preprocessing,
)
from phospy.preprocessing.planning.validation import (
    validate_resolved_preprocessing_plan,
)
from phospy.preprocessing.protein_aware.identity import (
    DEFAULT_TOTAL_PROTEIN_CORRECTION_IDENTITY_POLICY,
    TotalProteinCorrectionIdentityPolicy,
)

if TYPE_CHECKING:
    from phospy.preprocessing.contracts.configuration.dataset import (
        DatasetPreprocessingConfig,
    )
    from phospy.preprocessing.planning.boundary import (
        PreprocessingPlanInterpreter as PreprocessingPlanInterpreter,
    )


@dataclass(frozen=True, slots=True)
class PreprocessingPlan:
    """Execution-ready internal preprocessing plan."""

    intensity_transform_policy: IntensityTransformPolicy = (
        IntensityTransformPolicy.IDENTITY
    )
    intensity_transform_pseudocount: float = 1.0
    normalisation_policy: NormalisationPolicy = NormalisationPolicy.NONE
    missing_data_policy: MissingDataPolicy = MissingDataPolicy.FORBID
    missing_data_min_observed_values: int | None = None
    missing_data_q: float | None = None
    missing_data_width: float | None = None
    missing_data_seed: int | None = None
    missing_data_k: int | None = None
    missing_data_distance: str | None = None
    missing_data_max_missing_fraction_per_row: float | None = None
    missing_data_no_overlap_policy: str | None = None
    missing_data_group_column: str | None = None
    missing_data_min_partial_observed_fraction: float | None = None
    missing_data_min_reference_observed_fraction: float | None = None
    missing_data_input_scale: ImputationInputScale | None = None
    missing_data_input_scale_source: str | None = None
    missing_data_imputation_operation_order: str | None = None
    localisation_mode: LocalisationEligibilityMode = (
        LocalisationEligibilityMode.REQUIRE_THRESHOLD
    )
    localisation_min_confidence: float = 0.75
    localisation_confidence_column: str = "localisation_confidence"
    localisation_waiver_reason: str | None = None
    site_sequence_resolution_enabled: bool = False
    site_sequence_resolution_fasta_path: str | None = None
    site_sequence_resolution_mode: SiteSequenceResolutionMode = (
        SiteSequenceResolutionMode.VALIDATE_EXISTING_AND_FILL_MISSING
    )
    site_sequence_resolution_conflict_policy: SiteSequenceConflictPolicy = (
        SiteSequenceConflictPolicy.PRESERVE_EXISTING
    )
    site_sequence_resolution_flank_size: int = 7
    site_sequence_resolution_accession_column: str = "protein_accession"
    site_sequence_resolution_site_column: str = "site"
    group_coverage_filter_enabled: bool = False
    group_coverage_filter_group_column: str | None = None
    group_coverage_filter_min_finite_observations_per_group: int | None = None
    group_coverage_filter_min_finite_fraction_per_group: float | None = None
    group_coverage_filter_min_groups_passing_threshold: int = 1
    total_protein_correction_policy: TotalProteinCorrectionPolicy = (
        TotalProteinCorrectionPolicy.NONE
    )
    total_protein_correction_identity_policy: TotalProteinCorrectionIdentityPolicy = (
        DEFAULT_TOTAL_PROTEIN_CORRECTION_IDENTITY_POLICY
    )
    protein_aware_preparation_policy: DatasetProteinAwarePreparationPolicy = (
        DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED
    )
    protein_aware_preparation_mapping_policy: DatasetProteinAwarePreparationMappingPolicy = DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS
    site_matrix_policy: SiteMatrixPolicy = SiteMatrixPolicy.AS_INPUT
    comparison_building_policy: ComparisonBuildingPolicy = ComparisonBuildingPolicy.NONE
    site_matrix_duplicate_site_policy: SiteMatrixDuplicateSitePolicy = (
        SiteMatrixDuplicateSitePolicy.ERROR
    )
    site_matrix_missing_data_policy: SiteMatrixMissingDataPolicy = (
        SiteMatrixMissingDataPolicy.DROP_ANY_MISSING
    )
    site_matrix_minimum_observed_values: int | None = None
    comparison_sample_group_column: str = (
        DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN
    )
    comparison_pairs: tuple[DatasetComparisonPair, ...] | None = None
    ruv_readiness_enabled: bool = False
    ruv_readiness_control_feature_column: str = "is_control_feature"
    ruv_readiness_replicate_group_column: str = "replicate_group"
    ruv_readiness_batch_column: str | None = "batch"
    batch_correction_method: str = DATASET_BATCH_CORRECTION_METHOD_NONE
    batch_correction_batch_column: str = "batch"
    batch_correction_condition_column: str = "condition"
    batch_correction_condition_columns: tuple[str, ...] = ("condition",)
    batch_correction_replicate_column: str | None = None
    batch_correction_control_site_set: object | None = None
    batch_correction_missingness_policy: CorrectionMissingnessPolicy | None = None
    batch_correction_internal_request: InternalBatchCorrectionRequest | None = None
    batch_correction_preserve_condition_effects: bool = True
    stage_order: tuple[str, ...] = DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT
    stage_order_resolution: tuple[PreprocessingStageOrderResolution, ...] = ()

    def __post_init__(self) -> None:
        validate_resolved_preprocessing_plan(self)

    @property
    def resolved_batch_correction_method(self) -> str:
        """Return the core plan's normalised correction method."""

        method = str(self.batch_correction_method).strip()
        return method or DATASET_BATCH_CORRECTION_METHOD_NONE

    @property
    def batch_correction_enabled(self) -> bool:
        """Return whether the resolved core plan requests correction."""

        return (
            self.resolved_batch_correction_method
            != DATASET_BATCH_CORRECTION_METHOD_NONE
        )

    @classmethod
    def from_config(cls, config: DatasetPreprocessingConfig) -> PreprocessingPlan:
        from phospy.preprocessing.planning.boundary import (
            PreprocessingPlanInterpreter,
        )

        return PreprocessingPlanInterpreter(plan_type=cls).run(
            cast("DatasetPreprocessingConfig", config)
        )

    @classmethod
    def default(cls) -> PreprocessingPlan:
        from phospy.preprocessing.contracts.configuration.dataset import (
            DatasetPreprocessingConfig,
        )

        return cls.from_config(DatasetPreprocessingConfig())


def __getattr__(name: str) -> object:
    if name == "PreprocessingPlanInterpreter":
        from phospy.preprocessing.planning.boundary import (
            PreprocessingPlanInterpreter,
        )

        return PreprocessingPlanInterpreter
    raise AttributeError(name)


__all__ = [
    "DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION",
    "DATASET_PREPROCESSING_STAGE_COMPARISONS",
    "DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER",
    "DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM",
    "DATASET_PREPROCESSING_STAGE_LOCALISATION",
    "DATASET_PREPROCESSING_STAGE_MISSING_DATA",
    "DATASET_PREPROCESSING_STAGE_NORMALISATION",
    "DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT",
    "DATASET_PREPROCESSING_STAGE_PROTEIN_AWARE_PREPARATION",
    "DATASET_PREPROCESSING_STAGE_SITE_MATRIX",
    "DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION",
    "DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_BATCH_CORRECTION",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM",
    "PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA",
    "PreprocessingPlan",
    "PreprocessingPlanInterpreter",
    "PreprocessingStageOrderResolution",
    "PreprocessingStageOrderValidator",
    "TotalProteinCorrectionIdentityPolicy",
    "reject_external_corrected_output_after_downstream_preprocessing",
]
