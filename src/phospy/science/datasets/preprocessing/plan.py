"""Preprocessing stage-order planning and config interpretation."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table_tolerance
from phospy.science.configs.dataset import DatasetPreprocessingConfig
from phospy.science.configs.preprocessing import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICY_REQUIRE_UNAMBIGUOUS,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICY_DISABLED,
    DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICY_ERROR,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE,
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR,
    SPS_RUV_BATCH_CORRECTION_METHODS,
    CorrectionMissingnessPolicy,
    DatasetBatchCorrectionConfig,
    DatasetComparisonPair,
    DatasetProteinAwarePreparationMappingPolicy,
    DatasetProteinAwarePreparationPolicy,
    DatasetSiteSequenceConflictPolicy,
    DatasetTotalProteinCorrectionDuplicatePolicy,
    DatasetTotalProteinCorrectionIdentityConfig,
    DatasetTotalProteinCorrectionIdentityMode,
    DatasetTotalProteinCorrectionUnmatchedPolicy,
    InternalBatchCorrectionRequest,
    SpsRuvBatchCorrectionConfig,
)
from phospy.science.configs.preprocessing._validation import (
    reject_unsupported_ruv_iii_style_method,
    validate_group_coverage_filter_config,
    validate_protein_aware_preparation_config,
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
    TotalProteinCorrectionIdentityMatchingPolicy,
    TotalProteinCorrectionPolicy,
)

_UNSUPPORTED_BATCH_CORRECTION_METHOD_RUV_III_STYLE = "ruv_iii_style"

DATASET_PREPROCESSING_STAGE_MISSING_DATA = "missing_data"

DATASET_PREPROCESSING_STAGE_LOCALISATION = "localisation_confidence"

DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION = "site_sequence_resolution"

DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION = "total_protein_correction"

DATASET_PREPROCESSING_STAGE_PROTEIN_AWARE_PREPARATION = "protein_aware_preparation"

DATASET_PREPROCESSING_STAGE_SITE_MATRIX = "site_matrix"

DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM = "intensity_transform"

DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER = "group_coverage_filter"

DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION = "batch_correction"

DATASET_PREPROCESSING_STAGE_NORMALISATION = "normalisation"

DATASET_PREPROCESSING_STAGE_COMPARISONS = "comparisons"

DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT = (
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
)

PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM = (
    "impute_minprob requires log2 intensity space, so intensity_transform runs "
    "before missing_data."
)

PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA = (
    "impute_minprob is applied after log2 transformation because its left-censored "
    "sampling model operates on log2 intensities."
)

PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA = (
    "non-minprob missing-data policies run before optional log2 transformation."
)

PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM = (
    "optional log2 transformation runs after non-minprob missing-data handling."
)

PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER = (
    "group-aware coverage filtering runs before missing-data handling so it uses "
    "observed input values."
)

PREPROCESSING_STAGE_ORDER_RATIONALE_BATCH_CORRECTION = (
    "batch correction runs after any configured intensity transformation and "
    "before total-protein correction, site-matrix construction, normalisation, "
    "and comparison building."
)

PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE = (
    "stage included because preprocessing configuration enables it."
)

_BATCH_CORRECTION_DOWNSTREAM_BOUNDARY_STAGES = (
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
)

_EXTERNAL_CORRECTION_DOWNSTREAM_MATRIX_CONSUMING_STAGES = frozenset(
    {
        *_BATCH_CORRECTION_DOWNSTREAM_BOUNDARY_STAGES,
        "normalization",
        "differential_analysis_preparation",
        "kinase_activity_preparation",
        "enrichment_preparation",
        "signalome_preparation",
    }
)


class PreprocessingStageOrderValidator:
    """Validate scientific preprocessing stage-order constraints."""

    def run(
        self,
        *,
        stage_order: tuple[str, ...],
        batch_correction_requested: bool,
    ) -> None:
        stages = tuple(str(stage).strip() for stage in stage_order)
        blank_positions = [
            position for position, stage in enumerate(stages) if stage == ""
        ]
        if blank_positions:
            raise PhosPyInputError(
                "dataset preprocessing plan stage_order contains blank stage "
                f"entries at positions {blank_positions}"
            )
        duplicates = [
            stage for stage in dict.fromkeys(stages) if stages.count(stage) > 1
        ]
        if duplicates:
            raise PhosPyInputError(
                "dataset preprocessing plan stage_order contains duplicate stages: "
                + ", ".join(repr(stage) for stage in duplicates)
            )
        if not batch_correction_requested:
            return
        if DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION not in stages:
            raise PhosPyInputError(
                "dataset preprocessing plan requests batch correction but "
                "stage_order does not include 'batch_correction'. Build plans "
                "from DatasetPreprocessingConfig or include the batch_correction "
                "stage explicitly."
            )

        batch_position = stages.index(DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION)
        if (
            DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM in stages
            and batch_position
            < stages.index(DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM)
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan has unsupported stage_order: "
                "batch_correction must run after intensity_transform when both "
                "stages are configured"
            )

        downstream_before_batch = [
            stage
            for stage in _BATCH_CORRECTION_DOWNSTREAM_BOUNDARY_STAGES
            if stage in stages and stages.index(stage) < batch_position
        ]
        if downstream_before_batch:
            raise PhosPyInputError(
                "dataset preprocessing plan has unsupported stage_order: "
                "batch_correction cannot run after downstream stages have consumed "
                "the matrix because that would weaken the analysis-ready dataset "
                "boundary; downstream stages before batch_correction: "
                + ", ".join(downstream_before_batch)
            )


def reject_external_corrected_output_after_downstream_preprocessing(
    stage_order: tuple[str, ...],
) -> None:
    """Reject external correction when the plan has downstream matrix consumers."""

    configured = tuple(
        stage
        for stage in (str(item).strip() for item in stage_order)
        if stage in _EXTERNAL_CORRECTION_DOWNSTREAM_MATRIX_CONSUMING_STAGES
    )
    if not configured:
        return
    raise PhosPyInputError(
        "external corrected output cannot be integrated after downstream "
        "preprocessing stages. Configured downstream matrix-consuming "
        "preprocessing stages: "
        + ", ".join(configured)
        + ". Provide the corrected output as the only matrix-changing "
        "preprocessing input, or use native SpsRuvBatchCorrectionConfig inside "
        "the preprocessing pipeline."
    )


@dataclass(frozen=True, slots=True)
class TotalProteinCorrectionIdentityPolicy:
    """Resolved identity policy consumed by total/protein correction stages."""

    mode: DatasetTotalProteinCorrectionIdentityMode
    phosphosite_key: str
    total_protein_key: str
    matching_policy: TotalProteinCorrectionIdentityMatchingPolicy
    duplicate_policy: DatasetTotalProteinCorrectionDuplicatePolicy
    unmatched_policy: DatasetTotalProteinCorrectionUnmatchedPolicy
    mapping_table: tuple[tuple[str, str], ...] | None = None
    mapping_phosphosite_key: str | None = None
    mapping_total_protein_key: str | None = None
    mapping_table_fingerprint: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "matching_policy",
            TotalProteinCorrectionIdentityMatchingPolicy.parse(
                self.matching_policy,
                field_name=(
                    "dataset preprocessing plan total_protein_correction "
                    "identity matching_policy (internal model)"
                ),
            ),
        )


@dataclass(frozen=True, slots=True)
class PreprocessingStageOrderResolution:
    """Structured explanation of resolved preprocessing stage order."""

    stage: str
    order_index: int
    rationale: str


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
        TotalProteinCorrectionIdentityPolicy(
            mode=DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT,
            phosphosite_key="gene_symbol",
            total_protein_key="__index__",
            matching_policy=TotalProteinCorrectionIdentityMatchingPolicy.STRICT,
            duplicate_policy=DATASET_TOTAL_PROTEIN_CORRECTION_DUPLICATE_POLICY_ERROR,
            unmatched_policy=DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ERROR,
            mapping_table=None,
            mapping_phosphosite_key=None,
            mapping_total_protein_key=None,
            mapping_table_fingerprint=None,
        )
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
        object.__setattr__(
            self,
            "intensity_transform_policy",
            IntensityTransformPolicy.parse(
                self.intensity_transform_policy,
                field_name=(
                    "dataset preprocessing plan intensity_transform_policy "
                    "(internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "normalisation_policy",
            NormalisationPolicy.parse(
                self.normalisation_policy,
                field_name=(
                    "dataset preprocessing plan normalisation_policy (internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "missing_data_policy",
            MissingDataPolicy.parse(
                self.missing_data_policy,
                field_name=(
                    "dataset preprocessing plan missing_data_policy (internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "localisation_mode",
            LocalisationEligibilityMode.parse(
                self.localisation_mode,
                field_name=(
                    "dataset preprocessing plan localisation_mode (internal model)"
                ),
            ),
        )
        min_confidence = float(self.localisation_min_confidence)
        if not math.isfinite(min_confidence) or not (0.0 <= min_confidence <= 1.0):
            raise PhosPyInputError(
                "dataset preprocessing plan localisation_min_confidence "
                "(internal model) must be between 0.0 and 1.0"
            )
        object.__setattr__(self, "localisation_min_confidence", min_confidence)
        confidence_column = str(self.localisation_confidence_column).strip()
        if confidence_column == "":
            raise PhosPyInputError(
                "dataset preprocessing plan localisation_confidence_column "
                "(internal model) must be a non-empty string"
            )
        object.__setattr__(self, "localisation_confidence_column", confidence_column)
        if (
            self.localisation_mode
            is not LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER
            and self.localisation_waiver_reason is not None
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan localisation_waiver_reason "
                "(internal model) must be None unless localisation_mode="
                "'allow_missing_with_waiver'"
            )
        if (
            self.localisation_mode
            is LocalisationEligibilityMode.ALLOW_MISSING_WITH_WAIVER
        ):
            waiver_reason = (
                ""
                if self.localisation_waiver_reason is None
                else self.localisation_waiver_reason
            ).strip()
            if waiver_reason == "":
                raise PhosPyInputError(
                    "dataset preprocessing plan localisation_waiver_reason "
                    "(internal model) must be provided when localisation_mode="
                    "'allow_missing_with_waiver'"
                )
            object.__setattr__(self, "localisation_waiver_reason", waiver_reason)
        object.__setattr__(
            self,
            "site_sequence_resolution_mode",
            SiteSequenceResolutionMode.parse(
                self.site_sequence_resolution_mode,
                field_name=(
                    "dataset preprocessing plan site_sequence_resolution_mode "
                    "(internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "site_sequence_resolution_conflict_policy",
            SiteSequenceConflictPolicy.parse(
                self.site_sequence_resolution_conflict_policy,
                field_name=(
                    "dataset preprocessing plan "
                    "site_sequence_resolution_conflict_policy (internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "site_matrix_policy",
            SiteMatrixPolicy.parse(
                self.site_matrix_policy,
                field_name="dataset preprocessing plan site_matrix_policy (internal model)",
            ),
        )
        object.__setattr__(
            self,
            "comparison_building_policy",
            ComparisonBuildingPolicy.parse(
                self.comparison_building_policy,
                field_name=(
                    "dataset preprocessing plan comparison_building_policy "
                    "(internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "site_matrix_duplicate_site_policy",
            SiteMatrixDuplicateSitePolicy.parse(
                self.site_matrix_duplicate_site_policy,
                field_name=(
                    "dataset preprocessing plan site_matrix_duplicate_site_policy "
                    "(internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "site_matrix_missing_data_policy",
            SiteMatrixMissingDataPolicy.parse(
                self.site_matrix_missing_data_policy,
                field_name=(
                    "dataset preprocessing plan site_matrix_missing_data_policy "
                    "(internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "total_protein_correction_policy",
            TotalProteinCorrectionPolicy.parse(
                self.total_protein_correction_policy,
                field_name=(
                    "dataset preprocessing plan total_protein_correction_policy "
                    "(internal model)"
                ),
            ),
        )
        object.__setattr__(
            self,
            "protein_aware_preparation_policy",
            cast(
                DatasetProteinAwarePreparationPolicy,
                str(self.protein_aware_preparation_policy).strip(),
            ),
        )
        object.__setattr__(
            self,
            "protein_aware_preparation_mapping_policy",
            cast(
                DatasetProteinAwarePreparationMappingPolicy,
                str(self.protein_aware_preparation_mapping_policy).strip(),
            ),
        )
        batch_correction_method = str(self.batch_correction_method).strip()
        if not batch_correction_method:
            batch_correction_method = DATASET_BATCH_CORRECTION_METHOD_NONE
        if (
            batch_correction_method
            == _UNSUPPORTED_BATCH_CORRECTION_METHOD_RUV_III_STYLE
        ):
            reject_unsupported_ruv_iii_style_method(
                batch_correction_method,
                field_name=(
                    "dataset preprocessing plan batch_correction_method "
                    "(internal model)"
                ),
            )
        if (
            batch_correction_method
            not in {
                DATASET_BATCH_CORRECTION_METHOD_NONE,
                DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
            }
            and batch_correction_method not in SPS_RUV_BATCH_CORRECTION_METHODS
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan batch_correction_method "
                "(internal model) must be one of: none, "
                "linear_residualize_batch, sps_ruv_style"
            )
        object.__setattr__(
            self,
            "batch_correction_method",
            batch_correction_method,
        )
        batch_column = str(self.batch_correction_batch_column).strip()
        condition_column = str(self.batch_correction_condition_column).strip()
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
        object.__setattr__(self, "batch_correction_batch_column", batch_column)
        object.__setattr__(self, "batch_correction_condition_column", condition_column)
        condition_columns = tuple(
            str(column).strip() for column in self.batch_correction_condition_columns
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
        object.__setattr__(
            self,
            "batch_correction_condition_columns",
            condition_columns,
        )
        if self.batch_correction_replicate_column is not None:
            replicate_column = str(self.batch_correction_replicate_column).strip()
            if replicate_column == "":
                raise PhosPyInputError(
                    "dataset preprocessing plan batch_correction_replicate_column "
                    "(internal model) must be a non-empty string when provided"
                )
            object.__setattr__(
                self,
                "batch_correction_replicate_column",
                replicate_column,
            )
        if batch_correction_method in SPS_RUV_BATCH_CORRECTION_METHODS:
            if self.batch_correction_internal_request is None:
                raise PhosPyInputError(
                    "dataset preprocessing plan SPS/RUV-style batch correction "
                    "requires batch_correction_internal_request"
                )
            if self.batch_correction_control_site_set is None:
                raise PhosPyInputError(
                    "dataset preprocessing plan SPS/RUV-style batch correction "
                    "requires batch_correction_control_site_set"
                )
            if self.batch_correction_missingness_policy is None:
                raise PhosPyInputError(
                    "dataset preprocessing plan SPS/RUV-style batch correction "
                    "requires batch_correction_missingness_policy"
                )
        if (
            batch_correction_method != DATASET_BATCH_CORRECTION_METHOD_NONE
            and self.batch_correction_preserve_condition_effects is not True
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan "
                "batch_correction_preserve_condition_effects (internal model) "
                "must be True for linear_residualize_batch"
            )
        PreprocessingStageOrderValidator().run(
            stage_order=self.stage_order,
            batch_correction_requested=(
                batch_correction_method != DATASET_BATCH_CORRECTION_METHOD_NONE
            ),
        )
        validate_group_coverage_filter_config(
            enabled=self.group_coverage_filter_enabled,
            group_column=self.group_coverage_filter_group_column,
            min_finite_observations_per_group=(
                self.group_coverage_filter_min_finite_observations_per_group
            ),
            min_finite_fraction_per_group=(
                self.group_coverage_filter_min_finite_fraction_per_group
            ),
            min_groups_passing_threshold=(
                self.group_coverage_filter_min_groups_passing_threshold
            ),
        )
        if (
            self.group_coverage_filter_enabled
            and DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER
            not in self.stage_order
        ):
            raise PhosPyInputError(
                "dataset preprocessing plan requests group-aware coverage "
                "filtering but stage_order does not include "
                "'group_coverage_filter'. Build plans from "
                "DatasetPreprocessingConfig or include the group_coverage_filter "
                "stage explicitly."
            )
        object.__setattr__(
            self,
            "stage_order_resolution",
            _normalize_stage_order_resolution(
                stage_order=self.stage_order,
                stage_order_resolution=self.stage_order_resolution,
            ),
        )

    @classmethod
    def from_config(cls, config: DatasetPreprocessingConfig) -> PreprocessingPlan:
        return PreprocessingPlanInterpreter(plan_type=cls).run(config)

    @classmethod
    def default(cls) -> PreprocessingPlan:
        return cls.from_config(DatasetPreprocessingConfig())


class PreprocessingPlanInterpreter:
    """Convert public preprocessing config into an execution-ready plan."""

    def __init__(
        self,
        *,
        plan_type: type[PreprocessingPlan] = PreprocessingPlan,
    ) -> None:
        self._plan_type = plan_type

    def run(self, config: DatasetPreprocessingConfig) -> PreprocessingPlan:
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

        stage_order: list[str] = []
        stage_order_resolution: list[PreprocessingStageOrderResolution] = []

        def _append_stage(stage: str, *, rationale: str) -> None:
            stage_order.append(stage)
            stage_order_resolution.append(
                PreprocessingStageOrderResolution(
                    stage=stage,
                    order_index=len(stage_order) - 1,
                    rationale=rationale,
                )
            )

        site_sequence_resolution_enabled = (
            config.site_sequence_resolution.fasta_path is not None
        )
        intensity_transform_policy = IntensityTransformPolicy.parse(
            config.intensity_transform.policy,
            field_name="preprocessing_config.intensity_transform.policy",
        )
        normalisation_policy = NormalisationPolicy.parse(
            config.normalisation.policy,
            field_name="preprocessing_config.normalisation.policy",
        )
        site_matrix_policy = SiteMatrixPolicy.parse(
            config.site_matrix.policy,
            field_name="preprocessing_config.site_matrix.policy",
        )
        site_matrix_duplicate_site_policy = SiteMatrixDuplicateSitePolicy.parse(
            config.site_matrix.duplicate_site_policy,
            field_name="preprocessing_config.site_matrix.duplicate_site_policy",
        )
        site_matrix_missing_data_policy = SiteMatrixMissingDataPolicy.parse(
            config.site_matrix.missing_data_policy,
            field_name="preprocessing_config.site_matrix.missing_data_policy",
        )
        comparison_building_policy = ComparisonBuildingPolicy.parse(
            config.comparisons.policy,
            field_name="preprocessing_config.comparisons.policy",
        )
        site_sequence_resolution_mode = SiteSequenceResolutionMode.parse(
            config.site_sequence_resolution.mode,
            field_name="preprocessing_config.site_sequence_resolution.mode",
        )
        if site_sequence_resolution_enabled:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
            )
        localisation_mode = LocalisationEligibilityMode.parse(
            config.localisation.mode,
            field_name="preprocessing_config.localisation.mode",
        )
        if localisation_mode is not LocalisationEligibilityMode.IGNORE:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_LOCALISATION,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
            )
        missing_data_policy = MissingDataPolicy.parse(
            config.missing_data.policy,
            field_name="preprocessing_config.missing_data.policy",
        )
        batch_correction_method = str(config.batch_correction.method).strip()
        batch_correction_condition_column = _resolve_batch_condition_column(config)
        batch_correction_condition_columns = _resolve_batch_condition_columns(config)
        batch_correction_internal_request = _resolve_batch_internal_request(config)
        batch_correction_control_site_set = _resolve_batch_control_site_set(config)
        batch_correction_missingness_policy = _resolve_batch_missingness_policy(config)
        batch_correction_replicate_column = _resolve_batch_replicate_column(config)
        total_correction_policy = TotalProteinCorrectionPolicy.parse(
            config.total_protein_correction.policy,
            field_name="preprocessing_config.total_protein_correction.policy",
        )
        if config.group_coverage_filter.enabled:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER,
            )
        if missing_data_policy is MissingDataPolicy.IMPUTE_MINPROB:
            if intensity_transform_policy is not IntensityTransformPolicy.LOG2:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.policy="
                    "'impute_minprob' requires "
                    "preprocessing_config.intensity_transform.policy='log2'. "
                    "Set intensity_transform.policy='log2' or choose a different "
                    "missing_data policy."
                )
            _append_stage(
                DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
                rationale=(
                    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM
                ),
            )
            _append_stage(
                DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
            )
        else:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_MISSING_DATA,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
            )
            if intensity_transform_policy is not IntensityTransformPolicy.IDENTITY:
                _append_stage(
                    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
                    rationale=(
                        PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM
                    ),
                )
        if batch_correction_method != DATASET_BATCH_CORRECTION_METHOD_NONE:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_BATCH_CORRECTION,
            )
        if total_correction_policy is not TotalProteinCorrectionPolicy.NONE:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
            )
        if site_matrix_policy is not SiteMatrixPolicy.AS_INPUT:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
            )
        if normalisation_policy is not NormalisationPolicy.NONE:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_NORMALISATION,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
            )
        if comparison_building_policy is not ComparisonBuildingPolicy.NONE:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_COMPARISONS,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
            )
        return self._plan_type(
            intensity_transform_policy=intensity_transform_policy,
            intensity_transform_pseudocount=float(
                config.intensity_transform.pseudocount
            ),
            normalisation_policy=normalisation_policy,
            missing_data_policy=missing_data_policy,
            missing_data_min_observed_values=config.missing_data.min_observed_values,
            missing_data_q=(
                None if config.missing_data.q is None else float(config.missing_data.q)
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
            localisation_mode=localisation_mode,
            localisation_min_confidence=float(config.localisation.min_confidence),
            localisation_confidence_column=str(
                config.localisation.confidence_column
            ).strip(),
            localisation_waiver_reason=(
                None
                if config.localisation.waiver_reason is None
                else str(config.localisation.waiver_reason).strip()
            ),
            site_sequence_resolution_enabled=site_sequence_resolution_enabled,
            site_sequence_resolution_fasta_path=(
                config.site_sequence_resolution.fasta_path
            ),
            site_sequence_resolution_mode=site_sequence_resolution_mode,
            site_sequence_resolution_conflict_policy=(
                _resolve_site_sequence_resolution_conflict_policy(
                    mode=site_sequence_resolution_mode,
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
            group_coverage_filter_enabled=bool(config.group_coverage_filter.enabled),
            group_coverage_filter_group_column=(
                config.group_coverage_filter.group_column
            ),
            group_coverage_filter_min_finite_observations_per_group=(
                config.group_coverage_filter.min_finite_observations_per_group
            ),
            group_coverage_filter_min_finite_fraction_per_group=(
                None
                if config.group_coverage_filter.min_finite_fraction_per_group is None
                else float(config.group_coverage_filter.min_finite_fraction_per_group)
            ),
            group_coverage_filter_min_groups_passing_threshold=(
                config.group_coverage_filter.min_groups_passing_threshold
            ),
            total_protein_correction_policy=total_correction_policy,
            total_protein_correction_identity_policy=(
                _resolve_total_correction_identity_policy(
                    config.total_protein_correction.identity
                )
            ),
            protein_aware_preparation_policy=config.protein_aware_preparation.policy,
            protein_aware_preparation_mapping_policy=(
                config.protein_aware_preparation.protein_mapping_policy
            ),
            site_matrix_policy=site_matrix_policy,
            site_matrix_duplicate_site_policy=site_matrix_duplicate_site_policy,
            site_matrix_missing_data_policy=site_matrix_missing_data_policy,
            site_matrix_minimum_observed_values=(
                config.site_matrix.minimum_observed_values
            ),
            comparison_building_policy=comparison_building_policy,
            comparison_sample_group_column=config.comparisons.sample_group_column,
            comparison_pairs=(
                None
                if config.comparisons.pairs is None
                else tuple(config.comparisons.pairs)
            ),
            ruv_readiness_enabled=bool(config.ruv_readiness.enabled),
            ruv_readiness_control_feature_column=(
                config.ruv_readiness.control_feature_column
            ),
            ruv_readiness_replicate_group_column=(
                config.ruv_readiness.replicate_group_column
            ),
            ruv_readiness_batch_column=config.ruv_readiness.batch_column,
            batch_correction_method=batch_correction_method,
            batch_correction_batch_column=config.batch_correction.batch_column,
            batch_correction_condition_column=batch_correction_condition_column,
            batch_correction_condition_columns=batch_correction_condition_columns,
            batch_correction_replicate_column=batch_correction_replicate_column,
            batch_correction_control_site_set=batch_correction_control_site_set,
            batch_correction_missingness_policy=batch_correction_missingness_policy,
            batch_correction_internal_request=batch_correction_internal_request,
            batch_correction_preserve_condition_effects=(
                _resolve_batch_preserve_condition_effects(config)
            ),
            stage_order=tuple(stage_order),
            stage_order_resolution=tuple(stage_order_resolution),
        )


def _resolve_total_correction_identity_policy(
    config: DatasetTotalProteinCorrectionIdentityConfig,
) -> TotalProteinCorrectionIdentityPolicy:
    if config.mode == DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT:
        return TotalProteinCorrectionIdentityPolicy(
            mode=config.mode,
            phosphosite_key=str(config.phosphosite_key).strip(),
            total_protein_key=str(config.total_protein_key).strip(),
            matching_policy=TotalProteinCorrectionIdentityMatchingPolicy.parse(
                config.matching_policy,
                field_name=(
                    "preprocessing_config.total_protein_correction.identity."
                    "matching_policy"
                ),
            ),
            duplicate_policy=config.duplicate_policy,
            unmatched_policy=config.unmatched_policy,
            mapping_table=None,
            mapping_phosphosite_key=None,
            mapping_total_protein_key=None,
            mapping_table_fingerprint=None,
        )

    if config.mode != DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity contains an unsupported mode"
        )
    mapping_table = config.mapping_table
    if mapping_table is None:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.mapping_table is required when identity.mode='mapping_table'"
        )

    mapping_phosphosite_key = str(config.mapping_phosphosite_key).strip()
    mapping_total_protein_key = str(config.mapping_total_protein_key).strip()
    if mapping_phosphosite_key not in mapping_table.columns:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.mapping_table is missing column "
            f"{mapping_phosphosite_key!r}"
        )
    if mapping_total_protein_key not in mapping_table.columns:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.total_protein_correction."
            "identity.mapping_table is missing column "
            f"{mapping_total_protein_key!r}"
        )
    normalized_table = pd.DataFrame(
        {
            "phosphosite_id": mapping_table.loc[:, mapping_phosphosite_key]
            .astype("string")
            .str.strip(),
            "total_protein_id": mapping_table.loc[:, mapping_total_protein_key]
            .astype("string")
            .str.strip(),
        }
    )

    def _is_missing_mapping_value(value: object) -> bool:
        return bool(pd.Series((value,), dtype="object").isna().iat[0])

    mapping_rows = tuple(
        (
            ""
            if _is_missing_mapping_value(record.get("phosphosite_id"))
            else str(record.get("phosphosite_id")),
            ""
            if _is_missing_mapping_value(record.get("total_protein_id"))
            else str(record.get("total_protein_id")),
        )
        for record in normalized_table.to_dict(orient="records")
    )
    fingerprint_table = (
        normalized_table.fillna("<MISSING>")
        .sort_values(by=["phosphosite_id", "total_protein_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    return TotalProteinCorrectionIdentityPolicy(
        mode=config.mode,
        phosphosite_key=str(config.phosphosite_key).strip(),
        total_protein_key=str(config.total_protein_key).strip(),
        matching_policy=TotalProteinCorrectionIdentityMatchingPolicy.parse(
            config.matching_policy,
            field_name=(
                "preprocessing_config.total_protein_correction.identity.matching_policy"
            ),
        ),
        duplicate_policy=config.duplicate_policy,
        unmatched_policy=config.unmatched_policy,
        mapping_table=mapping_rows,
        mapping_phosphosite_key=mapping_phosphosite_key,
        mapping_total_protein_key=mapping_total_protein_key,
        mapping_table_fingerprint=hash_table_tolerance(
            fingerprint_table,
            name="total_protein_correction.identity.mapping_table",
        ),
    )


def _resolve_site_sequence_resolution_conflict_policy(
    *,
    mode: SiteSequenceResolutionMode,
    conflict_policy: DatasetSiteSequenceConflictPolicy | None,
) -> SiteSequenceConflictPolicy:
    if conflict_policy is not None:
        return SiteSequenceConflictPolicy.parse(
            conflict_policy,
            field_name="preprocessing_config.site_sequence_resolution.conflict_policy",
        )
    if mode is SiteSequenceResolutionMode.REPLACE_EXISTING:
        return SiteSequenceConflictPolicy.REPLACE_EXISTING
    return SiteSequenceConflictPolicy.PRESERVE_EXISTING


def _sps_ruv_config(
    config: DatasetPreprocessingConfig,
) -> SpsRuvBatchCorrectionConfig | None:
    batch_correction = config.batch_correction
    if isinstance(batch_correction, SpsRuvBatchCorrectionConfig):
        return batch_correction
    return None


def _resolve_batch_condition_column(config: DatasetPreprocessingConfig) -> str:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.condition_columns[0]
    batch_correction = cast(DatasetBatchCorrectionConfig, config.batch_correction)
    return batch_correction.condition_column


def _resolve_batch_condition_columns(
    config: DatasetPreprocessingConfig,
) -> tuple[str, ...]:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.condition_columns
    batch_correction = cast(DatasetBatchCorrectionConfig, config.batch_correction)
    return (batch_correction.condition_column,)


def _resolve_batch_replicate_column(config: DatasetPreprocessingConfig) -> str | None:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.replicate_column
    return None


def _resolve_batch_internal_request(
    config: DatasetPreprocessingConfig,
) -> InternalBatchCorrectionRequest | None:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.to_internal_request()
    return None


def _resolve_batch_control_site_set(
    config: DatasetPreprocessingConfig,
) -> object | None:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.control_site_set
    return None


def _resolve_batch_missingness_policy(
    config: DatasetPreprocessingConfig,
) -> CorrectionMissingnessPolicy | None:
    sps_config = _sps_ruv_config(config)
    if sps_config is not None:
        return sps_config.missingness_policy
    return None


def _resolve_batch_preserve_condition_effects(
    config: DatasetPreprocessingConfig,
) -> bool:
    if _sps_ruv_config(config) is not None:
        return True
    batch_correction = cast(DatasetBatchCorrectionConfig, config.batch_correction)
    return cast(bool, batch_correction.preserve_condition_effects)


def _normalize_stage_order_resolution(
    *,
    stage_order: tuple[str, ...],
    stage_order_resolution: tuple[PreprocessingStageOrderResolution, ...],
) -> tuple[PreprocessingStageOrderResolution, ...]:
    if not stage_order:
        return ()
    if (
        len(stage_order_resolution) == len(stage_order)
        and tuple(item.stage for item in stage_order_resolution) == stage_order
    ):
        return tuple(
            PreprocessingStageOrderResolution(
                stage=stage,
                order_index=index,
                rationale=(
                    str(stage_order_resolution[index].rationale).strip()
                    or PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE
                ),
            )
            for index, stage in enumerate(stage_order)
        )
    return tuple(
        PreprocessingStageOrderResolution(
            stage=stage,
            order_index=index,
            rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
        )
        for index, stage in enumerate(stage_order)
    )
