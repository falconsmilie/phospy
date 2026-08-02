"""Preprocessing stage-order planning and config interpretation."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Any

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table_tolerance
from phospy.science.configs.dataset import DatasetPreprocessingConfig
from phospy.science.configs.preprocessing import (
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
    CorrectionMissingnessPolicy,
    DatasetComparisonPair,
    DatasetProteinAwarePreparationMappingPolicy,
    DatasetProteinAwarePreparationPolicy,
    DatasetTotalProteinCorrectionDuplicatePolicy,
    DatasetTotalProteinCorrectionIdentityConfig,
    DatasetTotalProteinCorrectionIdentityMode,
    DatasetTotalProteinCorrectionUnmatchedPolicy,
    InternalBatchCorrectionRequest,
)
from phospy.science.configs.preprocessing._validation import (
    validate_protein_aware_preparation_config,
)
from phospy.science.datasets.preprocessing.imputation_scale_policy import (
    reject_incompatible_imputation_stage_order,
    reject_unestablished_log2_imputation_input_scale,
    resolve_imputation_scale_policy,
)
from phospy.science.datasets.preprocessing.plan_config_resolution import (
    PreprocessingConfigPolicyResolver,
    resolve_site_sequence_resolution_conflict_policy,
)
from phospy.science.datasets.preprocessing.plan_constants import (
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
from phospy.science.datasets.preprocessing.plan_rules import (
    PreprocessingBatchCorrectionPlanRuleFamily,
    PreprocessingCorePlanPolicyRuleFamily,
    PreprocessingDownstreamPlanPolicyRuleFamily,
    PreprocessingGroupCoveragePlanRuleFamily,
    PreprocessingLocalisationPlanRuleFamily,
)
from phospy.science.datasets.preprocessing.plan_stage_order import (
    PreprocessingStageOrderPlanner,
    PreprocessingStageOrderResolution,
    PreprocessingStageOrderValidator,
    normalize_stage_order_resolution,
    reject_external_corrected_output_after_downstream_preprocessing,
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
    TotalProteinCorrectionIdentityMatchingPolicy,
    TotalProteinCorrectionPolicy,
)
from phospy.science.transformations.models import IntensityScaleKind


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


def _set_resolved_plan_fields(plan: object, resolved_fields: Any) -> None:
    for field in fields(resolved_fields):
        object.__setattr__(plan, field.name, getattr(resolved_fields, field.name))


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
        core = PreprocessingCorePlanPolicyRuleFamily().run(
            intensity_transform_policy=self.intensity_transform_policy,
            normalisation_policy=self.normalisation_policy,
            missing_data_policy=self.missing_data_policy,
        )
        _set_resolved_plan_fields(self, core)
        imputation_scale = resolve_imputation_scale_policy(
            missing_data_policy=core.missing_data_policy,
            requested_input_scale=self.missing_data_input_scale,
            intensity_transform_policy=core.intensity_transform_policy,
        )
        object.__setattr__(
            self,
            "missing_data_input_scale",
            imputation_scale.input_scale,
        )
        object.__setattr__(
            self,
            "missing_data_input_scale_source",
            imputation_scale.input_scale_source,
        )
        object.__setattr__(
            self,
            "missing_data_imputation_operation_order",
            imputation_scale.operation_order,
        )
        reject_incompatible_imputation_stage_order(
            stage_order=self.stage_order,
            resolved_policy=imputation_scale,
        )

        localisation = PreprocessingLocalisationPlanRuleFamily().run(
            localisation_mode=self.localisation_mode,
            localisation_min_confidence=self.localisation_min_confidence,
            localisation_confidence_column=self.localisation_confidence_column,
            localisation_waiver_reason=self.localisation_waiver_reason,
        )
        _set_resolved_plan_fields(self, localisation)

        downstream = PreprocessingDownstreamPlanPolicyRuleFamily().run(
            site_sequence_resolution_mode=self.site_sequence_resolution_mode,
            site_sequence_resolution_conflict_policy=(
                self.site_sequence_resolution_conflict_policy
            ),
            site_matrix_policy=self.site_matrix_policy,
            comparison_building_policy=self.comparison_building_policy,
            site_matrix_duplicate_site_policy=self.site_matrix_duplicate_site_policy,
            site_matrix_missing_data_policy=self.site_matrix_missing_data_policy,
            total_protein_correction_policy=self.total_protein_correction_policy,
            protein_aware_preparation_policy=self.protein_aware_preparation_policy,
            protein_aware_preparation_mapping_policy=(
                self.protein_aware_preparation_mapping_policy
            ),
        )
        _set_resolved_plan_fields(self, downstream)

        batch_correction = PreprocessingBatchCorrectionPlanRuleFamily().run(
            batch_correction_method=self.batch_correction_method,
            batch_correction_batch_column=self.batch_correction_batch_column,
            batch_correction_condition_column=self.batch_correction_condition_column,
            batch_correction_condition_columns=self.batch_correction_condition_columns,
            batch_correction_replicate_column=self.batch_correction_replicate_column,
            batch_correction_control_site_set=self.batch_correction_control_site_set,
            batch_correction_missingness_policy=(
                self.batch_correction_missingness_policy
            ),
            batch_correction_internal_request=self.batch_correction_internal_request,
            batch_correction_preserve_condition_effects=(
                self.batch_correction_preserve_condition_effects
            ),
        )
        _set_resolved_plan_fields(self, batch_correction)
        PreprocessingStageOrderValidator().run(
            stage_order=self.stage_order,
            batch_correction_requested=(
                batch_correction.batch_correction_method
                != DATASET_BATCH_CORRECTION_METHOD_NONE
            ),
        )
        PreprocessingGroupCoveragePlanRuleFamily().run(
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
            stage_order=self.stage_order,
        )
        object.__setattr__(
            self,
            "stage_order_resolution",
            normalize_stage_order_resolution(
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

    def run(
        self,
        config: DatasetPreprocessingConfig,
        *,
        declared_input_scale_kind: IntensityScaleKind | None = None,
    ) -> PreprocessingPlan:
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
        imputation_scale = resolve_imputation_scale_policy(
            missing_data_policy=policies.missing_data_policy,
            requested_input_scale=config.missing_data.input_scale,
            intensity_transform_policy=policies.intensity_transform_policy,
        )
        reject_unestablished_log2_imputation_input_scale(
            resolved_policy=imputation_scale,
            intensity_transform_policy=policies.intensity_transform_policy,
            declared_input_scale_kind=declared_input_scale_kind,
        )
        stage_plan = PreprocessingStageOrderPlanner().run(
            site_sequence_resolution_enabled=(
                policies.site_sequence_resolution_enabled
            ),
            intensity_transform_policy=policies.intensity_transform_policy,
            normalisation_policy=policies.normalisation_policy,
            site_matrix_policy=policies.site_matrix_policy,
            comparison_building_policy=policies.comparison_building_policy,
            localisation_mode=policies.localisation_mode,
            missing_data_policy=policies.missing_data_policy,
            missing_data_input_scale=imputation_scale.input_scale,
            batch_correction_method=policies.batch_correction_method,
            total_correction_policy=policies.total_correction_policy,
            group_coverage_filter_enabled=bool(config.group_coverage_filter.enabled),
        )
        return self._plan_type(
            intensity_transform_policy=policies.intensity_transform_policy,
            intensity_transform_pseudocount=float(
                config.intensity_transform.pseudocount
            ),
            normalisation_policy=policies.normalisation_policy,
            missing_data_policy=policies.missing_data_policy,
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
            missing_data_input_scale=imputation_scale.input_scale,
            missing_data_input_scale_source=imputation_scale.input_scale_source,
            missing_data_imputation_operation_order=(imputation_scale.operation_order),
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
            total_protein_correction_policy=policies.total_correction_policy,
            total_protein_correction_identity_policy=(
                _resolve_total_correction_identity_policy(
                    config.total_protein_correction.identity
                )
            ),
            protein_aware_preparation_policy=config.protein_aware_preparation.policy,
            protein_aware_preparation_mapping_policy=(
                config.protein_aware_preparation.protein_mapping_policy
            ),
            site_matrix_policy=policies.site_matrix_policy,
            site_matrix_duplicate_site_policy=(
                policies.site_matrix_duplicate_site_policy
            ),
            site_matrix_missing_data_policy=policies.site_matrix_missing_data_policy,
            site_matrix_minimum_observed_values=(
                config.site_matrix.minimum_observed_values
            ),
            comparison_building_policy=policies.comparison_building_policy,
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
            stage_order=stage_plan.stage_order,
            stage_order_resolution=stage_plan.stage_order_resolution,
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
