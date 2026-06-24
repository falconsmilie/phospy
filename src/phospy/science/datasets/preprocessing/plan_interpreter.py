"""Interpret public preprocessing config into internal preprocessing plans."""

from __future__ import annotations

from typing import cast

import pandas as pd

from phospy.contracts.configs import (
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    DATASET_PROTEIN_AWARE_PREPARATION_MAPPING_POLICIES,
    DATASET_PROTEIN_AWARE_PREPARATION_POLICIES,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_DIRECT,
    DATASET_TOTAL_PROTEIN_CORRECTION_IDENTITY_MODE_MAPPING_TABLE,
    CorrectionMissingnessPolicy,
    DatasetBatchCorrectionConfig,
    DatasetPreprocessingConfig,
    DatasetSiteSequenceConflictPolicy,
    DatasetTotalProteinCorrectionIdentityConfig,
    SpsRuvBatchCorrectionConfig,
)
from phospy.contracts.configs.preprocessing import InternalBatchCorrectionRequest
from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table_tolerance
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_LOCALISATION,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
    DATASET_PREPROCESSING_STAGE_NORMALISATION,
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
    PreprocessingPlan,
    PreprocessingStageOrderResolution,
    TotalProteinCorrectionIdentityPolicy,
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
from phospy.validation.configs.preprocessing import (
    validate_protein_aware_preparation_config,
)


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
            site_sequence_resolution_fasta_path=config.site_sequence_resolution.fasta_path,
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
            total_protein_correction_identity_policy=_resolve_total_correction_identity_policy(
                config.total_protein_correction.identity
            ),
            protein_aware_preparation_policy=config.protein_aware_preparation.policy,
            protein_aware_preparation_mapping_policy=(
                config.protein_aware_preparation.protein_mapping_policy
            ),
            site_matrix_policy=site_matrix_policy,
            site_matrix_duplicate_site_policy=site_matrix_duplicate_site_policy,
            site_matrix_missing_data_policy=site_matrix_missing_data_policy,
            site_matrix_minimum_observed_values=config.site_matrix.minimum_observed_values,
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


__all__ = ["PreprocessingPlanInterpreter"]
