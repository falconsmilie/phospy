"""Stage-order validation and planning for dataset preprocessing."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from phospy.errors.input import PhosPyInputError
from phospy.science.configs.preprocessing import DATASET_BATCH_CORRECTION_METHOD_NONE
from phospy.science.datasets.preprocessing.plan_constants import (
    BATCH_CORRECTION_DOWNSTREAM_BOUNDARY_STAGES,
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
    EXTERNAL_CORRECTION_DOWNSTREAM_MATRIX_CONSUMING_STAGES,
    PREPROCESSING_STAGE_ORDER_RATIONALE_BATCH_CORRECTION,
    PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
    PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER,
    PREPROCESSING_STAGE_ORDER_RATIONALE_LOG2_IMPUTATION_INTENSITY_TRANSFORM,
    PREPROCESSING_STAGE_ORDER_RATIONALE_LOG2_IMPUTATION_MISSING_DATA,
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM,
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
    ImputationInputScale,
    IntensityTransformPolicy,
    LocalisationEligibilityMode,
    MissingDataPolicy,
    NormalisationPolicy,
    SiteMatrixPolicy,
    TotalProteinCorrectionPolicy,
)


@dataclass(frozen=True, slots=True)
class PreprocessingStageOrderResolution:
    """Structured explanation of resolved preprocessing stage order."""

    stage: str
    order_index: int
    rationale: str


@dataclass(frozen=True, slots=True)
class PreprocessingStageOrderPlan:
    """Resolved preprocessing stage order and deterministic rationale rows."""

    stage_order: tuple[str, ...]
    stage_order_resolution: tuple[PreprocessingStageOrderResolution, ...]


class _StageAppender(Protocol):
    def __call__(self, stage: str, *, rationale: str) -> None: ...


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
            for stage in BATCH_CORRECTION_DOWNSTREAM_BOUNDARY_STAGES
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


class PreprocessingStageOrderPlanner:
    """Resolve configured preprocessing policies into ordered execution stages."""

    def run(
        self,
        *,
        site_sequence_resolution_enabled: bool,
        intensity_transform_policy: IntensityTransformPolicy,
        normalisation_policy: NormalisationPolicy,
        site_matrix_policy: SiteMatrixPolicy,
        comparison_building_policy: ComparisonBuildingPolicy,
        localisation_mode: LocalisationEligibilityMode,
        missing_data_policy: MissingDataPolicy,
        batch_correction_method: str,
        total_correction_policy: TotalProteinCorrectionPolicy,
        group_coverage_filter_enabled: bool,
        missing_data_input_scale: ImputationInputScale | None = None,
    ) -> PreprocessingStageOrderPlan:
        stage_order: list[str] = []
        stage_order_resolution: list[PreprocessingStageOrderResolution] = []
        missing_data_policy = MissingDataPolicy.parse(
            missing_data_policy,
            field_name="dataset preprocessing stage-order missing_data_policy",
        )
        intensity_transform_policy = IntensityTransformPolicy.parse(
            intensity_transform_policy,
            field_name=("dataset preprocessing stage-order intensity_transform_policy"),
        )

        def _append_stage(stage: str, *, rationale: str) -> None:
            stage_order.append(stage)
            stage_order_resolution.append(
                PreprocessingStageOrderResolution(
                    stage=stage,
                    order_index=len(stage_order) - 1,
                    rationale=rationale,
                )
            )

        if site_sequence_resolution_enabled:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_SITE_SEQUENCE_RESOLUTION,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
            )
        if localisation_mode is not LocalisationEligibilityMode.IGNORE:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_LOCALISATION,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_CONFIGURED_STAGE,
            )
        if group_coverage_filter_enabled:
            _append_stage(
                DATASET_PREPROCESSING_STAGE_GROUP_COVERAGE_FILTER,
                rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_GROUP_COVERAGE_FILTER,
            )
        _append_missing_data_and_transform_stages(
            append_stage=_append_stage,
            missing_data_policy=missing_data_policy,
            missing_data_input_scale=missing_data_input_scale,
            intensity_transform_policy=intensity_transform_policy,
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
        return PreprocessingStageOrderPlan(
            stage_order=tuple(stage_order),
            stage_order_resolution=tuple(stage_order_resolution),
        )


def _append_missing_data_and_transform_stages(
    *,
    append_stage: _StageAppender,
    missing_data_policy: MissingDataPolicy,
    missing_data_input_scale: ImputationInputScale | None,
    intensity_transform_policy: IntensityTransformPolicy,
) -> None:
    if missing_data_policy is MissingDataPolicy.FORBID:
        _append_linear_or_strict_missing_data_path(
            append_stage=append_stage,
            intensity_transform_policy=intensity_transform_policy,
        )
        return
    input_scale = _require_imputation_input_scale(
        missing_data_policy=missing_data_policy,
        missing_data_input_scale=missing_data_input_scale,
    )
    if input_scale is ImputationInputScale.LOG2:
        _append_log2_imputation_path(
            append_stage=append_stage,
            missing_data_policy=missing_data_policy,
            intensity_transform_policy=intensity_transform_policy,
        )
        return
    if missing_data_policy is MissingDataPolicy.IMPUTE_MINPROB:
        raise PhosPyInputError(
            "dataset preprocessing plan has unsupported stage-order policy: "
            "missing_data.policy='impute_minprob' cannot run on "
            "missing_data.input_scale='linear'"
        )
    _append_linear_or_strict_missing_data_path(
        append_stage=append_stage,
        intensity_transform_policy=intensity_transform_policy,
    )


def _append_linear_or_strict_missing_data_path(
    *,
    append_stage: _StageAppender,
    intensity_transform_policy: IntensityTransformPolicy,
) -> None:
    append_stage(
        DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
    )
    if intensity_transform_policy is not IntensityTransformPolicy.IDENTITY:
        append_stage(
            DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
            rationale=PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM,
        )


def _append_log2_imputation_path(
    *,
    append_stage: _StageAppender,
    missing_data_policy: MissingDataPolicy,
    intensity_transform_policy: IntensityTransformPolicy,
) -> None:
    if intensity_transform_policy is IntensityTransformPolicy.LOG2:
        append_stage(
            DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
            rationale=_log2_transform_rationale(missing_data_policy),
        )
    append_stage(
        DATASET_PREPROCESSING_STAGE_MISSING_DATA,
        rationale=_log2_missing_data_rationale(missing_data_policy),
    )


def _require_imputation_input_scale(
    *,
    missing_data_policy: MissingDataPolicy,
    missing_data_input_scale: ImputationInputScale | None,
) -> ImputationInputScale:
    if missing_data_input_scale is None:
        raise PhosPyInputError(
            "dataset preprocessing stage-order planning requires "
            "missing_data.input_scale for imputation policy "
            f"{missing_data_policy.value!r}"
        )
    return ImputationInputScale.parse(
        missing_data_input_scale,
        field_name="dataset preprocessing plan missing_data_input_scale",
    )


def _log2_transform_rationale(missing_data_policy: MissingDataPolicy) -> str:
    if missing_data_policy is MissingDataPolicy.IMPUTE_MINPROB:
        return PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM
    return PREPROCESSING_STAGE_ORDER_RATIONALE_LOG2_IMPUTATION_INTENSITY_TRANSFORM


def _log2_missing_data_rationale(missing_data_policy: MissingDataPolicy) -> str:
    if missing_data_policy is MissingDataPolicy.IMPUTE_MINPROB:
        return PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA
    return PREPROCESSING_STAGE_ORDER_RATIONALE_LOG2_IMPUTATION_MISSING_DATA


def reject_external_corrected_output_after_downstream_preprocessing(
    stage_order: tuple[str, ...],
) -> None:
    """Reject external correction when the plan has downstream matrix consumers."""

    configured = tuple(
        stage
        for stage in (str(item).strip() for item in stage_order)
        if stage in EXTERNAL_CORRECTION_DOWNSTREAM_MATRIX_CONSUMING_STAGES
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


def normalize_stage_order_resolution(
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


__all__ = [
    "PreprocessingStageOrderPlan",
    "PreprocessingStageOrderPlanner",
    "PreprocessingStageOrderResolution",
    "PreprocessingStageOrderValidator",
    "normalize_stage_order_resolution",
    "reject_external_corrected_output_after_downstream_preprocessing",
]
