"""Stage-order validation and planning for dataset preprocessing."""

from __future__ import annotations

from dataclasses import dataclass

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
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_INTENSITY_TRANSFORM,
    PREPROCESSING_STAGE_ORDER_RATIONALE_MINPROB_MISSING_DATA,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_INTENSITY_TRANSFORM,
    PREPROCESSING_STAGE_ORDER_RATIONALE_NON_MINPROB_MISSING_DATA,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ComparisonBuildingPolicy,
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
    ) -> PreprocessingStageOrderPlan:
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
        return PreprocessingStageOrderPlan(
            stage_order=tuple(stage_order),
            stage_order_resolution=tuple(stage_order_resolution),
        )


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
