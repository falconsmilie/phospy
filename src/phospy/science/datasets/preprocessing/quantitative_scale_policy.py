"""Private quantitative-scale policy for additive preprocessing operations."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from typing import Final

from phospy.errors.input import PhosPyInputError
from phospy.science.configs import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    DATASET_BATCH_CORRECTION_METHOD_SPS_RUV_STYLE,
    SPS_RUV_BATCH_CORRECTION_METHODS,
)
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.science.datasets.preprocessing.plan import PreprocessingPlan
from phospy.science.datasets.preprocessing.policy_models import (
    IntensityTransformPolicy,
    NormalisationPolicy,
)
from phospy.science.transformations.models import IntensityScaleKind

SUPPORTED_ADDITIVE_PREPROCESSING_SCALE_KINDS: Final[frozenset[IntensityScaleKind]] = (
    frozenset({IntensityScaleKind.LOG2})
)
_REQUIRED_ADDITIVE_SCALE = IntensityScaleKind.LOG2


class AdditivePreprocessingOperation(str, Enum):
    """Structured additive preprocessing operations requiring scale policy."""

    MEDIAN_CENTER = NormalisationPolicy.MEDIAN_CENTER.value
    LINEAR_RESIDUALIZE_BATCH = DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH
    SPS_RUV_STYLE = DATASET_BATCH_CORRECTION_METHOD_SPS_RUV_STYLE


@dataclass(frozen=True, slots=True)
class AdditivePreprocessingScaleGuard:
    """Reject additive preprocessing unless the current scale is supported."""

    def run(
        self,
        *,
        preprocessing_plan: PreprocessingPlan,
        declared_input_scale_kind: IntensityScaleKind | None,
        corrected_preprocessing_output: CorrectedPreprocessingOutput | None = None,
    ) -> None:
        current_scale_kind = _resolve_downstream_scale_kind(
            preprocessing_plan=preprocessing_plan,
            declared_input_scale_kind=declared_input_scale_kind,
        )
        for operation in _iter_additive_operations(
            preprocessing_plan=preprocessing_plan,
            corrected_preprocessing_output=corrected_preprocessing_output,
        ):
            _require_supported_additive_scale(
                operation=operation,
                current_scale_kind=current_scale_kind,
            )


def _iter_additive_operations(
    *,
    preprocessing_plan: PreprocessingPlan,
    corrected_preprocessing_output: CorrectedPreprocessingOutput | None,
) -> Iterable[AdditivePreprocessingOperation]:
    if preprocessing_plan.normalisation_policy is NormalisationPolicy.MEDIAN_CENTER:
        yield AdditivePreprocessingOperation.MEDIAN_CENTER

    method = preprocessing_plan.batch_correction_method
    if method == DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH:
        yield AdditivePreprocessingOperation.LINEAR_RESIDUALIZE_BATCH
    elif method in SPS_RUV_BATCH_CORRECTION_METHODS:
        yield AdditivePreprocessingOperation.SPS_RUV_STYLE

    if (
        corrected_preprocessing_output is None
        or method != DATASET_BATCH_CORRECTION_METHOD_NONE
    ):
        return
    corrected_method = corrected_preprocessing_output.batch_correction_report.method
    if corrected_method == DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH:
        yield AdditivePreprocessingOperation.LINEAR_RESIDUALIZE_BATCH
    elif corrected_method in SPS_RUV_BATCH_CORRECTION_METHODS:
        yield AdditivePreprocessingOperation.SPS_RUV_STYLE


def _resolve_downstream_scale_kind(
    *,
    preprocessing_plan: PreprocessingPlan,
    declared_input_scale_kind: IntensityScaleKind | None,
) -> IntensityScaleKind | None:
    if preprocessing_plan.intensity_transform_policy is IntensityTransformPolicy.LOG2:
        return IntensityScaleKind.LOG2
    return declared_input_scale_kind


def _require_supported_additive_scale(
    *,
    operation: AdditivePreprocessingOperation,
    current_scale_kind: IntensityScaleKind | None,
) -> None:
    if current_scale_kind in SUPPORTED_ADDITIVE_PREPROCESSING_SCALE_KINDS:
        return
    raise PhosPyInputError(
        "dataset build request additive preprocessing scale policy rejected "
        f"operation='{operation.value}'; "
        f"current_scale='{_scale_label(current_scale_kind)}'; "
        f"required_scale='{_REQUIRED_ADDITIVE_SCALE.value}'; "
        "required_state='established log2 phosphosite abundance'; "
        'corrective_action="configure '
        "preprocessing_config.intensity_transform.policy='log2' before this "
        "additive operation, or provide already-log2 abundance data with "
        "input_intensity_scale='log2'. Linear abundance cannot be corrected by "
        "subtraction/residualisation; use no additive preprocessing for linear "
        "input until a separate multiplicative median-scaling operation is "
        'implemented."'
    )


def _scale_label(scale_kind: IntensityScaleKind | None) -> str:
    if scale_kind is None:
        return "unestablished"
    return scale_kind.value


__all__ = [
    "AdditivePreprocessingOperation",
    "AdditivePreprocessingScaleGuard",
    "SUPPORTED_ADDITIVE_PREPROCESSING_SCALE_KINDS",
]
