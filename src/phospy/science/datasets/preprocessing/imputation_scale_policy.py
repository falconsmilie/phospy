"""Scientific policy for imputation input scale and preprocessing order."""

from __future__ import annotations

from dataclasses import dataclass

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.plan_constants import (
    DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM,
    DATASET_PREPROCESSING_STAGE_MISSING_DATA,
)
from phospy.science.datasets.preprocessing.policy_models import (
    ImputationInputScale,
    IntensityTransformPolicy,
    MissingDataPolicy,
)
from phospy.science.transformations.models import IntensityScaleKind

IMPUTATION_INPUT_SCALE_SOURCE_CALLER_SELECTED = "caller_selected"
IMPUTATION_INPUT_SCALE_SOURCE_METHOD_REQUIRED = "method_required"

IMPUTATION_OPERATION_ORDER_BEFORE_INTENSITY_TRANSFORM = "before_intensity_transform"
IMPUTATION_OPERATION_ORDER_AFTER_INTENSITY_TRANSFORM = "after_intensity_transform"
IMPUTATION_OPERATION_ORDER_NO_INTENSITY_TRANSFORM = "no_intensity_transform"
IMPUTATION_OPERATION_ORDERS = frozenset(
    {
        IMPUTATION_OPERATION_ORDER_BEFORE_INTENSITY_TRANSFORM,
        IMPUTATION_OPERATION_ORDER_AFTER_INTENSITY_TRANSFORM,
        IMPUTATION_OPERATION_ORDER_NO_INTENSITY_TRANSFORM,
    }
)


@dataclass(frozen=True, slots=True)
class ResolvedImputationScalePolicy:
    """Resolved imputation input domain and operation-order policy."""

    input_scale: ImputationInputScale | None
    input_scale_source: str | None
    operation_order: str | None


def resolve_imputation_scale_policy(
    *,
    missing_data_policy: MissingDataPolicy,
    requested_input_scale: object | None,
    intensity_transform_policy: IntensityTransformPolicy,
) -> ResolvedImputationScalePolicy:
    """Resolve the scale policy owned by missing-data preprocessing."""

    policy = MissingDataPolicy.parse(
        missing_data_policy,
        field_name="dataset preprocessing plan missing_data_policy",
    )
    transform_policy = IntensityTransformPolicy.parse(
        intensity_transform_policy,
        field_name="dataset preprocessing plan intensity_transform_policy",
    )
    input_scale = _resolve_requested_input_scale(requested_input_scale)
    if policy is MissingDataPolicy.FORBID:
        _reject_forbid_input_scale(input_scale)
        return ResolvedImputationScalePolicy(
            input_scale=None,
            input_scale_source=None,
            operation_order=None,
        )
    if policy in {MissingDataPolicy.IMPUTE_ROW_MEDIAN, MissingDataPolicy.IMPUTE_KNN}:
        input_scale = _require_caller_selected_input_scale(
            policy=policy,
            input_scale=input_scale,
        )
        return ResolvedImputationScalePolicy(
            input_scale=input_scale,
            input_scale_source=IMPUTATION_INPUT_SCALE_SOURCE_CALLER_SELECTED,
            operation_order=_resolve_operation_order(
                input_scale=input_scale,
                intensity_transform_policy=transform_policy,
            ),
        )
    if policy is MissingDataPolicy.IMPUTE_MINPROB:
        input_scale = _resolve_minprob_input_scale(input_scale)
        return ResolvedImputationScalePolicy(
            input_scale=input_scale,
            input_scale_source=IMPUTATION_INPUT_SCALE_SOURCE_METHOD_REQUIRED,
            operation_order=_resolve_operation_order(
                input_scale=input_scale,
                intensity_transform_policy=transform_policy,
            ),
        )
    raise PhosPyInputError(
        "dataset preprocessing plan missing_data_policy contains an unsupported "
        "imputation policy"
    )


def reject_unestablished_log2_imputation_input_scale(
    *,
    resolved_policy: ResolvedImputationScalePolicy,
    intensity_transform_policy: IntensityTransformPolicy,
    declared_input_scale_kind: IntensityScaleKind | None,
) -> None:
    """Reject log2 imputation without an upstream transform or input declaration."""

    if resolved_policy.input_scale is not ImputationInputScale.LOG2:
        return
    transform_policy = IntensityTransformPolicy.parse(
        intensity_transform_policy,
        field_name="dataset preprocessing plan intensity_transform_policy",
    )
    if transform_policy is IntensityTransformPolicy.LOG2:
        return
    declared_scale = _resolve_declared_input_scale_kind(declared_input_scale_kind)
    if declared_scale is IntensityScaleKind.LOG2:
        return
    raise PhosPyInputError(
        "dataset build request preprocessing_config.missing_data.input_scale='log2' "
        "requires either preprocessing_config.intensity_transform.policy='log2' "
        "or dataset build request input_intensity_scale='log2'. This prevents an "
        "imputer from silently changing scale without a recorded transition."
    )


def reject_incompatible_imputation_stage_order(
    *,
    stage_order: tuple[str, ...],
    resolved_policy: ResolvedImputationScalePolicy,
) -> None:
    """Reject manual stage orders that contradict resolved imputation scale."""

    if resolved_policy.input_scale is None:
        return
    stages = tuple(str(stage) for stage in stage_order)
    if (
        DATASET_PREPROCESSING_STAGE_MISSING_DATA not in stages
        or DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM not in stages
    ):
        return
    missing_data_index = stages.index(DATASET_PREPROCESSING_STAGE_MISSING_DATA)
    transform_index = stages.index(DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM)
    if (
        resolved_policy.input_scale is ImputationInputScale.LINEAR
        and transform_index < missing_data_index
    ):
        raise PhosPyInputError(
            "dataset preprocessing plan stage_order is incompatible with "
            "missing_data.input_scale='linear': missing_data must run before "
            "intensity_transform so imputation consumes linear values."
        )
    if (
        resolved_policy.input_scale is ImputationInputScale.LOG2
        and missing_data_index < transform_index
    ):
        raise PhosPyInputError(
            "dataset preprocessing plan stage_order is incompatible with "
            "missing_data.input_scale='log2': missing_data must run after "
            "intensity_transform or on input declared as already log2."
        )


def imputation_input_scale_kind(
    input_scale: ImputationInputScale | str | None,
) -> IntensityScaleKind:
    """Return the quantitative scale kind for an imputation input-scale policy."""

    scale = _require_input_scale_for_execution(input_scale)
    if scale is ImputationInputScale.LOG2:
        return IntensityScaleKind.LOG2
    return IntensityScaleKind.LINEAR


def _resolve_requested_input_scale(value: object | None) -> ImputationInputScale | None:
    if value is None:
        return None
    return ImputationInputScale.parse(
        value,
        field_name="dataset build request preprocessing_config.missing_data.input_scale",
    )


def _resolve_declared_input_scale_kind(
    value: IntensityScaleKind | str | None,
) -> IntensityScaleKind | None:
    if value is None:
        return None
    if isinstance(value, IntensityScaleKind):
        return value
    try:
        return IntensityScaleKind(str(value))
    except ValueError as exc:
        supported = ", ".join(member.value for member in IntensityScaleKind)
        raise PhosPyInputError(
            f"dataset build request input_intensity_scale must be one of: {supported}"
        ) from exc


def _reject_forbid_input_scale(input_scale: ImputationInputScale | None) -> None:
    if input_scale is None:
        return
    raise PhosPyInputError(
        "dataset build request preprocessing_config.missing_data.input_scale must "
        "be None when missing_data.policy='forbid'"
    )


def _require_caller_selected_input_scale(
    *,
    policy: MissingDataPolicy,
    input_scale: ImputationInputScale | None,
) -> ImputationInputScale:
    if input_scale is not None:
        return input_scale
    raise PhosPyInputError(
        "ambiguous missing-data imputation configuration: "
        "preprocessing_config.missing_data.input_scale is required when "
        f"missing_data.policy='{policy.value}'. Set input_scale='linear' to "
        "impute before log2 transformation, or input_scale='log2' to impute after "
        "a recorded log2 transform / declared log2 input."
    )


def _resolve_minprob_input_scale(
    input_scale: ImputationInputScale | None,
) -> ImputationInputScale:
    if input_scale in {None, ImputationInputScale.LOG2}:
        return ImputationInputScale.LOG2
    raise PhosPyInputError(
        "invalid missing-data imputation scale configuration: "
        "missing_data.policy='impute_minprob' requires "
        "preprocessing_config.missing_data.input_scale='log2' because MinProb's "
        "left-censored sampling model operates on log2 intensities."
    )


def _resolve_operation_order(
    *,
    input_scale: ImputationInputScale,
    intensity_transform_policy: IntensityTransformPolicy,
) -> str:
    if intensity_transform_policy is IntensityTransformPolicy.IDENTITY:
        return IMPUTATION_OPERATION_ORDER_NO_INTENSITY_TRANSFORM
    if input_scale is ImputationInputScale.LOG2:
        return IMPUTATION_OPERATION_ORDER_AFTER_INTENSITY_TRANSFORM
    return IMPUTATION_OPERATION_ORDER_BEFORE_INTENSITY_TRANSFORM


def _require_input_scale_for_execution(
    input_scale: ImputationInputScale | str | None,
) -> ImputationInputScale:
    if input_scale is None:
        raise PhosPyInputError(
            "dataset preprocessing plan missing_data_input_scale is required "
            "for imputation execution"
        )
    return ImputationInputScale.parse(
        input_scale,
        field_name="dataset preprocessing plan missing_data_input_scale",
    )


__all__ = [
    "IMPUTATION_INPUT_SCALE_SOURCE_CALLER_SELECTED",
    "IMPUTATION_INPUT_SCALE_SOURCE_METHOD_REQUIRED",
    "IMPUTATION_OPERATION_ORDER_AFTER_INTENSITY_TRANSFORM",
    "IMPUTATION_OPERATION_ORDER_BEFORE_INTENSITY_TRANSFORM",
    "IMPUTATION_OPERATION_ORDER_NO_INTENSITY_TRANSFORM",
    "IMPUTATION_OPERATION_ORDERS",
    "ResolvedImputationScalePolicy",
    "imputation_input_scale_kind",
    "reject_incompatible_imputation_stage_order",
    "reject_unestablished_log2_imputation_input_scale",
    "resolve_imputation_scale_policy",
]
