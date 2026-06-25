"""Stage-order validation wrapper for batch-correction workflows."""

from __future__ import annotations

from phospy.contracts.configs.preprocessing import (
    INTERNAL_BATCH_CORRECTION_STAGE_ORDERS,
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
    SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER,
    InternalBatchCorrectionRequest,
    InternalBatchCorrectionStageOrder,
)
from phospy.errors.input import PhosPyInputError

SUPPORTED_BATCH_CORRECTION_STAGE_ORDER_POLICY = (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER
)
SUPPORTED_BATCH_CORRECTION_EXECUTED_STAGE_ORDER = (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER
)


class BatchCorrectionWorkflowStageOrderValidator:
    """Validate internal batch-correction stage-order policy selection."""

    def run(self, *, config: InternalBatchCorrectionRequest) -> None:
        resolve_supported_batch_correction_stage_order(config.stage_order)


def resolve_supported_batch_correction_stage_order(
    stage_order: object,
) -> tuple[str, ...]:
    """Return the executable dataset correction stage order or reject it."""

    try:
        resolved_stage_order = InternalBatchCorrectionStageOrder.parse(
            stage_order,
            field_name="batch-correction workflow stage_order",
        )
    except PhosPyInputError as exc:
        raise PhosPyInputError(
            "batch-correction workflow stage_order is not a supported internal "
            "batch-correction stage-order policy"
        ) from exc

    if resolved_stage_order not in INTERNAL_BATCH_CORRECTION_STAGE_ORDERS:
        raise PhosPyInputError(
            "batch-correction workflow stage_order is not a supported internal "
            "batch-correction stage-order policy"
        )
    if resolved_stage_order is SUPPORTED_BATCH_CORRECTION_STAGE_ORDER_POLICY:
        return SUPPORTED_BATCH_CORRECTION_EXECUTED_STAGE_ORDER

    raise PhosPyInputError(
        "batch-correction workflow stage_order="
        f"{resolved_stage_order.value!r} is unsupported by the current dataset "
        "preprocessing pipeline; requested stage order implies "
        f"{' -> '.join(_implied_stage_order(resolved_stage_order))}; supported "
        "stage order is "
        f"{' -> '.join(SUPPORTED_BATCH_CORRECTION_EXECUTED_STAGE_ORDER)}; "
        "provenance must match the actual executed pipeline"
    )


def _implied_stage_order(
    stage_order: InternalBatchCorrectionStageOrder,
) -> tuple[str, ...]:
    if stage_order is (
        InternalBatchCorrectionStageOrder.AFTER_INTENSITY_TRANSFORM_BEFORE_MISSING_DATA
    ):
        return ("intensity_transform", "batch_correction", "missing_data")
    if stage_order is (
        InternalBatchCorrectionStageOrder.AFTER_TOTAL_PROTEIN_CORRECTION_BEFORE_DOWNSTREAM
    ):
        return (
            "total_protein_correction",
            "batch_correction",
            "downstream_workflows",
        )
    return SUPPORTED_BATCH_CORRECTION_EXECUTED_STAGE_ORDER


__all__ = [
    "BatchCorrectionWorkflowStageOrderValidator",
    "SUPPORTED_BATCH_CORRECTION_EXECUTED_STAGE_ORDER",
    "SUPPORTED_BATCH_CORRECTION_STAGE_ORDER_POLICY",
    "resolve_supported_batch_correction_stage_order",
]
