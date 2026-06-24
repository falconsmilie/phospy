"""Stage-order validation wrapper for batch-correction workflows."""

from __future__ import annotations

from phospy.contracts.configs.preprocessing import (
    INTERNAL_BATCH_CORRECTION_STAGE_ORDERS,
    InternalBatchCorrectionRequest,
)
from phospy.errors.input import PhosPyInputError


class BatchCorrectionWorkflowStageOrderValidator:
    """Validate internal batch-correction stage-order policy selection."""

    def run(self, *, config: InternalBatchCorrectionRequest) -> None:
        if config.stage_order not in INTERNAL_BATCH_CORRECTION_STAGE_ORDERS:
            raise PhosPyInputError(
                "batch-correction workflow stage_order is not a supported internal "
                "batch-correction stage-order policy"
            )


__all__ = ["BatchCorrectionWorkflowStageOrderValidator"]
