"""Missingness validation wrapper for batch-correction workflows."""

from __future__ import annotations

from phospy.contracts.configs.preprocessing import (
    CorrectionMissingnessCompatibilityValidator,
    CorrectionMissingnessPolicy,
)
from phospy.workflows.batch_correction.contracts import BatchCorrectionWorkflowRequest


class BatchCorrectionWorkflowMissingnessValidator:
    """Validate correction missingness policy and provide a complete-data default."""

    def __init__(
        self,
        *,
        compatibility_validator: CorrectionMissingnessCompatibilityValidator
        | None = None,
    ) -> None:
        self._compatibility_validator = (
            compatibility_validator or CorrectionMissingnessCompatibilityValidator()
        )

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
    ) -> CorrectionMissingnessPolicy:
        policy = request.missingness_policy
        self._compatibility_validator.run(
            phospho=request.phospho,
            policy=policy,
            context="batch-correction workflow missingness validation",
        )
        if policy is None:
            return CorrectionMissingnessPolicy()
        return policy


__all__ = ["BatchCorrectionWorkflowMissingnessValidator"]
