"""Differential model fitting collaborator."""

from __future__ import annotations

from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest as DifferentialComputationRequest,
)
from phospy.science.differential.models import DifferentialComputationResult


class DifferentialModelFitter:
    """Run the scientific differential computation on eligible inputs."""

    def __init__(
        self,
        *,
        computation_executor: DifferentialComputationExecutor | None = None,
    ) -> None:
        self._computation_executor = (
            computation_executor or DifferentialComputationExecutor()
        )

    def run(
        self,
        computation_request: DifferentialComputationRequest,
    ) -> DifferentialComputationResult:
        return self._computation_executor.run(computation_request)


__all__ = ["DifferentialModelFitter"]
