"""Internal executor for differential workflow requests."""

from __future__ import annotations

from phospy.api.requests import MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG
from phospy.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.differential.models import DifferentialAnalysisResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
)


class DifferentialAnalysisExecutor:
    """Run differential computation for interpreted execution inputs."""

    def __init__(
        self,
        *,
        computation_executor: DifferentialComputationExecutor | None = None,
    ) -> None:
        self._computation_executor = (
            computation_executor or DifferentialComputationExecutor()
        )

    def run(
        self, request: InterpretedDifferentialAnalysisRequest
    ) -> DifferentialAnalysisResult:
        if not isinstance(request, InterpretedDifferentialAnalysisRequest):
            raise WorkflowBoundaryError(
                seam="differential.executor.interpreted_request_type",
                next_action=(
                    "pass interpreter output into DifferentialAnalysisExecutor.run"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if (
            request.multiple_testing.method
            != MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG
        ):
            raise WorkflowBoundaryError(
                seam="differential.executor.multiple_testing_method",
                next_action=(
                    "use multiple_testing.method='benjamini_hochberg' for this release"
                ),
                details={"method": request.multiple_testing.method},
                message_prefix="differential workflow boundary validation failed",
            )
        return self._computation_executor.run(request.computation_request)


__all__ = ["DifferentialAnalysisExecutor"]
