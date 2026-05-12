"""Public differential-analysis entrypoint."""

from __future__ import annotations

import warnings

from phospy.api.requests import DifferentialAnalysisRequest
from phospy.differential.models import DifferentialAnalysisResult
from phospy.workflows.differential.public import DifferentialAnalysisWorkflow


class DifferentialAnalysis:
    """Deprecated legacy shell for differential workflow execution."""

    def __init__(
        self,
        *,
        workflow: DifferentialAnalysisWorkflow | None = None,
    ) -> None:
        warnings.warn(
            (
                "DifferentialAnalysis is deprecated; use "
                "DifferentialAnalysisWorkflow from top-level phospy "
                "(from phospy import DifferentialAnalysisWorkflow)."
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        self._workflow = workflow or DifferentialAnalysisWorkflow()

    def run(self, request: DifferentialAnalysisRequest) -> DifferentialAnalysisResult:
        return self._workflow.run(request)
