"""Compatibility wrapper for differential workflow execution."""

from __future__ import annotations

from phospy.api.requests import DifferentialAnalysisRequest
from phospy.differential.models import DifferentialAnalysisResult
from phospy.workflows.differential.public import DifferentialAnalysisWorkflow


class DifferentialAnalysis:
    """Compatibility alias for `DifferentialAnalysisWorkflow`."""

    def __init__(
        self,
        *,
        workflow: DifferentialAnalysisWorkflow | None = None,
    ) -> None:
        self._workflow = workflow or DifferentialAnalysisWorkflow()

    def run(self, request: DifferentialAnalysisRequest) -> DifferentialAnalysisResult:
        return self._workflow.run(request)
