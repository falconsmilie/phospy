"""Public differential-analysis entrypoint."""

from __future__ import annotations

import warnings
from typing import Protocol

from phospy.science.differential.models import DifferentialAnalysisResult


class DifferentialAnalysisRequestProtocol(Protocol):
    """Request shape accepted by the delegated differential workflow."""

    dataset: object
    design: object
    contrasts: tuple[object, ...]
    config: object


class DifferentialAnalysisWorkflowContract(Protocol):
    """Workflow collaborator required by the deprecated science shell."""

    def run(
        self,
        request: DifferentialAnalysisRequestProtocol,
    ) -> DifferentialAnalysisResult:
        """Execute differential analysis for a validated request."""
        ...


class DifferentialAnalysis:
    """Deprecated legacy shell for differential workflow execution."""

    def __init__(
        self,
        *,
        workflow: DifferentialAnalysisWorkflowContract,
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
        self._workflow = workflow

    def run(
        self,
        request: DifferentialAnalysisRequestProtocol,
    ) -> DifferentialAnalysisResult:
        return self._workflow.run(request)
