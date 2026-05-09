"""Public differential-analysis shell."""

from __future__ import annotations

from phospy.differential.executor import DifferentialAnalysisExecutor
from phospy.differential.models import (
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
)
from phospy.differential.validator import DifferentialAnalysisRequestValidator


class DifferentialAnalysis:
    """Public entrypoint for limma-style moderated differential analysis."""

    def __init__(
        self,
        *,
        validator: DifferentialAnalysisRequestValidator | None = None,
        executor: DifferentialAnalysisExecutor | None = None,
    ) -> None:
        self._validator = validator or DifferentialAnalysisRequestValidator()
        self._executor = executor or DifferentialAnalysisExecutor()

    def run(self, request: DifferentialAnalysisRequest) -> DifferentialAnalysisResult:
        validated = self._validator.run(request)
        return self._executor.run(validated)
