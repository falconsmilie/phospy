"""Public differential workflow shell."""

from __future__ import annotations

from phospy.contracts.requests import DifferentialAnalysisRequest
from phospy.science.differential.models import DifferentialAnalysisResult
from phospy.workflows.differential.executor import DifferentialAnalysisExecutor
from phospy.workflows.differential.interpreter import DifferentialAnalysisInterpreter
from phospy.workflows.differential.models import (
    DifferentialAnalysisExecutorContract,
    DifferentialAnalysisInterpreterContract,
    DifferentialAnalysisValidatorContract,
)
from phospy.workflows.differential.validator import DifferentialAnalysisValidator


class DifferentialAnalysisWorkflow:
    """Public entrypoint for differential analysis workflow execution."""

    def __init__(self) -> None:
        self._init_components(
            validator=None,
            interpreter=None,
            executor=None,
        )

    @classmethod
    def _with_components(
        cls,
        *,
        validator: DifferentialAnalysisValidatorContract | None = None,
        interpreter: DifferentialAnalysisInterpreterContract | None = None,
        executor: DifferentialAnalysisExecutorContract | None = None,
    ) -> DifferentialAnalysisWorkflow:
        workflow = cls.__new__(cls)
        workflow._init_components(
            validator=validator,
            interpreter=interpreter,
            executor=executor,
        )
        return workflow

    def _init_components(
        self,
        *,
        validator: DifferentialAnalysisValidatorContract | None,
        interpreter: DifferentialAnalysisInterpreterContract | None,
        executor: DifferentialAnalysisExecutorContract | None,
    ) -> None:
        self._validator = validator or DifferentialAnalysisValidator()
        self._interpreter = interpreter or DifferentialAnalysisInterpreter()
        self._executor = executor or DifferentialAnalysisExecutor()

    def run(self, request: DifferentialAnalysisRequest) -> DifferentialAnalysisResult:
        """Validate, interpret, and execute differential analysis."""
        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return self._executor.run(interpreted)


__all__ = ["DifferentialAnalysisWorkflow"]
