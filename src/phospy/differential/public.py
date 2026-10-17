"""Public differential workflow shell."""

from __future__ import annotations

from phospy.differential.contracts import DifferentialAnalysisRequest
from phospy.differential.execution.workflow import DifferentialAnalysisExecutor
from phospy.differential.planning.interpreter import DifferentialAnalysisInterpreter
from phospy.differential.planning.models import (
    DifferentialAnalysisExecutorContract,
    DifferentialAnalysisInterpreterContract,
    DifferentialAnalysisValidatorContract,
)
from phospy.differential.results import DifferentialAnalysisResult
from phospy.differential.validation.workflow import DifferentialAnalysisValidator


class DifferentialAnalysisWorkflow:
    """Public entrypoint for differential analysis workflow execution.

    By default, this workflow rejects datasets containing cells imputed during
    upstream dataset preparation. This boundary is controlled by
    ``DifferentialAnalysisConfig.imputed_value_policy="reject"``. Callers who
    deliberately want to analyse an imputed dataset must explicitly select
    ``"withhold_imputed_features"`` or provide a non-imputed dataset.

    The non-default policy performs feature-level withholding, not
    observed-only fitting. A retained tested feature can still contain imputed
    cells within the configured threshold, and those cells participate in
    fitting.
    """

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
