"""Public kinase workflow shell."""

from __future__ import annotations

from phospy.kinase.contracts import (
    KinaseWorkflowExecutorContract,
    KinaseWorkflowInterpreterContract,
    KinaseWorkflowValidatorContract,
)
from phospy.kinase.executor import KinaseWorkflowExecutor
from phospy.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.kinase.requests import KinaseWorkflowRequest
from phospy.kinase.results import KinaseWorkflowResult
from phospy.kinase.validator import KinaseWorkflowValidator


class KinaseWorkflow:
    """Public entrypoint for the kinase workflow."""

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
        validator: KinaseWorkflowValidatorContract | None = None,
        interpreter: KinaseWorkflowInterpreterContract | None = None,
        executor: KinaseWorkflowExecutorContract | None = None,
    ) -> KinaseWorkflow:
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
        validator: KinaseWorkflowValidatorContract | None,
        interpreter: KinaseWorkflowInterpreterContract | None,
        executor: KinaseWorkflowExecutorContract | None,
    ) -> None:
        self._validator = validator or KinaseWorkflowValidator()
        self._interpreter = interpreter or KinaseWorkflowInterpreter()
        self._executor = executor or KinaseWorkflowExecutor()

    def run(self, request: KinaseWorkflowRequest) -> KinaseWorkflowResult:
        """Validate, interpret, and execute the kinase workflow."""
        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return self._executor.run(interpreted)
