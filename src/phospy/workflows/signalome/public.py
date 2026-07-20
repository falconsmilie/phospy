"""Public signalome workflow shell."""

from __future__ import annotations

from phospy.contracts.requests import SignalomeWorkflowRequest
from phospy.contracts.results import SignalomeWorkflowResult
from phospy.workflows.signalome.contracts import (
    SignalomeWorkflowExecutorContract,
    SignalomeWorkflowInterpreterContract,
    SignalomeWorkflowValidatorContract,
)
from phospy.workflows.signalome.executor import SignalomeWorkflowExecutor
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator


class SignalomeWorkflow:
    """Public entrypoint for the signalome workflow."""

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
        validator: SignalomeWorkflowValidatorContract | None = None,
        interpreter: SignalomeWorkflowInterpreterContract | None = None,
        executor: SignalomeWorkflowExecutorContract | None = None,
    ) -> SignalomeWorkflow:
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
        validator: SignalomeWorkflowValidatorContract | None,
        interpreter: SignalomeWorkflowInterpreterContract | None,
        executor: SignalomeWorkflowExecutorContract | None,
    ) -> None:
        self._validator = validator or SignalomeWorkflowValidator()
        self._interpreter = interpreter or SignalomeWorkflowInterpreter()
        self._executor = executor or SignalomeWorkflowExecutor()

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowResult:
        """Validate, interpret, and execute the signalome workflow."""
        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return self._executor.run(interpreted)
