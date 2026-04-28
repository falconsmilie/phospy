"""Public signalome workflow shell."""

from __future__ import annotations

from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import SignalomeWorkflowResult
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

    def __init__(
        self,
        *,
        validator: SignalomeWorkflowValidatorContract | None = None,
        interpreter: SignalomeWorkflowInterpreterContract | None = None,
        executor: SignalomeWorkflowExecutorContract | None = None,
    ) -> None:
        self._validator = validator or SignalomeWorkflowValidator()
        self._interpreter = interpreter or SignalomeWorkflowInterpreter()
        self._executor = executor or SignalomeWorkflowExecutor()

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowResult:
        """Validate, interpret, and execute the signalome workflow."""
        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return self._executor.run(interpreted)
