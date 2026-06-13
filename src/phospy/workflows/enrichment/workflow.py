"""Public enrichment workflow shell."""

from __future__ import annotations

from phospy.contracts.requests import EnrichmentWorkflowRequest
from phospy.contracts.results import EnrichmentWorkflowResult
from phospy.workflows.enrichment.executor import EnrichmentWorkflowExecutor
from phospy.workflows.enrichment.interpreter import EnrichmentWorkflowInterpreter
from phospy.workflows.enrichment.models import (
    EnrichmentWorkflowExecutorContract,
    EnrichmentWorkflowInterpreterContract,
    EnrichmentWorkflowValidatorContract,
)
from phospy.workflows.enrichment.validator import EnrichmentWorkflowValidator


class EnrichmentWorkflow:
    """Public entrypoint for enrichment workflow execution."""

    def __init__(
        self,
        *,
        validator: EnrichmentWorkflowValidatorContract | None = None,
        interpreter: EnrichmentWorkflowInterpreterContract | None = None,
        executor: EnrichmentWorkflowExecutorContract | None = None,
    ) -> None:
        self._validator = validator or EnrichmentWorkflowValidator()
        self._interpreter = interpreter or EnrichmentWorkflowInterpreter()
        self._executor = executor or EnrichmentWorkflowExecutor()

    def run(self, request: EnrichmentWorkflowRequest) -> EnrichmentWorkflowResult:
        """Validate, interpret, and execute native enrichment."""
        validated = self._validator.run(request)
        interpreted = self._interpreter.run(validated)
        return self._executor.run(interpreted)


__all__ = ["EnrichmentWorkflow"]
