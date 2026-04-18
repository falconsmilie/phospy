"""Public kinase workflow shell."""

from __future__ import annotations

from phospy.api.requests import KinaseWorkflowRequest
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.references import PhosPyReferenceError
from phospy.errors.validation import WorkflowValidationError
from phospy.errors.workflows import PhosPyWorkflowError, WorkflowStageError
from phospy.workflows.kinase.contracts import (
    KinaseWorkflowExecutorContract,
    KinaseWorkflowInterpreterContract,
    KinaseWorkflowValidatorContract,
)
from phospy.workflows.kinase.executor import KinaseWorkflowExecutor
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.validator import KinaseWorkflowValidator


class KinaseWorkflow:
    """Public entrypoint for the kinase workflow."""

    def __init__(
        self,
        *,
        validator: KinaseWorkflowValidatorContract | None = None,
        interpreter: KinaseWorkflowInterpreterContract | None = None,
        executor: KinaseWorkflowExecutorContract | None = None,
    ) -> None:
        self._validator = validator or KinaseWorkflowValidator()
        self._interpreter = interpreter or KinaseWorkflowInterpreter()
        self._executor = executor or KinaseWorkflowExecutor()

    def run(self, request: KinaseWorkflowRequest) -> KinaseWorkflowResult:
        """Validate, interpret, and execute the kinase workflow."""

        try:
            validated = self._validator.run(request)
            interpreted = self._interpreter.run(validated)
            return self._executor.run(interpreted)
        except (
            PhosPyWorkflowError,
            WorkflowValidationError,
            PhosPyReferenceError,
            WorkflowStageError,
        ):
            raise
        except Exception as exc:  # pragma: no cover - defensive boundary translation
            raise PhosPyWorkflowError("kinase workflow execution failed") from exc
