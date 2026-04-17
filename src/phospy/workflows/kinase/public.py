"""Public kinase workflow shell."""

from __future__ import annotations

from phospy.api.requests import SimpleKinaseWorkflowRequest
from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.errors.references import PhosPyReferenceError
from phospy.errors.validation import WorkflowValidationError
from phospy.errors.workflows import PhosPyWorkflowError, WorkflowStageError
from phospy.workflows.kinase.contracts import (
    KinaseWorkflowExecutorContract,
    KinaseWorkflowInterpreterContract,
    KinaseWorkflowValidatorContract,
)
from phospy.workflows.kinase.executor import SimpleKinaseWorkflowExecutor
from phospy.workflows.kinase.interpreter import SimpleKinaseWorkflowInterpreter
from phospy.workflows.kinase.validator import SimpleKinaseWorkflowValidator


class KinaseWorkflow:
    """Public entrypoint for the kinase workflow."""

    def __init__(
        self,
        *,
        validator: KinaseWorkflowValidatorContract | None = None,
        interpreter: KinaseWorkflowInterpreterContract | None = None,
        executor: KinaseWorkflowExecutorContract | None = None,
    ) -> None:
        self._validator = validator or SimpleKinaseWorkflowValidator()
        self._interpreter = interpreter or SimpleKinaseWorkflowInterpreter()
        self._executor = executor or SimpleKinaseWorkflowExecutor()

    def run(self, request: SimpleKinaseWorkflowRequest) -> SimpleKinaseWorkflowResult:
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
