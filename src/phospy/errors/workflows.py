"""Workflow execution exceptions."""

from phospy.errors.base import PhosPyError
from phospy.errors.validation import WorkflowValidationError


class PhosPyWorkflowError(PhosPyError):
    """Workflow-stage orchestration or execution failure."""


class WorkflowBoundaryError(WorkflowValidationError):
    """Workflow boundary validation failed with actionable diagnostics."""


class WorkflowStageError(PhosPyWorkflowError):
    """An internal workflow stage failed after request interpretation."""
