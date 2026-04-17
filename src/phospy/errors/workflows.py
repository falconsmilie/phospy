"""Workflow execution exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyWorkflowError(PhosPyError):
    """Workflow-stage orchestration or execution failure."""


class WorkflowStageError(PhosPyWorkflowError):
    """An internal workflow stage failed after request interpretation."""
