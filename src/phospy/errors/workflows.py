"""Workflow execution exceptions."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.errors.base import PhosPyError
from phospy.errors.validation import WorkflowValidationError


class PhosPyWorkflowError(PhosPyError):
    """Workflow-stage orchestration or execution failure."""


class WorkflowBoundaryError(WorkflowValidationError):
    """Workflow boundary validation failed with actionable diagnostics."""

    def __init__(
        self,
        message: str | None = None,
        *,
        seam: str | None = None,
        next_action: str | None = None,
        details: Mapping[str, object] | None = None,
        message_prefix: str = "workflow boundary validation failed",
    ) -> None:
        self.seam = seam
        self.next_action = next_action
        self.details = dict(details) if details is not None else {}
        resolved_message = self._resolve_message(
            message=message,
            message_prefix=message_prefix,
        )
        super().__init__(resolved_message)

    def _resolve_message(
        self,
        *,
        message: str | None,
        message_prefix: str,
    ) -> str:
        if message is not None:
            return message
        if self.seam is None and self.next_action is None:
            return message_prefix
        seam_text = self.seam if self.seam is not None else "unknown"
        next_action_text = (
            self.next_action
            if self.next_action is not None
            else "no next_action provided"
        )
        details_text = ", ".join(
            f"{key}={value}" for key, value in self.details.items()
        )
        if details_text == "":
            return (
                f"{message_prefix} at seam={seam_text}; next_action={next_action_text}"
            )
        return (
            f"{message_prefix} at seam={seam_text}; "
            f"{details_text}; next_action={next_action_text}"
        )


class WorkflowStageError(PhosPyWorkflowError):
    """An internal workflow stage failed after request interpretation."""


class SignalomeScaleError(WorkflowValidationError):
    """Signalome scale guard blocked an unsafe expensive execution path."""


class SignalomeModuleCountValidationError(WorkflowValidationError):
    """Signalome module-count request violates clustering-site constraints."""
