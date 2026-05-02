"""Signalome workflow boundary error helpers."""

from __future__ import annotations

from typing import NoReturn

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.workflows.signalome.constants import (
    SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
)


def raise_signalome_boundary_error(
    *,
    seam: str,
    next_action: str,
    **details: object,
) -> NoReturn:
    raise WorkflowBoundaryError(
        seam=seam,
        next_action=next_action,
        details=details,
        message_prefix=SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
    )


def raise_wrapped_signalome_boundary_error(
    *,
    stage_name: str,
    seam: str,
    field_name: str,
    operation: str,
    next_action: str,
    original_error: Exception,
    **details: object,
) -> NoReturn:
    original_message = " ".join(str(original_error).split())
    message = (
        f"{stage_name} failed while {operation} for {field_name}. "
        f"Original error: {type(original_error).__name__}: {original_message}. "
        f"Next action: {next_action}"
    )
    raise WorkflowBoundaryError(
        message=message,
        seam=seam,
        next_action=next_action,
        details={
            "field_name": field_name,
            "operation": operation,
            "original_error_type": type(original_error).__name__,
            "original_error_message": original_message,
            **details,
        },
        message_prefix=SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
    ) from original_error


__all__ = [
    "raise_signalome_boundary_error",
    "raise_wrapped_signalome_boundary_error",
]
