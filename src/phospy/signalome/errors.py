"""Signalome workflow exceptions."""

from phospy.errors.validation import WorkflowValidationError


class SignalomeScaleError(WorkflowValidationError):
    """Signalome scale guard blocked an unsafe expensive execution path."""


class SignalomeModuleCountValidationError(WorkflowValidationError):
    """Signalome module-count request violates clustering-site constraints."""


__all__ = ["SignalomeModuleCountValidationError", "SignalomeScaleError"]
