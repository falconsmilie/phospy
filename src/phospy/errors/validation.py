"""Validation exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyValidationError(PhosPyError):
    """Validation failure raised for invalid public or internal inputs."""


class DatasetValidationError(PhosPyValidationError):
    """Dataset contract validation failed."""


class ReferenceValidationError(PhosPyValidationError):
    """Reference contract validation failed."""


class TransformationValidationError(PhosPyValidationError):
    """Transformation-state validation failed."""


class WorkflowValidationError(PhosPyValidationError):
    """Workflow-level validation failed."""
