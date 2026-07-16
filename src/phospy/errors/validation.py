"""Validation exceptions."""

from __future__ import annotations

from phospy.errors.base import PhosPyError


class PhosPyValidationError(PhosPyError):
    """Validation failure raised for invalid public or internal inputs."""


class DatasetValidationError(PhosPyValidationError):
    """Dataset contract validation failed."""


class ReferenceValidationError(PhosPyValidationError):
    """Reference contract validation failed."""


class ReferenceIdentifierNormalisationValidationError(ReferenceValidationError):
    """Reference validation error carrying identifier-normalisation provenance."""

    message: str
    identifier_normalisation_report: object

    def __init__(
        self,
        message: str,
        identifier_normalisation_report: object,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.identifier_normalisation_report = identifier_normalisation_report


class TransformationValidationError(PhosPyValidationError):
    """Intensity-scale-state validation failed."""


class WorkflowValidationError(PhosPyValidationError):
    """Workflow-level validation failed."""
