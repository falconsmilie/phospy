"""Transformation-domain exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyTransformationError(PhosPyError):
    """Intensity scale state or transformer failure."""


class InvalidTransformationStateError(PhosPyTransformationError):
    """Intensity scale state violates domain invariants."""


class TransformerExecutionError(PhosPyTransformationError):
    """A transformer failed while establishing intensity scale state."""


class TransformationStateEstablishmentError(PhosPyTransformationError):
    """Intensity scale state could not be established with supported evidence."""
