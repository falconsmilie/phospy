"""Transformation-domain exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyTransformationError(PhosPyError):
    """Transformation state or transformer failure."""


class InvalidTransformationStateError(PhosPyTransformationError):
    """Transformation state violates domain invariants."""


class TransformerExecutionError(PhosPyTransformationError):
    """A transformer failed while establishing transformation state."""
