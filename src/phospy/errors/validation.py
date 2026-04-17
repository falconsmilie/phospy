"""Validation exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyValidationError(PhosPyError):
    """Validation failure raised for invalid public or internal inputs."""
