"""Reference-domain exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyReferenceError(PhosPyError):
    """Reference resolution or compatibility failure."""


class ReferenceResolutionError(PhosPyReferenceError):
    """Failed to resolve a concrete reference bundle."""


class ReferenceCompatibilityError(PhosPyReferenceError):
    """Dataset and reference organisms are incompatible."""


class UnsupportedOrganismError(PhosPyReferenceError):
    """The requested organism is not supported by the selected reference source."""
