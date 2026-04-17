"""Reference-domain exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyReferenceError(PhosPyError):
    """Reference resolution or compatibility failure."""
