"""Builder-boundary exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyBuildError(PhosPyError):
    """Dataset build or preprocessing failure."""
