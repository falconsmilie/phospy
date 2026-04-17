"""Input-boundary exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyInputError(PhosPyError):
    """Input read/interpretation failure."""


class UnsupportedInputFormatError(PhosPyInputError):
    """Input shape or format is not supported by the current ingestion contract."""
