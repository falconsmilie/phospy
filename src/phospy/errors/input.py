"""Input-boundary exceptions."""

from phospy.errors.base import PhosPyError


class PhosPyInputError(PhosPyError):
    """Input read/interpretation failure."""

    def __init__(self, message: str, *, diagnostics: object | None = None) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


class DatasetProcessingStateError(PhosPyInputError):
    """Dataset processing-state invariants are internally inconsistent."""


class UnsupportedInputFormatError(PhosPyInputError):
    """Input shape or format is not supported by the current ingestion contract."""
