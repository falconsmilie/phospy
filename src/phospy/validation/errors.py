from __future__ import annotations

from collections.abc import Sequence

from pydantic import ValidationError


class PhospyValidationError(ValueError):
    """Base class for phospy validation failures."""


class RequestValidationError(PhospyValidationError):
    """Raised when a validated request object cannot be created."""

    @classmethod
    def from_pydantic(
        cls,
        *,
        context: str,
        error: ValidationError,
    ) -> RequestValidationError:
        details = _format_pydantic_errors(error)
        return cls(f"{context}: {details}")


class TableSchemaError(PhospyValidationError):
    """Raised when a tabular input fails schema validation."""


class InputCompatibilityError(PhospyValidationError):
    """Raised when otherwise valid inputs are incompatible with each other."""


def _format_pydantic_errors(error: ValidationError) -> str:
    parts: list[str] = []
    for item in error.errors(include_url=False):
        location = _format_error_location(item.get("loc", ()))
        message = str(item.get("msg", "Invalid value"))
        parts.append(f"{location}: {message}" if location else message)
    return "; ".join(parts) if parts else str(error)


def _format_error_location(location: Sequence[object]) -> str:
    return ".".join(str(part) for part in location if part not in {None, ""})
