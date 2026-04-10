from __future__ import annotations

from pathlib import Path

from ..errors import RequestValidationError


def validate_existing_file_path(
    path: str | Path,
    *,
    context: str,
) -> Path:
    """Validate that ``path`` exists and refers to a file."""

    resolved = Path(path)
    if not resolved.exists():
        msg = f"Invalid {context}: Path does not exist: {resolved}"
        raise RequestValidationError(msg)
    if not resolved.is_file():
        msg = f"Invalid {context}: Path is not a file: {resolved}"
        raise RequestValidationError(msg)
    return resolved


__all__ = ["validate_existing_file_path"]
