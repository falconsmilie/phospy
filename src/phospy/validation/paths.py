from __future__ import annotations

from pathlib import Path

from .errors import RequestValidationError


def validate_existing_file_path(
    path: str | Path,
    *,
    context: str,
) -> Path:
    """Validate that ``path`` exists and refers to a readable file boundary.

    Public file-based entry points should use this helper so missing or invalid
    file paths fail with a consistent package-level validation error.
    """

    resolved = Path(path)
    if not resolved.exists():
        msg = f"Invalid {context}: Path does not exist: {resolved}"
        raise RequestValidationError(msg)
    if not resolved.is_file():
        msg = f"Invalid {context}: Path is not a file: {resolved}"
        raise RequestValidationError(msg)
    return resolved
