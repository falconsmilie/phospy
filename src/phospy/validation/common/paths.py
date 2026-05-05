"""Shared filesystem-path validation helpers."""

from __future__ import annotations

from phospy.validation.common.config_values import (
    ValidationErrorType,
    require_non_empty_string,
)


def require_local_filesystem_path(
    value: object,
    *,
    field_name: str,
    error_type: ValidationErrorType,
    when_provided: bool = False,
) -> str:
    """Require one local filesystem path string and reject remote URLs."""

    path = require_non_empty_string(
        value,
        field_name=field_name,
        error_type=error_type,
        when_provided=when_provided,
    )
    if "://" in path.lower():
        raise error_type(
            f"{field_name} must be a local filesystem path; remote URLs are not "
            "supported"
        )
    return path
