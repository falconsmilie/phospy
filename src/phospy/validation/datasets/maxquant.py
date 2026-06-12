"""Validation helpers for MaxQuant phosphosite importer options."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError

MAXQUANT_FLAG_POLICY_REMOVE = "remove"
MAXQUANT_FLAG_POLICY_FLAG = "flag"
MAXQUANT_FLAG_POLICY_ERROR = "error"
SUPPORTED_MAXQUANT_FLAG_POLICIES: tuple[str, ...] = (
    MAXQUANT_FLAG_POLICY_REMOVE,
    MAXQUANT_FLAG_POLICY_FLAG,
    MAXQUANT_FLAG_POLICY_ERROR,
)


def validate_maxquant_flag_policy(value: object, *, field_name: str) -> str:
    if isinstance(value, str) and value.strip() in SUPPORTED_MAXQUANT_FLAG_POLICIES:
        return value.strip()
    supported = ", ".join(repr(item) for item in SUPPORTED_MAXQUANT_FLAG_POLICIES)
    raise PhosPyInputError(f"{field_name} must be one of: {supported}")


def validate_optional_maxquant_column_name(
    value: object,
    *,
    field_name: str,
) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string or None")
    return value.strip()


__all__ = [
    "MAXQUANT_FLAG_POLICY_ERROR",
    "MAXQUANT_FLAG_POLICY_FLAG",
    "MAXQUANT_FLAG_POLICY_REMOVE",
    "SUPPORTED_MAXQUANT_FLAG_POLICIES",
    "validate_maxquant_flag_policy",
    "validate_optional_maxquant_column_name",
]
