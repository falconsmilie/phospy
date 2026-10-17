"""Validation helpers for FragPipe/PTMProphet importer options."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError

FRAGPIPE_FLAG_POLICY_REMOVE = "remove"
FRAGPIPE_FLAG_POLICY_FLAG = "flag"
FRAGPIPE_FLAG_POLICY_ERROR = "error"
SUPPORTED_FRAGPIPE_FLAG_POLICIES: tuple[str, ...] = (
    FRAGPIPE_FLAG_POLICY_REMOVE,
    FRAGPIPE_FLAG_POLICY_FLAG,
    FRAGPIPE_FLAG_POLICY_ERROR,
)

FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE = "peptide"
FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PROTEIN = "protein"
SUPPORTED_FRAGPIPE_PTMPROPHET_POSITION_REFERENCES: tuple[str, ...] = (
    FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE,
    FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PROTEIN,
)


def validate_fragpipe_flag_policy(value: object, *, field_name: str) -> str:
    if isinstance(value, str) and value.strip() in SUPPORTED_FRAGPIPE_FLAG_POLICIES:
        return value.strip()
    supported = ", ".join(repr(item) for item in SUPPORTED_FRAGPIPE_FLAG_POLICIES)
    raise PhosPyInputError(f"{field_name} must be one of: {supported}")


def validate_optional_fragpipe_column_name(
    value: object,
    *,
    field_name: str,
) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string or None")
    return value.strip()


def validate_ptmprophet_position_reference(
    value: object,
    *,
    field_name: str,
) -> str:
    if (
        isinstance(value, str)
        and value.strip() in SUPPORTED_FRAGPIPE_PTMPROPHET_POSITION_REFERENCES
    ):
        return value.strip()
    supported = ", ".join(
        repr(item) for item in SUPPORTED_FRAGPIPE_PTMPROPHET_POSITION_REFERENCES
    )
    raise PhosPyInputError(f"{field_name} must be one of: {supported}")


__all__ = [
    "FRAGPIPE_FLAG_POLICY_ERROR",
    "FRAGPIPE_FLAG_POLICY_FLAG",
    "FRAGPIPE_FLAG_POLICY_REMOVE",
    "FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PEPTIDE",
    "FRAGPIPE_PTMPROPHET_POSITION_REFERENCE_PROTEIN",
    "SUPPORTED_FRAGPIPE_FLAG_POLICIES",
    "SUPPORTED_FRAGPIPE_PTMPROPHET_POSITION_REFERENCES",
    "validate_fragpipe_flag_policy",
    "validate_optional_fragpipe_column_name",
    "validate_ptmprophet_position_reference",
]
