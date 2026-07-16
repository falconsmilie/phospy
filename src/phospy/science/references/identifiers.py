"""Reference identifier normalisation boundary and provenance records."""

from __future__ import annotations

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.reference_identifiers import (
    REFERENCE_IDENTIFIER_NORMALISATION_SCHEMA_VERSION,
    ReferenceIdentifierNormalisationRecord,
    ReferenceIdentifierNormalisationReport,
    build_reference_identifier_normalisation_report,
    merge_reference_identifier_normalisation_reports,
)
from phospy.science.sites.identifiers import canonicalize_site_identifier


def normalise_reference_site_id(
    value: object,
    *,
    table_name: str,
    column_name: str,
    row_position: int,
) -> ReferenceIdentifierNormalisationRecord:
    original_value = _stringify_identifier(value)
    field_name = f"{table_name}.{column_name}[{row_position}]"
    try:
        normalised = canonicalize_site_identifier(
            value,
            field_name=field_name,
            error_type=ReferenceValidationError,
        )
    except ReferenceValidationError as exc:
        return ReferenceIdentifierNormalisationRecord(
            table_name=table_name,
            column_name=column_name,
            row_position=row_position,
            identifier_kind="site_id",
            original_value=original_value,
            normalised_value=None,
            status="invalid",
            reason=str(exc),
        )
    status = "unchanged" if normalised == original_value else "normalised"
    return ReferenceIdentifierNormalisationRecord(
        table_name=table_name,
        column_name=column_name,
        row_position=row_position,
        identifier_kind="site_id",
        original_value=original_value,
        normalised_value=normalised,
        status=status,
        reason=None,
    )


def normalise_reference_kinase_id(
    value: object,
    *,
    table_name: str,
    column_name: str,
    row_position: int,
) -> ReferenceIdentifierNormalisationRecord:
    original_value = _stringify_identifier(value)
    if _is_missing(value):
        return ReferenceIdentifierNormalisationRecord(
            table_name=table_name,
            column_name=column_name,
            row_position=row_position,
            identifier_kind="kinase",
            original_value=original_value,
            normalised_value=None,
            status="invalid",
            reason=(
                f"{table_name}.{column_name}[{row_position}] "
                "must not contain missing values"
            ),
        )
    if not isinstance(value, str):
        return ReferenceIdentifierNormalisationRecord(
            table_name=table_name,
            column_name=column_name,
            row_position=row_position,
            identifier_kind="kinase",
            original_value=original_value,
            normalised_value=None,
            status="invalid",
            reason=(
                f"{table_name}.{column_name}[{row_position}] "
                "must contain non-empty string values"
            ),
        )
    stripped = value.strip()
    if stripped == "":
        return ReferenceIdentifierNormalisationRecord(
            table_name=table_name,
            column_name=column_name,
            row_position=row_position,
            identifier_kind="kinase",
            original_value=original_value,
            normalised_value=None,
            status="invalid",
            reason=(
                f"{table_name}.{column_name}[{row_position}] "
                "must contain non-empty string values"
            ),
        )
    normalised = stripped.upper()
    status = "unchanged" if normalised == original_value else "normalised"
    return ReferenceIdentifierNormalisationRecord(
        table_name=table_name,
        column_name=column_name,
        row_position=row_position,
        identifier_kind="kinase",
        original_value=original_value,
        normalised_value=normalised,
        status=status,
        reason=None,
    )


def normalise_reference_protein_accession(
    value: object,
    *,
    table_name: str,
    column_name: str,
    row_position: int,
) -> ReferenceIdentifierNormalisationRecord:
    original_value = _stringify_identifier(value)
    if _is_missing(value):
        return ReferenceIdentifierNormalisationRecord(
            table_name=table_name,
            column_name=column_name,
            row_position=row_position,
            identifier_kind="protein_accession",
            original_value=original_value,
            normalised_value=None,
            status="invalid",
            reason=(
                f"{table_name}.{column_name}[{row_position}] "
                "must not contain missing values"
            ),
        )
    if not isinstance(value, str):
        return ReferenceIdentifierNormalisationRecord(
            table_name=table_name,
            column_name=column_name,
            row_position=row_position,
            identifier_kind="protein_accession",
            original_value=original_value,
            normalised_value=None,
            status="invalid",
            reason=(
                f"{table_name}.{column_name}[{row_position}] "
                "must contain non-empty string values"
            ),
        )
    normalised = value.strip()
    if normalised == "":
        return ReferenceIdentifierNormalisationRecord(
            table_name=table_name,
            column_name=column_name,
            row_position=row_position,
            identifier_kind="protein_accession",
            original_value=original_value,
            normalised_value=None,
            status="invalid",
            reason=(
                f"{table_name}.{column_name}[{row_position}] "
                "must contain non-empty string values"
            ),
        )
    status = "unchanged" if normalised == original_value else "normalised"
    return ReferenceIdentifierNormalisationRecord(
        table_name=table_name,
        column_name=column_name,
        row_position=row_position,
        identifier_kind="protein_accession",
        original_value=original_value,
        normalised_value=normalised,
        status=status,
        reason=None,
    )


def _stringify_identifier(value: object) -> str:
    if value is None:
        return ""
    return str(value)


def _is_missing(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


__all__ = [
    "REFERENCE_IDENTIFIER_NORMALISATION_SCHEMA_VERSION",
    "ReferenceIdentifierNormalisationRecord",
    "ReferenceIdentifierNormalisationReport",
    "build_reference_identifier_normalisation_report",
    "merge_reference_identifier_normalisation_reports",
    "normalise_reference_kinase_id",
    "normalise_reference_protein_accession",
    "normalise_reference_site_id",
]
