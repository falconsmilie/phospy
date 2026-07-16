"""Leaf provenance records for reference identifier normalisation."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

REFERENCE_IDENTIFIER_NORMALISATION_SCHEMA_VERSION = 1
_IDENTIFIER_KINDS = frozenset({"site_id", "kinase", "protein_accession"})
_STATUSES = frozenset(
    {
        "unchanged",
        "normalised",
        "invalid",
        "duplicate_after_normalisation",
        "conflict_after_normalisation",
    }
)


@dataclass(frozen=True, slots=True)
class ReferenceIdentifierNormalisationRecord:
    table_name: str
    column_name: str
    row_position: int
    identifier_kind: str
    original_value: str
    normalised_value: str | None
    status: str
    reason: str | None

    def __post_init__(self) -> None:
        if self.identifier_kind not in _IDENTIFIER_KINDS:
            raise ValueError(
                "reference identifier normalisation record identifier_kind "
                f"must be one of {sorted(_IDENTIFIER_KINDS)}"
            )
        if self.status not in _STATUSES:
            raise ValueError(
                "reference identifier normalisation record status "
                f"must be one of {sorted(_STATUSES)}"
            )


@dataclass(frozen=True, slots=True)
class ReferenceIdentifierNormalisationReport:
    schema_version: int
    original_row_count: int
    normalised_row_count: int
    invalid_identifier_count: int
    changed_identifier_count: int
    duplicate_identifier_count: int
    conflict_count: int
    records: tuple[ReferenceIdentifierNormalisationRecord, ...]


def build_reference_identifier_normalisation_report(
    *,
    original_row_count: int,
    normalised_row_count: int,
    records: Iterable[ReferenceIdentifierNormalisationRecord],
    duplicate_identifier_count: int = 0,
    conflict_count: int = 0,
) -> ReferenceIdentifierNormalisationReport:
    resolved_records = tuple(records)
    return ReferenceIdentifierNormalisationReport(
        schema_version=REFERENCE_IDENTIFIER_NORMALISATION_SCHEMA_VERSION,
        original_row_count=int(original_row_count),
        normalised_row_count=int(normalised_row_count),
        invalid_identifier_count=sum(
            1 for record in resolved_records if record.status == "invalid"
        ),
        changed_identifier_count=sum(
            1 for record in resolved_records if record.status == "normalised"
        ),
        duplicate_identifier_count=int(duplicate_identifier_count),
        conflict_count=int(conflict_count),
        records=resolved_records,
    )


def merge_reference_identifier_normalisation_reports(
    reports: Iterable[ReferenceIdentifierNormalisationReport],
) -> ReferenceIdentifierNormalisationReport | None:
    resolved = tuple(reports)
    if not resolved:
        return None
    return ReferenceIdentifierNormalisationReport(
        schema_version=REFERENCE_IDENTIFIER_NORMALISATION_SCHEMA_VERSION,
        original_row_count=sum(item.original_row_count for item in resolved),
        normalised_row_count=sum(item.normalised_row_count for item in resolved),
        invalid_identifier_count=sum(
            item.invalid_identifier_count for item in resolved
        ),
        changed_identifier_count=sum(
            item.changed_identifier_count for item in resolved
        ),
        duplicate_identifier_count=sum(
            item.duplicate_identifier_count for item in resolved
        ),
        conflict_count=sum(item.conflict_count for item in resolved),
        records=tuple(record for item in resolved for record in item.records),
    )


__all__ = [
    "REFERENCE_IDENTIFIER_NORMALISATION_SCHEMA_VERSION",
    "ReferenceIdentifierNormalisationRecord",
    "ReferenceIdentifierNormalisationReport",
    "build_reference_identifier_normalisation_report",
    "merge_reference_identifier_normalisation_reports",
]
