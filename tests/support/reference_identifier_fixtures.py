from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.validation import (
    ReferenceIdentifierNormalisationValidationError,
)
from phospy.science.references.identifiers import (
    ReferenceIdentifierNormalisationRecord,
    ReferenceIdentifierNormalisationReport,
    build_reference_identifier_normalisation_report,
    normalise_reference_protein_accession,
)

_PROTEIN_ACCESSION_COLUMN = "protein_accession"


@dataclass(frozen=True, slots=True)
class ProteinAccessionReferenceFixtureResult:
    frame: pd.DataFrame
    report: ReferenceIdentifierNormalisationReport


def normalise_explicit_reference_protein_accession_fixture(
    frame: pd.DataFrame,
    *,
    table_name: str = "tests.reference.protein_accession_fixture",
) -> ProteinAccessionReferenceFixtureResult:
    """Test-only explicit-reference fixture for protein accession identifiers.

    Production ingestion is owned by ``ProteinAccessionReference`` in
    ``phospy.tables.references``. This fixture remains for focused collaborator
    tests that need direct helper-path control.
    """

    records: list[ReferenceIdentifierNormalisationRecord] = []
    normalised_values: list[str | None] = []
    for row_position, raw_value in enumerate(
        frame.loc[:, _PROTEIN_ACCESSION_COLUMN].tolist()
    ):
        record = normalise_reference_protein_accession(
            raw_value,
            table_name=table_name,
            column_name=_PROTEIN_ACCESSION_COLUMN,
            row_position=row_position,
        )
        records.append(record)
        normalised_values.append(record.normalised_value)

    valid_row_count = int(sum(value is not None for value in normalised_values))
    report = build_reference_identifier_normalisation_report(
        original_row_count=int(frame.shape[0]),
        normalised_row_count=valid_row_count,
        records=records,
    )
    invalid_records = [record for record in records if record.status == "invalid"]
    if invalid_records:
        raise ReferenceIdentifierNormalisationValidationError(
            message=invalid_records[0].reason or "invalid identifier",
            identifier_normalisation_report=report,
        )

    normalised_frame = frame.copy(deep=True)
    normalised_frame.loc[:, _PROTEIN_ACCESSION_COLUMN] = pd.Series(
        [value for value in normalised_values if value is not None],
        index=normalised_frame.index.copy(),
        dtype="string",
    )

    duplicated = normalised_frame.duplicated(
        subset=[_PROTEIN_ACCESSION_COLUMN],
        keep=False,
    )
    if not bool(duplicated.any()):
        return ProteinAccessionReferenceFixtureResult(
            frame=normalised_frame,
            report=report,
        )

    duplicate_reasons_by_row: dict[int, str] = {}
    conflict_reasons_by_row: dict[int, str] = {}
    duplicate_rows = normalised_frame.loc[duplicated, :].copy()
    duplicate_rows.loc[:, "_row_position"] = [
        int(position)
        for position, is_duplicate in enumerate(duplicated.tolist())
        if is_duplicate
    ]
    payload_columns = [
        column
        for column in normalised_frame.columns.tolist()
        if column != _PROTEIN_ACCESSION_COLUMN
    ]
    for accession_value, grouped in duplicate_rows.groupby(
        _PROTEIN_ACCESSION_COLUMN,
        sort=False,
    ):
        row_positions = grouped.loc[:, "_row_position"].astype(int).tolist()
        if payload_columns:
            payload_rows = grouped.loc[:, payload_columns]
            anchor = payload_rows.iloc[0]
            has_conflict = any(
                not payload_rows.iloc[row_position].equals(anchor)
                for row_position in range(1, int(payload_rows.shape[0]))
            )
        else:
            has_conflict = False
        if has_conflict:
            reason = (
                "conflicting payload for protein_accession after normalisation: "
                f"{str(accession_value)!r}"
            )
            for row_position in row_positions:
                conflict_reasons_by_row[int(row_position)] = reason
            continue
        reason = (
            f"duplicate protein_accession after normalisation: {str(accession_value)!r}"
        )
        for row_position in row_positions:
            duplicate_reasons_by_row[int(row_position)] = reason

    latest_by_row = {record.row_position: record for record in records}
    duplicate_records = _build_fixture_classification_records(
        reasons_by_row=duplicate_reasons_by_row,
        latest_by_row=latest_by_row,
        table_name=table_name,
        status="duplicate_after_normalisation",
    )
    conflict_records = _build_fixture_classification_records(
        reasons_by_row=conflict_reasons_by_row,
        latest_by_row=latest_by_row,
        table_name=table_name,
        status="conflict_after_normalisation",
    )
    report_with_duplicates = build_reference_identifier_normalisation_report(
        original_row_count=int(normalised_frame.shape[0]),
        normalised_row_count=int(normalised_frame.shape[0]),
        records=[*records, *duplicate_records, *conflict_records],
        duplicate_identifier_count=len(duplicate_records),
        conflict_count=len(conflict_records),
    )
    if conflict_records:
        raise ReferenceIdentifierNormalisationValidationError(
            message=conflict_records[0].reason or "conflicting payload after trim",
            identifier_normalisation_report=report_with_duplicates,
        )
    raise ReferenceIdentifierNormalisationValidationError(
        message=duplicate_records[0].reason or "duplicate identifiers after trim",
        identifier_normalisation_report=report_with_duplicates,
    )


def _build_fixture_classification_records(
    *,
    reasons_by_row: dict[int, str],
    latest_by_row: dict[int, ReferenceIdentifierNormalisationRecord],
    table_name: str,
    status: str,
) -> tuple[ReferenceIdentifierNormalisationRecord, ...]:
    classified_records: list[ReferenceIdentifierNormalisationRecord] = []
    for row_position, reason in reasons_by_row.items():
        source = latest_by_row[row_position]
        classified_records.append(
            ReferenceIdentifierNormalisationRecord(
                table_name=table_name,
                column_name=_PROTEIN_ACCESSION_COLUMN,
                row_position=row_position,
                identifier_kind=source.identifier_kind,
                original_value=source.original_value,
                normalised_value=source.normalised_value,
                status=status,
                reason=reason,
            )
        )
    return tuple(classified_records)


__all__ = [
    "ProteinAccessionReferenceFixtureResult",
    "normalise_explicit_reference_protein_accession_fixture",
]
