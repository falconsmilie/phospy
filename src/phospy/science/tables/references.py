"""Reference scientific table wrappers."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import NoReturn

import pandas as pd

from phospy.errors.validation import (
    ReferenceIdentifierNormalisationValidationError,
    ReferenceValidationError,
)
from phospy.frames.table_schema import TableSchema
from phospy.frames.validation import (
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_unique_index,
)
from phospy.science.references.identifiers import (
    ReferenceIdentifierNormalisationRecord,
    ReferenceIdentifierNormalisationReport,
    build_reference_identifier_normalisation_report,
    normalise_reference_kinase_id,
    normalise_reference_protein_accession,
    normalise_reference_site_id,
)


def _raise_with_identifier_normalisation_report(
    *,
    message: str,
    report: ReferenceIdentifierNormalisationReport,
) -> NoReturn:
    raise ReferenceIdentifierNormalisationValidationError(
        message=message,
        identifier_normalisation_report=report,
    )


@dataclass(frozen=True, slots=True)
class KinaseSubstrateReference(TableSchema):
    """Schema wrapper for ``references.kinase_substrate_map``.

    Validation here is structural and content-based. It does not enforce
    workflow-scale limits; performance behavior is owned by workflow overlap
    filtering and scoring-lane selection.
    """

    _field_name = "references.kinase_substrate_map"
    _error_type = ReferenceValidationError
    identifier_normalisation: ReferenceIdentifierNormalisationReport | None = field(
        init=False,
        default=None,
    )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        require_columns(
            frame,
            field_name=self._field_name,
            required_columns=("kinase", "substrate_site"),
            error_type=self._error_type,
        )
        records: list[ReferenceIdentifierNormalisationRecord] = []
        canonical_kinase: list[str | None] = []
        canonical_substrate_site: list[str | None] = []
        for row_position, (kinase_raw, site_raw) in enumerate(
            zip(
                frame.loc[:, "kinase"].tolist(),
                frame.loc[:, "substrate_site"].tolist(),
                strict=False,
            )
        ):
            kinase_record = normalise_reference_kinase_id(
                kinase_raw,
                table_name=self._field_name,
                column_name="kinase",
                row_position=row_position,
            )
            site_record = normalise_reference_site_id(
                site_raw,
                table_name=self._field_name,
                column_name="substrate_site",
                row_position=row_position,
            )
            records.extend((kinase_record, site_record))
            canonical_kinase.append(kinase_record.normalised_value)
            canonical_substrate_site.append(site_record.normalised_value)
        valid_rows = int(
            sum(
                1
                for kinase_value, site_value in zip(
                    canonical_kinase,
                    canonical_substrate_site,
                    strict=False,
                )
                if kinase_value is not None and site_value is not None
            )
        )
        report = build_reference_identifier_normalisation_report(
            original_row_count=int(frame.shape[0]),
            normalised_row_count=valid_rows,
            records=records,
        )
        object.__setattr__(self, "identifier_normalisation", report)
        invalid_records = [record for record in records if record.status == "invalid"]
        if invalid_records:
            _raise_with_identifier_normalisation_report(
                message=invalid_records[0].reason or "invalid identifier",
                report=report,
            )
        frame = frame.copy()
        frame.loc[:, "kinase"] = pd.Series(
            [value for value in canonical_kinase if value is not None],
            index=frame.index.copy(),
            dtype="string",
        )
        frame.loc[:, "substrate_site"] = pd.Series(
            [value for value in canonical_substrate_site if value is not None],
            index=frame.index.copy(),
            dtype="string",
        )
        duplicated = frame.duplicated(
            subset=["kinase", "substrate_site"],
            keep=False,
        )
        if bool(duplicated.any()):
            duplicate_records, conflict_records, conflicting_pairs = (
                _classify_duplicate_and_conflicting_pair_records(
                    frame=frame,
                    duplicated=duplicated,
                    existing_records=records,
                    table_name=self._field_name,
                )
            )
            report = build_reference_identifier_normalisation_report(
                original_row_count=int(frame.shape[0]),
                normalised_row_count=int(frame.shape[0]),
                records=[*records, *duplicate_records, *conflict_records],
                duplicate_identifier_count=len(duplicate_records),
                conflict_count=len(conflict_records),
            )
            object.__setattr__(self, "identifier_normalisation", report)
            if conflicting_pairs:
                preview = ", ".join(repr(pair) for pair in conflicting_pairs[:5])
                suffix = "" if len(conflicting_pairs) <= 5 else " ..."
                _raise_with_identifier_normalisation_report(
                    message=(
                        f"{self._field_name} contains conflicting payload rows for "
                        f"normalised (kinase, substrate_site) pairs: {preview}{suffix}"
                    ),
                    report=report,
                )
            duplicate_pairs = list(
                frame.loc[duplicated, ["kinase", "substrate_site"]]
                .drop_duplicates()
                .itertuples(index=False, name=None)
            )
            duplicate_preview = ", ".join(repr(pair) for pair in duplicate_pairs[:5])
            suffix = "" if len(duplicate_pairs) <= 5 else " ..."
            _raise_with_identifier_normalisation_report(
                message=(
                    f"{self._field_name} contains duplicate (kinase, substrate_site) "
                    f"pairs: {duplicate_preview}{suffix}"
                ),
                report=report,
            )
        return frame


@dataclass(frozen=True, slots=True)
class SiteSequenceReference(TableSchema):
    """Schema wrapper for ``references.site_sequences``."""

    _field_name = "references.site_sequences"
    _error_type = ReferenceValidationError
    identifier_normalisation: ReferenceIdentifierNormalisationReport | None = field(
        init=False,
        default=None,
    )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        require_columns(
            frame,
            field_name=self._field_name,
            required_columns=("site_sequence",),
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="site_sequence",
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name="site_sequence",
            error_type=self._error_type,
        )
        records: list[ReferenceIdentifierNormalisationRecord] = []
        canonical_index: list[str | None] = []
        for row_position, raw_value in enumerate(frame.index.tolist()):
            record = normalise_reference_site_id(
                raw_value,
                table_name=self._field_name,
                column_name="index",
                row_position=row_position,
            )
            records.append(record)
            canonical_index.append(record.normalised_value)
        valid_rows = int(sum(1 for value in canonical_index if value is not None))
        report = build_reference_identifier_normalisation_report(
            original_row_count=int(frame.shape[0]),
            normalised_row_count=valid_rows,
            records=records,
        )
        object.__setattr__(self, "identifier_normalisation", report)
        invalid_records = [record for record in records if record.status == "invalid"]
        if invalid_records:
            _raise_with_identifier_normalisation_report(
                message=invalid_records[0].reason or "invalid identifier",
                report=report,
            )
        frame = frame.copy()
        frame.index = pd.Index(
            [value for value in canonical_index if value is not None],
            name=frame.index.name,
        )
        duplicated = pd.Series(
            frame.index.duplicated(keep=False),
            index=frame.index.copy(),
            dtype="bool",
        )
        if bool(duplicated.any()):
            duplicate_records, conflict_records, duplicate_values, conflict_values = (
                _classify_duplicate_and_conflicting_index_records(
                    frame=frame,
                    duplicated=duplicated,
                    existing_records=records,
                    table_name=self._field_name,
                )
            )
            report = build_reference_identifier_normalisation_report(
                original_row_count=int(frame.shape[0]),
                normalised_row_count=int(frame.shape[0]),
                records=[*records, *duplicate_records, *conflict_records],
                duplicate_identifier_count=len(duplicate_records),
                conflict_count=len(conflict_records),
            )
            object.__setattr__(self, "identifier_normalisation", report)
            if conflict_values:
                preview = ", ".join(conflict_values[:5])
                suffix = "" if len(conflict_values) <= 5 else " ..."
                _raise_with_identifier_normalisation_report(
                    message=(
                        f"{self._field_name}.index contains conflicting site_sequence "
                        f"values after normalisation: {preview}{suffix}"
                    ),
                    report=report,
                )
            duplicate_preview = ", ".join(duplicate_values[:5])
            suffix = "" if len(duplicate_values) <= 5 else " ..."
            _raise_with_identifier_normalisation_report(
                message=(
                    f"{self._field_name}.index contains duplicate site identifiers "
                    f"after normalisation: {duplicate_preview}{suffix}"
                ),
                report=report,
            )
        require_unique_index(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        return frame


@dataclass(frozen=True, slots=True)
class ProteinAccessionReference(TableSchema):
    """Schema wrapper for ``references.protein_accessions``."""

    _field_name = "references.protein_accessions"
    _error_type = ReferenceValidationError
    identifier_normalisation: ReferenceIdentifierNormalisationReport | None = field(
        init=False,
        default=None,
    )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=False,
            error_type=self._error_type,
        )
        require_columns(
            frame,
            field_name=self._field_name,
            required_columns=("protein_accession",),
            error_type=self._error_type,
        )
        records: list[ReferenceIdentifierNormalisationRecord] = []
        canonical_accessions: list[str | None] = []
        for row_position, raw_value in enumerate(
            frame.loc[:, "protein_accession"].tolist()
        ):
            record = normalise_reference_protein_accession(
                raw_value,
                table_name=self._field_name,
                column_name="protein_accession",
                row_position=row_position,
            )
            records.append(record)
            canonical_accessions.append(record.normalised_value)
        valid_rows = int(sum(1 for value in canonical_accessions if value is not None))
        report = build_reference_identifier_normalisation_report(
            original_row_count=int(frame.shape[0]),
            normalised_row_count=valid_rows,
            records=records,
        )
        object.__setattr__(self, "identifier_normalisation", report)
        invalid_records = [record for record in records if record.status == "invalid"]
        if invalid_records:
            _raise_with_identifier_normalisation_report(
                message=invalid_records[0].reason or "invalid identifier",
                report=report,
            )
        frame = frame.copy(deep=True)
        frame.loc[:, "protein_accession"] = pd.Series(
            [value for value in canonical_accessions if value is not None],
            index=frame.index.copy(),
            dtype="string",
        )

        duplicated = frame.duplicated(
            subset=["protein_accession"],
            keep=False,
        )
        if not bool(duplicated.any()):
            return frame

        duplicate_records, conflict_records = (
            _classify_duplicate_and_conflicting_protein_accession_records(
                frame=frame,
                duplicated=duplicated,
                existing_records=records,
                table_name=self._field_name,
            )
        )
        report = build_reference_identifier_normalisation_report(
            original_row_count=int(frame.shape[0]),
            normalised_row_count=int(frame.shape[0]),
            records=[*records, *duplicate_records, *conflict_records],
            duplicate_identifier_count=len(duplicate_records),
            conflict_count=len(conflict_records),
        )
        object.__setattr__(self, "identifier_normalisation", report)
        if conflict_records:
            _raise_with_identifier_normalisation_report(
                message=conflict_records[0].reason or "conflicting payload",
                report=report,
            )
        _raise_with_identifier_normalisation_report(
            message=duplicate_records[0].reason or "duplicate identifier",
            report=report,
        )


def _classify_duplicate_and_conflicting_pair_records(
    *,
    frame: pd.DataFrame,
    duplicated: pd.Series,
    existing_records: list[ReferenceIdentifierNormalisationRecord],
    table_name: str,
) -> tuple[
    tuple[ReferenceIdentifierNormalisationRecord, ...],
    tuple[ReferenceIdentifierNormalisationRecord, ...],
    tuple[tuple[str, str], ...],
]:
    duplicate_reasons_by_row: dict[int, str] = {}
    conflict_reasons_by_row: dict[int, str] = {}
    conflict_pairs: list[tuple[str, str]] = []
    duplicate_rows = frame.loc[duplicated, :].copy()
    duplicate_rows.loc[:, "_row_position"] = [
        int(position)
        for position, is_duplicate in enumerate(duplicated.tolist())
        if is_duplicate
    ]
    payload_columns = [
        column
        for column in frame.columns.tolist()
        if column not in ("kinase", "substrate_site")
    ]
    for _, grouped in duplicate_rows.groupby(["kinase", "substrate_site"], sort=False):
        kinase = str(grouped.iloc[0]["kinase"])
        substrate_site = str(grouped.iloc[0]["substrate_site"])
        row_positions = grouped.loc[:, "_row_position"].astype(int).tolist()
        has_conflicting_payload = (
            _group_has_conflicting_payload_rows(grouped, payload_columns)
            if payload_columns
            else False
        )
        if has_conflicting_payload:
            conflict_pairs.append((kinase, substrate_site))
            reason = (
                "conflicting payload for (kinase, substrate_site) pair after "
                f"normalisation: ({kinase!r}, {substrate_site!r})"
            )
            for row_position in row_positions:
                conflict_reasons_by_row[int(row_position)] = reason
            continue
        reason = (
            "duplicate (kinase, substrate_site) pair after normalisation: "
            f"({kinase!r}, {substrate_site!r})"
        )
        for row_position in row_positions:
            duplicate_reasons_by_row[int(row_position)] = reason
    latest_by_key = {
        (record.row_position, record.column_name): record for record in existing_records
    }
    duplicate_records = _build_pair_classification_records(
        reasons_by_row=duplicate_reasons_by_row,
        existing_records=latest_by_key,
        table_name=table_name,
        status="duplicate_after_normalisation",
    )
    conflict_records = _build_pair_classification_records(
        reasons_by_row=conflict_reasons_by_row,
        existing_records=latest_by_key,
        table_name=table_name,
        status="conflict_after_normalisation",
    )
    return (
        duplicate_records,
        conflict_records,
        tuple(conflict_pairs),
    )


def _group_has_conflicting_payload_rows(
    grouped: pd.DataFrame,
    payload_columns: Sequence[str],
) -> bool:
    payload_rows = grouped.loc[:, payload_columns]
    anchor = payload_rows.iloc[0]
    for row_position in range(1, int(payload_rows.shape[0])):
        if not payload_rows.iloc[row_position].equals(anchor):
            return True
    return False


def _build_pair_classification_records(
    *,
    reasons_by_row: dict[int, str],
    existing_records: dict[tuple[int, str], ReferenceIdentifierNormalisationRecord],
    table_name: str,
    status: str,
) -> tuple[ReferenceIdentifierNormalisationRecord, ...]:
    classified_records: list[ReferenceIdentifierNormalisationRecord] = []
    for row_position, reason in reasons_by_row.items():
        for column_name in ("kinase", "substrate_site"):
            source = existing_records[(row_position, column_name)]
            classified_records.append(
                ReferenceIdentifierNormalisationRecord(
                    table_name=table_name,
                    column_name=column_name,
                    row_position=row_position,
                    identifier_kind=source.identifier_kind,
                    original_value=source.original_value,
                    normalised_value=source.normalised_value,
                    status=status,
                    reason=reason,
                )
            )
    return tuple(classified_records)


def _classify_duplicate_and_conflicting_index_records(
    *,
    frame: pd.DataFrame,
    duplicated: pd.Series,
    existing_records: list[ReferenceIdentifierNormalisationRecord],
    table_name: str,
) -> tuple[
    tuple[ReferenceIdentifierNormalisationRecord, ...],
    tuple[ReferenceIdentifierNormalisationRecord, ...],
    tuple[str, ...],
    tuple[str, ...],
]:
    duplicate_reasons_by_row: dict[int, str] = {}
    conflict_reasons_by_row: dict[int, str] = {}
    duplicate_values: list[str] = []
    conflict_values: list[str] = []
    duplicate_rows = frame.loc[duplicated, ["site_sequence"]].copy()
    duplicate_rows.loc[:, "_row_position"] = [
        int(position)
        for position, is_duplicate in enumerate(duplicated.tolist())
        if is_duplicate
    ]
    duplicate_rows.loc[:, "_site_id"] = frame.index[duplicated].tolist()
    for site_id, grouped in duplicate_rows.groupby("_site_id", sort=False):
        row_positions = grouped.loc[:, "_row_position"].astype(int).tolist()
        site_sequences = grouped.loc[:, "site_sequence"].drop_duplicates().tolist()
        normalised_site_id = str(site_id)
        if len(site_sequences) > 1:
            conflict_values.append(normalised_site_id)
            reason = (
                "conflicting site_sequence values for site identifier after "
                f"normalisation: {normalised_site_id!r}"
            )
            for row_position in row_positions:
                conflict_reasons_by_row[int(row_position)] = reason
            continue
        duplicate_values.append(normalised_site_id)
        reason = (
            "duplicate site identifier after normalisation with identical "
            f"site_sequence: {normalised_site_id!r}"
        )
        for row_position in row_positions:
            duplicate_reasons_by_row[int(row_position)] = reason
    latest_by_row = {
        record.row_position: record
        for record in existing_records
        if record.column_name == "index"
    }
    duplicate_records = _build_index_classification_records(
        reasons_by_row=duplicate_reasons_by_row,
        latest_by_row=latest_by_row,
        table_name=table_name,
        frame=frame,
        status="duplicate_after_normalisation",
    )
    conflict_records = _build_index_classification_records(
        reasons_by_row=conflict_reasons_by_row,
        latest_by_row=latest_by_row,
        table_name=table_name,
        frame=frame,
        status="conflict_after_normalisation",
    )
    return (
        duplicate_records,
        conflict_records,
        tuple(dict.fromkeys(duplicate_values)),
        tuple(dict.fromkeys(conflict_values)),
    )


def _build_index_classification_records(
    *,
    reasons_by_row: dict[int, str],
    latest_by_row: dict[int, ReferenceIdentifierNormalisationRecord],
    table_name: str,
    frame: pd.DataFrame,
    status: str,
) -> tuple[ReferenceIdentifierNormalisationRecord, ...]:
    classified_records: list[ReferenceIdentifierNormalisationRecord] = []
    for row_position, reason in reasons_by_row.items():
        source = latest_by_row[row_position]
        classified_records.append(
            ReferenceIdentifierNormalisationRecord(
                table_name=table_name,
                column_name="index",
                row_position=row_position,
                identifier_kind="site_id",
                original_value=source.original_value,
                normalised_value=str(frame.index[row_position]),
                status=status,
                reason=reason,
            )
        )
    return tuple(classified_records)


def _classify_duplicate_and_conflicting_protein_accession_records(
    *,
    frame: pd.DataFrame,
    duplicated: pd.Series,
    existing_records: list[ReferenceIdentifierNormalisationRecord],
    table_name: str,
) -> tuple[
    tuple[ReferenceIdentifierNormalisationRecord, ...],
    tuple[ReferenceIdentifierNormalisationRecord, ...],
]:
    duplicate_reasons_by_row: dict[int, str] = {}
    conflict_reasons_by_row: dict[int, str] = {}
    duplicate_rows = frame.loc[duplicated, :].copy()
    duplicate_rows.loc[:, "_row_position"] = [
        int(position)
        for position, is_duplicate in enumerate(duplicated.tolist())
        if is_duplicate
    ]
    payload_columns = [
        column for column in frame.columns.tolist() if column != "protein_accession"
    ]
    for accession_value, grouped in duplicate_rows.groupby(
        "protein_accession",
        sort=False,
    ):
        row_positions = grouped.loc[:, "_row_position"].astype(int).tolist()
        has_conflicting_payload = (
            _group_has_conflicting_payload_rows(grouped, payload_columns)
            if payload_columns
            else False
        )
        if has_conflicting_payload:
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

    latest_by_row = {
        record.row_position: record
        for record in existing_records
        if record.column_name == "protein_accession"
    }
    duplicate_records = _build_protein_accession_classification_records(
        reasons_by_row=duplicate_reasons_by_row,
        latest_by_row=latest_by_row,
        table_name=table_name,
        status="duplicate_after_normalisation",
    )
    conflict_records = _build_protein_accession_classification_records(
        reasons_by_row=conflict_reasons_by_row,
        latest_by_row=latest_by_row,
        table_name=table_name,
        status="conflict_after_normalisation",
    )
    return duplicate_records, conflict_records


def _build_protein_accession_classification_records(
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
                column_name="protein_accession",
                row_position=row_position,
                identifier_kind=source.identifier_kind,
                original_value=source.original_value,
                normalised_value=source.normalised_value,
                status=status,
                reason=reason,
            )
        )
    return tuple(classified_records)
