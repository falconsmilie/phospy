"""Reference scientific table wrappers."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.references.identifiers import (
    ReferenceIdentifierNormalisationRecord,
    ReferenceIdentifierNormalisationReport,
    build_reference_identifier_normalisation_report,
    normalise_reference_kinase_id,
    normalise_reference_site_id,
)
from phospy.tables.base import TableSchema
from phospy.validation.common.dataframes import (
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_unique_index,
    require_unique_row_pairs,
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
            raise self._error_type(invalid_records[0].reason or "invalid identifier")
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
            duplicate_records = _duplicate_pair_records(
                frame=frame,
                duplicated=duplicated,
                existing_records=records,
                table_name=self._field_name,
            )
            report = build_reference_identifier_normalisation_report(
                original_row_count=int(frame.shape[0]),
                normalised_row_count=int(frame.shape[0]),
                records=[*records, *duplicate_records],
                duplicate_identifier_count=len(duplicate_records),
            )
            object.__setattr__(self, "identifier_normalisation", report)
        require_unique_row_pairs(
            frame,
            field_name=self._field_name,
            column_names=("kinase", "substrate_site"),
            error_type=self._error_type,
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
            raise self._error_type(invalid_records[0].reason or "invalid identifier")
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
            duplicate_records = _duplicate_index_records(
                frame=frame,
                duplicated=duplicated,
                existing_records=records,
                table_name=self._field_name,
            )
            report = build_reference_identifier_normalisation_report(
                original_row_count=int(frame.shape[0]),
                normalised_row_count=int(frame.shape[0]),
                records=[*records, *duplicate_records],
                duplicate_identifier_count=len(duplicate_records),
            )
            object.__setattr__(self, "identifier_normalisation", report)
            duplicate_values = list(dict.fromkeys(frame.index[duplicated].tolist()))
            preview = ", ".join(duplicate_values[:5])
            suffix = "" if len(duplicate_values) <= 5 else " ..."
            raise self._error_type(
                f"{self._field_name}.index contains duplicate site identifiers after "
                f"canonicalization: {preview}{suffix}"
            )
        require_unique_index(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        return frame


def _duplicate_pair_records(
    *,
    frame: pd.DataFrame,
    duplicated: pd.Series,
    existing_records: list[ReferenceIdentifierNormalisationRecord],
    table_name: str,
) -> tuple[ReferenceIdentifierNormalisationRecord, ...]:
    reasons_by_row: dict[int, str] = {}
    duplicate_rows = frame.loc[duplicated, ["kinase", "substrate_site"]].copy()
    duplicate_rows.loc[:, "_row_position"] = [
        int(position)
        for position, is_duplicate in enumerate(duplicated.tolist())
        if is_duplicate
    ]
    for _, grouped in duplicate_rows.groupby(["kinase", "substrate_site"], sort=False):
        kinase = str(grouped.iloc[0]["kinase"])
        substrate_site = str(grouped.iloc[0]["substrate_site"])
        reason = (
            "duplicate (kinase, substrate_site) pair after normalisation: "
            f"({kinase!r}, {substrate_site!r})"
        )
        for row_position in grouped.loc[:, "_row_position"].tolist():
            reasons_by_row[int(row_position)] = reason
    latest_by_key = {
        (record.row_position, record.column_name): record for record in existing_records
    }
    duplicate_records: list[ReferenceIdentifierNormalisationRecord] = []
    for row_position, reason in reasons_by_row.items():
        for column_name in ("kinase", "substrate_site"):
            source = latest_by_key[(row_position, column_name)]
            duplicate_records.append(
                ReferenceIdentifierNormalisationRecord(
                    table_name=table_name,
                    column_name=column_name,
                    row_position=row_position,
                    identifier_kind=source.identifier_kind,
                    original_value=source.original_value,
                    normalised_value=source.normalised_value,
                    status="duplicate_after_normalisation",
                    reason=reason,
                )
            )
    return tuple(duplicate_records)


def _duplicate_index_records(
    *,
    frame: pd.DataFrame,
    duplicated: pd.Series,
    existing_records: list[ReferenceIdentifierNormalisationRecord],
    table_name: str,
) -> tuple[ReferenceIdentifierNormalisationRecord, ...]:
    duplicate_records: list[ReferenceIdentifierNormalisationRecord] = []
    latest_by_row = {
        record.row_position: record
        for record in existing_records
        if record.column_name == "index"
    }
    for row_position, is_duplicate in enumerate(duplicated.tolist()):
        if not is_duplicate:
            continue
        normalised_value = frame.index[row_position]
        source = latest_by_row[row_position]
        duplicate_records.append(
            ReferenceIdentifierNormalisationRecord(
                table_name=table_name,
                column_name="index",
                row_position=row_position,
                identifier_kind="site_id",
                original_value=source.original_value,
                normalised_value=str(normalised_value),
                status="duplicate_after_normalisation",
                reason=(
                    "duplicate site identifier after normalisation: "
                    f"{normalised_value!r}"
                ),
            )
        )
    return tuple(duplicate_records)
