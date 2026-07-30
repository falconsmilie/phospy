"""Table and row-count provenance models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import freeze_optional_json_mapping
from phospy.provenance.models._shared import (
    JsonValue,
    _provenance_string_tuple,
    _required_non_negative_row_count,
    _required_provenance_text,
    _required_provenance_text_tuple,
)


@dataclass(frozen=True, slots=True)
class RowAttritionRecord:
    """Standard count-only provenance for rows removed at one workflow stage."""

    stage: str
    input_rows: int
    output_rows: int
    removed_rows: int
    reason: str
    examples: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        input_rows = _required_non_negative_row_count(
            self.input_rows,
            field_name="row_attrition_record.input_rows",
        )
        output_rows = _required_non_negative_row_count(
            self.output_rows,
            field_name="row_attrition_record.output_rows",
        )
        removed_rows = _required_non_negative_row_count(
            self.removed_rows,
            field_name="row_attrition_record.removed_rows",
        )
        if output_rows > input_rows:
            raise PhosPyInputError(
                "row_attrition_record.output_rows must be less than or equal "
                "to input_rows"
            )
        expected_removed_rows = input_rows - output_rows
        if removed_rows != expected_removed_rows:
            raise PhosPyInputError(
                "row_attrition_record.removed_rows must equal input_rows - output_rows"
            )
        object.__setattr__(
            self,
            "stage",
            _required_provenance_text(
                self.stage,
                field_name="row_attrition_record.stage",
            ),
        )
        object.__setattr__(self, "input_rows", input_rows)
        object.__setattr__(self, "output_rows", output_rows)
        object.__setattr__(self, "removed_rows", removed_rows)
        object.__setattr__(
            self,
            "reason",
            _required_provenance_text(
                self.reason,
                field_name="row_attrition_record.reason",
            ),
        )
        object.__setattr__(
            self,
            "examples",
            _required_provenance_text_tuple(
                self.examples,
                field_name="row_attrition_record.examples",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible row-attrition record payload."""

        return {
            "stage": self.stage,
            "input_rows": int(self.input_rows),
            "output_rows": int(self.output_rows),
            "removed_rows": int(self.removed_rows),
            "reason": self.reason,
            "examples": list(self.examples),
        }


@dataclass(frozen=True, slots=True)
class RowAttritionReport:
    """Ordered row-attrition provenance across workflow stages."""

    records: tuple[RowAttritionRecord, ...]
    input_rows: int
    final_rows: int

    def __post_init__(self) -> None:
        records = _required_row_attrition_record_tuple(self.records)
        input_rows = _required_non_negative_row_count(
            self.input_rows,
            field_name="row_attrition_report.input_rows",
        )
        final_rows = _required_non_negative_row_count(
            self.final_rows,
            field_name="row_attrition_report.final_rows",
        )
        if final_rows > input_rows:
            raise PhosPyInputError(
                "row_attrition_report.final_rows must be less than or equal "
                "to input_rows"
            )
        if records:
            if records[0].input_rows != input_rows:
                raise PhosPyInputError(
                    "row_attrition_report.input_rows must match the first "
                    "record input_rows"
                )
            if records[-1].output_rows != final_rows:
                raise PhosPyInputError(
                    "row_attrition_report.final_rows must match the last "
                    "record output_rows"
                )
            for previous, current in zip(records[:-1], records[1:], strict=True):
                if current.input_rows != previous.output_rows:
                    raise PhosPyInputError(
                        "row_attrition_report records must form a continuous "
                        "row-count chain"
                    )
        elif final_rows != input_rows:
            raise PhosPyInputError(
                "row_attrition_report.final_rows must equal input_rows when "
                "records is empty"
            )
        object.__setattr__(self, "records", records)
        object.__setattr__(self, "input_rows", input_rows)
        object.__setattr__(self, "final_rows", final_rows)

    @classmethod
    def from_records(
        cls,
        records: Sequence[RowAttritionRecord],
    ) -> RowAttritionReport:
        """Create a report using the first input and last output row counts."""

        record_tuple = _required_row_attrition_record_tuple(records)
        if not record_tuple:
            raise PhosPyInputError(
                "row_attrition_report.records must contain at least one record"
            )
        return cls(
            records=record_tuple,
            input_rows=record_tuple[0].input_rows,
            final_rows=record_tuple[-1].output_rows,
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible row-attrition report payload."""

        return {
            "records": [record.to_payload() for record in self.records],
            "input_rows": int(self.input_rows),
            "final_rows": int(self.final_rows),
        }


@dataclass(frozen=True, slots=True)
class TableFingerprint:
    """Deterministic table fingerprint metadata."""

    name: str
    rows: int
    columns: int
    index_name: str | None
    column_names: tuple[str, ...]
    dtypes: tuple[str, ...]
    exact_hash_algorithm: str
    exact_hash_value: str
    tolerance_hash_algorithm: str
    tolerance_hash_value: str
    index_structure: Mapping[str, JsonValue] | None = None
    column_index_structure: Mapping[str, JsonValue] | None = None

    def __post_init__(self) -> None:
        name = _required_provenance_text(self.name, field_name="table_fingerprint.name")
        rows = _required_non_negative_row_count(
            self.rows,
            field_name="table_fingerprint.rows",
        )
        columns = _required_non_negative_row_count(
            self.columns,
            field_name="table_fingerprint.columns",
        )
        column_names = _provenance_string_tuple(
            self.column_names,
            field_name="table_fingerprint.column_names",
        )
        dtypes = _provenance_string_tuple(
            self.dtypes,
            field_name="table_fingerprint.dtypes",
        )
        if len(column_names) != columns:
            raise PhosPyInputError(
                "table_fingerprint.column_names length must match "
                "table_fingerprint.columns"
            )
        if len(dtypes) != columns:
            raise PhosPyInputError(
                "table_fingerprint.dtypes length must match table_fingerprint.columns"
            )
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "rows", rows)
        object.__setattr__(self, "columns", columns)
        object.__setattr__(
            self,
            "index_name",
            None if self.index_name is None else str(self.index_name),
        )
        object.__setattr__(self, "column_names", column_names)
        object.__setattr__(self, "dtypes", dtypes)
        object.__setattr__(
            self,
            "exact_hash_algorithm",
            _required_provenance_text(
                self.exact_hash_algorithm,
                field_name="table_fingerprint.exact_hash_algorithm",
            ),
        )
        object.__setattr__(
            self,
            "exact_hash_value",
            _required_provenance_text(
                self.exact_hash_value,
                field_name="table_fingerprint.exact_hash_value",
            ),
        )
        object.__setattr__(
            self,
            "tolerance_hash_algorithm",
            _required_provenance_text(
                self.tolerance_hash_algorithm,
                field_name="table_fingerprint.tolerance_hash_algorithm",
            ),
        )
        object.__setattr__(
            self,
            "tolerance_hash_value",
            _required_provenance_text(
                self.tolerance_hash_value,
                field_name="table_fingerprint.tolerance_hash_value",
            ),
        )
        object.__setattr__(
            self,
            "index_structure",
            freeze_optional_json_mapping(
                self.index_structure,
                field_name="table_fingerprint.index_structure",
            ),
        )
        object.__setattr__(
            self,
            "column_index_structure",
            freeze_optional_json_mapping(
                self.column_index_structure,
                field_name="table_fingerprint.column_index_structure",
            ),
        )


def _required_row_attrition_record_tuple(
    records: object,
) -> tuple[RowAttritionRecord, ...]:
    if isinstance(records, (str, bytes, bytearray)) or not isinstance(
        records,
        Sequence,
    ):
        raise PhosPyInputError(
            "row_attrition_report.records must be a sequence of "
            "RowAttritionRecord values"
        )
    record_tuple = tuple(records)
    for record in record_tuple:
        if not isinstance(record, RowAttritionRecord):
            raise PhosPyInputError(
                "row_attrition_report.records must contain only "
                "RowAttritionRecord values"
            )
    return record_tuple


def _required_table_fingerprint_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    fingerprints = tuple(values)
    for fingerprint in fingerprints:
        if not isinstance(fingerprint, TableFingerprint):
            raise PhosPyInputError(
                f"{field_name} must contain only TableFingerprint values"
            )
    return fingerprints
