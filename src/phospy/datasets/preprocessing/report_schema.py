"""Central schema ownership for dataset preprocessing report tables."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass, fields
from typing import ClassVar, Protocol, TypeVar, cast

import pandas as pd

PREPROCESSING_REPORT_ROW_COUNTS_TABLE = "row_counts"
PREPROCESSING_REPORT_OPERATIONS_TABLE = "operations"
PREPROCESSING_REPORT_ROW_AUDIT_TABLE = "row_audit"
PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE = "duplicate_site_resolution"
PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE = "metadata_conflicts"
PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE = "comparison_group_stats"
PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE = "comparison_pair_stats"


class _DataclassInstance(Protocol):
    __dataclass_fields__: ClassVar[dict[str, object]]


def _columns_for_dataclass(row_type: type[_DataclassInstance]) -> tuple[str, ...]:
    return tuple(
        field.name
        for field in fields(
            cast(type[object], row_type)  # pyright: ignore[reportArgumentType] - validated by dataclass decorators on row types
        )
    )


@dataclass(frozen=True, slots=True)
class PreprocessingRowCountRow:
    stage: str
    input_rows: int
    output_rows: int
    dropped_rows: int

    @classmethod
    def columns(cls) -> tuple[str, ...]:
        return _columns_for_dataclass(cls)


@dataclass(frozen=True, slots=True)
class PreprocessingOperationRow:
    step_order: int
    stage: str
    operation: str
    parameters: object
    input_rows: int
    output_rows: int
    notes: object

    @classmethod
    def columns(cls) -> tuple[str, ...]:
        return _columns_for_dataclass(cls)


@dataclass(frozen=True, slots=True)
class PreprocessingRowAuditRow:
    stage: str
    action: str
    reason: str
    source_row_id: str
    site_id: str
    retained: bool
    retained_row_id: object
    source_rows: object
    retained_row: object
    parameter_snapshot: object

    @classmethod
    def columns(cls) -> tuple[str, ...]:
        return _columns_for_dataclass(cls)


@dataclass(frozen=True, slots=True)
class DuplicateSiteResolutionRow:
    site_id: str
    source_row_id: str
    retained: bool
    resolution_policy: str
    retained_reason: object
    dropped_reason: object
    observed_values: object
    mean_signal: object
    n_source_rows: int
    n_aggregated_rows: object
    source_protein_id: object
    source_gene_symbol: object
    source_site: object
    source_site_sequence: object
    metadata_conflict_detected: bool

    @classmethod
    def columns(cls) -> tuple[str, ...]:
        return _columns_for_dataclass(cls)


@dataclass(frozen=True, slots=True)
class MetadataConflictRow:
    site_id: str
    field: str
    values: object
    n_distinct_values: int
    source_row_ids: object

    @classmethod
    def columns(cls) -> tuple[str, ...]:
        return _columns_for_dataclass(cls)


@dataclass(frozen=True, slots=True)
class ComparisonGroupStatsRow:
    site_id: str
    group: str
    n: int
    mean: object
    sd: object
    sem: object
    median: object
    min: object
    max: object
    sample_ids: object

    @classmethod
    def columns(cls) -> tuple[str, ...]:
        return _columns_for_dataclass(cls)


@dataclass(frozen=True, slots=True)
class ComparisonPairStatsRow:
    site_id: str
    comparison: str
    left_group: str
    right_group: str
    left_n: int
    right_n: int
    left_mean: object
    right_mean: object
    left_sd: object
    right_sd: object
    left_sem: object
    right_sem: object
    effect_size: object
    left_median: object
    right_median: object
    left_min: object
    right_min: object
    left_max: object
    right_max: object

    @classmethod
    def columns(cls) -> tuple[str, ...]:
        return _columns_for_dataclass(cls)


ROW_COUNTS_COLUMNS = PreprocessingRowCountRow.columns()
OPERATIONS_COLUMNS = PreprocessingOperationRow.columns()
ROW_AUDIT_COLUMNS = PreprocessingRowAuditRow.columns()
DUPLICATE_SITE_RESOLUTION_COLUMNS = DuplicateSiteResolutionRow.columns()
METADATA_CONFLICT_COLUMNS = MetadataConflictRow.columns()
COMPARISON_GROUP_STATS_COLUMNS = ComparisonGroupStatsRow.columns()
COMPARISON_PAIR_STATS_COLUMNS = ComparisonPairStatsRow.columns()

PREPROCESSING_REPORT_TABLE_COLUMNS = {
    PREPROCESSING_REPORT_ROW_COUNTS_TABLE: ROW_COUNTS_COLUMNS,
    PREPROCESSING_REPORT_OPERATIONS_TABLE: OPERATIONS_COLUMNS,
    PREPROCESSING_REPORT_ROW_AUDIT_TABLE: ROW_AUDIT_COLUMNS,
    PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE: (
        DUPLICATE_SITE_RESOLUTION_COLUMNS
    ),
    PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE: METADATA_CONFLICT_COLUMNS,
    PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE: COMPARISON_GROUP_STATS_COLUMNS,
    PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE: COMPARISON_PAIR_STATS_COLUMNS,
}

_ReportRowT = TypeVar("_ReportRowT", bound=_DataclassInstance)


def missing_columns(
    frame: pd.DataFrame, *, expected_columns: tuple[str, ...]
) -> tuple[str, ...]:
    return tuple(column for column in expected_columns if column not in frame.columns)


def reorder_columns(
    frame: pd.DataFrame, *, expected_columns: tuple[str, ...]
) -> pd.DataFrame:
    ordered = list(expected_columns)
    ordered.extend(
        column for column in frame.columns.tolist() if column not in expected_columns
    )
    return frame.loc[:, ordered]


def dataframe_from_rows(
    rows: Sequence[_ReportRowT],
    *,
    row_type: type[_ReportRowT],
) -> pd.DataFrame:
    columns = tuple(
        field.name
        for field in fields(
            cast(type[object], row_type)  # pyright: ignore[reportArgumentType] - row_type is constrained to dataclass-backed report rows
        )
    )
    if not rows:
        return pd.DataFrame.from_records([], columns=columns)
    records: list[dict[str, object]] = []
    for row in rows:
        if not isinstance(row, row_type):
            raise TypeError(
                f"expected rows of type {row_type.__name__}; got {type(row).__name__}"
            )
        records.append(
            asdict(
                cast(_DataclassInstance, row)  # pyright: ignore[reportArgumentType] - runtime isinstance(row, row_type) guarantees dataclass instance
            )
        )
    return pd.DataFrame.from_records(records, columns=columns)


def rows_from_dataframe(
    frame: pd.DataFrame | None,
    *,
    row_type: type[_ReportRowT],
) -> tuple[_ReportRowT, ...]:
    if frame is None:
        return ()
    columns = tuple(
        field.name
        for field in fields(
            cast(type[object], row_type)  # pyright: ignore[reportArgumentType] - row_type is constrained to dataclass-backed report rows
        )
    )
    missing = missing_columns(frame, expected_columns=columns)
    if missing:
        joined = ", ".join(missing)
        raise ValueError(
            f"cannot parse {row_type.__name__} rows; missing columns: {joined}"
        )
    if frame.empty:
        return ()
    records = frame.loc[:, list(columns)].to_dict(orient="records")
    return tuple(row_type(**record) for record in records)


def dataframe_from_row_count_rows(
    rows: Sequence[PreprocessingRowCountRow],
) -> pd.DataFrame:
    return dataframe_from_rows(rows, row_type=PreprocessingRowCountRow)


def dataframe_from_operation_rows(
    rows: Sequence[PreprocessingOperationRow],
) -> pd.DataFrame:
    return dataframe_from_rows(rows, row_type=PreprocessingOperationRow)


def dataframe_from_row_audit_rows(
    rows: Sequence[PreprocessingRowAuditRow],
) -> pd.DataFrame:
    return dataframe_from_rows(rows, row_type=PreprocessingRowAuditRow)


def dataframe_from_duplicate_site_resolution_rows(
    rows: Sequence[DuplicateSiteResolutionRow],
) -> pd.DataFrame:
    return dataframe_from_rows(rows, row_type=DuplicateSiteResolutionRow)


def dataframe_from_metadata_conflict_rows(
    rows: Sequence[MetadataConflictRow],
) -> pd.DataFrame:
    return dataframe_from_rows(rows, row_type=MetadataConflictRow)


def dataframe_from_comparison_group_stats_rows(
    rows: Sequence[ComparisonGroupStatsRow],
) -> pd.DataFrame:
    return dataframe_from_rows(rows, row_type=ComparisonGroupStatsRow)


def dataframe_from_comparison_pair_stats_rows(
    rows: Sequence[ComparisonPairStatsRow],
) -> pd.DataFrame:
    return dataframe_from_rows(rows, row_type=ComparisonPairStatsRow)


def row_count_rows_from_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[PreprocessingRowCountRow, ...]:
    return rows_from_dataframe(frame, row_type=PreprocessingRowCountRow)


def operation_rows_from_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[PreprocessingOperationRow, ...]:
    return rows_from_dataframe(frame, row_type=PreprocessingOperationRow)


def row_audit_rows_from_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[PreprocessingRowAuditRow, ...]:
    return rows_from_dataframe(frame, row_type=PreprocessingRowAuditRow)


def duplicate_site_resolution_rows_from_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[DuplicateSiteResolutionRow, ...]:
    return rows_from_dataframe(frame, row_type=DuplicateSiteResolutionRow)


def metadata_conflict_rows_from_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[MetadataConflictRow, ...]:
    return rows_from_dataframe(frame, row_type=MetadataConflictRow)


def comparison_group_stats_rows_from_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[ComparisonGroupStatsRow, ...]:
    return rows_from_dataframe(frame, row_type=ComparisonGroupStatsRow)


def comparison_pair_stats_rows_from_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[ComparisonPairStatsRow, ...]:
    return rows_from_dataframe(frame, row_type=ComparisonPairStatsRow)


__all__ = [
    "COMPARISON_GROUP_STATS_COLUMNS",
    "COMPARISON_PAIR_STATS_COLUMNS",
    "ComparisonGroupStatsRow",
    "ComparisonPairStatsRow",
    "DUPLICATE_SITE_RESOLUTION_COLUMNS",
    "DuplicateSiteResolutionRow",
    "METADATA_CONFLICT_COLUMNS",
    "MetadataConflictRow",
    "OPERATIONS_COLUMNS",
    "PreprocessingOperationRow",
    "PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE",
    "PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE",
    "PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE",
    "PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE",
    "PREPROCESSING_REPORT_OPERATIONS_TABLE",
    "PREPROCESSING_REPORT_ROW_AUDIT_TABLE",
    "PREPROCESSING_REPORT_ROW_COUNTS_TABLE",
    "PREPROCESSING_REPORT_TABLE_COLUMNS",
    "PreprocessingRowAuditRow",
    "PreprocessingRowCountRow",
    "ROW_AUDIT_COLUMNS",
    "ROW_COUNTS_COLUMNS",
    "comparison_group_stats_rows_from_dataframe",
    "comparison_pair_stats_rows_from_dataframe",
    "dataframe_from_comparison_group_stats_rows",
    "dataframe_from_comparison_pair_stats_rows",
    "dataframe_from_duplicate_site_resolution_rows",
    "dataframe_from_metadata_conflict_rows",
    "dataframe_from_operation_rows",
    "dataframe_from_row_audit_rows",
    "dataframe_from_row_count_rows",
    "dataframe_from_rows",
    "duplicate_site_resolution_rows_from_dataframe",
    "metadata_conflict_rows_from_dataframe",
    "missing_columns",
    "operation_rows_from_dataframe",
    "reorder_columns",
    "row_audit_rows_from_dataframe",
    "row_count_rows_from_dataframe",
    "rows_from_dataframe",
]
