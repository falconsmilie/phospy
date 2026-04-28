"""Typed stage-owned preprocessing report-row composition helpers."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from phospy.datasets.preprocessing.models import (
    PreprocessingReportRow,
    StageOwnedPreprocessingReportValue,
)
from phospy.datasets.preprocessing.report_schema import (
    PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE,
    PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE,
    PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE,
    PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE,
    PREPROCESSING_REPORT_ROW_AUDIT_TABLE,
    ComparisonGroupStatsRow,
    ComparisonPairStatsRow,
    DuplicateSiteResolutionRow,
    MetadataConflictRow,
    PreprocessingRowAuditRow,
    comparison_group_stats_rows_from_dataframe,
    comparison_pair_stats_rows_from_dataframe,
    dataframe_from_comparison_group_stats_rows,
    dataframe_from_comparison_pair_stats_rows,
    dataframe_from_duplicate_site_resolution_rows,
    dataframe_from_metadata_conflict_rows,
    dataframe_from_row_audit_rows,
    duplicate_site_resolution_rows_from_dataframe,
    metadata_conflict_rows_from_dataframe,
    row_audit_rows_from_dataframe,
)
from phospy.errors.build import DatasetBuildError

_SUPPORTED_STAGE_REPORT_TABLES: dict[str, type[StageOwnedPreprocessingReportValue]] = {
    PREPROCESSING_REPORT_ROW_AUDIT_TABLE: PreprocessingRowAuditRow,
    PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE: DuplicateSiteResolutionRow,
    PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE: MetadataConflictRow,
    PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE: ComparisonGroupStatsRow,
    PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE: ComparisonPairStatsRow,
}


@dataclass(frozen=True, slots=True)
class StageOwnedPreprocessingReportTables:
    """Composed report tables contributed directly by preprocessing stages."""

    row_audit: pd.DataFrame
    duplicate_site_resolution: pd.DataFrame
    metadata_conflicts: pd.DataFrame
    comparison_group_stats: pd.DataFrame
    comparison_pair_stats: pd.DataFrame


def validate_preprocessing_report_row(
    row: PreprocessingReportRow,
) -> PreprocessingReportRow:
    """Validate that a stage-owned report row targets a supported typed table."""

    if not isinstance(row, PreprocessingReportRow):
        raise DatasetBuildError(
            "dataset preprocessing stage returned an invalid report row payload"
        )
    expected_type = _SUPPORTED_STAGE_REPORT_TABLES.get(str(row.table))
    if expected_type is None:
        raise DatasetBuildError(
            "dataset preprocessing stage emitted report rows for unsupported table: "
            f"{row.table!r}"
        )
    if not isinstance(row.values, expected_type):
        raise DatasetBuildError(
            "dataset preprocessing stage emitted report row values with invalid type "
            f"for table {row.table!r}: expected {expected_type.__name__}, got "
            f"{type(row.values).__name__}"
        )
    return row


def compose_stage_owned_report_tables(
    rows: Sequence[PreprocessingReportRow],
) -> StageOwnedPreprocessingReportTables:
    """Compose stage-owned report rows into the public report sidecar tables."""

    row_audit_rows: list[PreprocessingRowAuditRow] = []
    duplicate_site_resolution_rows: list[DuplicateSiteResolutionRow] = []
    metadata_conflict_rows: list[MetadataConflictRow] = []
    comparison_group_stats_rows: list[ComparisonGroupStatsRow] = []
    comparison_pair_stats_rows: list[ComparisonPairStatsRow] = []

    for raw_row in rows:
        row = validate_preprocessing_report_row(raw_row)
        if row.table == PREPROCESSING_REPORT_ROW_AUDIT_TABLE:
            if not isinstance(row.values, PreprocessingRowAuditRow):
                raise DatasetBuildError(
                    "invalid row-audit preprocessing report row type"
                )
            row_audit_rows.append(row.values)
            continue
        if row.table == PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE:
            if not isinstance(row.values, DuplicateSiteResolutionRow):
                raise DatasetBuildError(
                    "invalid duplicate-site-resolution preprocessing report row type"
                )
            duplicate_site_resolution_rows.append(row.values)
            continue
        if row.table == PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE:
            if not isinstance(row.values, MetadataConflictRow):
                raise DatasetBuildError(
                    "invalid metadata-conflict preprocessing report row type"
                )
            metadata_conflict_rows.append(row.values)
            continue
        if row.table == PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE:
            if not isinstance(row.values, ComparisonGroupStatsRow):
                raise DatasetBuildError(
                    "invalid comparison-group-stats preprocessing report row type"
                )
            comparison_group_stats_rows.append(row.values)
            continue
        if row.table == PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE:
            if not isinstance(row.values, ComparisonPairStatsRow):
                raise DatasetBuildError(
                    "invalid comparison-pair-stats preprocessing report row type"
                )
            comparison_pair_stats_rows.append(row.values)
            continue
        raise DatasetBuildError(
            "dataset preprocessing stage emitted report rows for unsupported table: "
            f"{row.table!r}"
        )

    return StageOwnedPreprocessingReportTables(
        row_audit=dataframe_from_row_audit_rows(tuple(row_audit_rows)),
        duplicate_site_resolution=dataframe_from_duplicate_site_resolution_rows(
            tuple(duplicate_site_resolution_rows)
        ),
        metadata_conflicts=dataframe_from_metadata_conflict_rows(
            tuple(metadata_conflict_rows)
        ),
        comparison_group_stats=dataframe_from_comparison_group_stats_rows(
            tuple(comparison_group_stats_rows)
        ),
        comparison_pair_stats=dataframe_from_comparison_pair_stats_rows(
            tuple(comparison_pair_stats_rows)
        ),
    )


def report_rows_from_row_audit_rows(
    rows: Sequence[PreprocessingRowAuditRow],
) -> tuple[PreprocessingReportRow, ...]:
    return tuple(
        PreprocessingReportRow(
            table=PREPROCESSING_REPORT_ROW_AUDIT_TABLE,
            values=row,
        )
        for row in rows
    )


def report_rows_from_duplicate_site_resolution_rows(
    rows: Sequence[DuplicateSiteResolutionRow],
) -> tuple[PreprocessingReportRow, ...]:
    return tuple(
        PreprocessingReportRow(
            table=PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE,
            values=row,
        )
        for row in rows
    )


def report_rows_from_metadata_conflict_rows(
    rows: Sequence[MetadataConflictRow],
) -> tuple[PreprocessingReportRow, ...]:
    return tuple(
        PreprocessingReportRow(
            table=PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE,
            values=row,
        )
        for row in rows
    )


def report_rows_from_comparison_group_stats_rows(
    rows: Sequence[ComparisonGroupStatsRow],
) -> tuple[PreprocessingReportRow, ...]:
    return tuple(
        PreprocessingReportRow(
            table=PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE,
            values=row,
        )
        for row in rows
    )


def report_rows_from_comparison_pair_stats_rows(
    rows: Sequence[ComparisonPairStatsRow],
) -> tuple[PreprocessingReportRow, ...]:
    return tuple(
        PreprocessingReportRow(
            table=PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE,
            values=row,
        )
        for row in rows
    )


def report_rows_from_row_audit_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[PreprocessingReportRow, ...]:
    return report_rows_from_row_audit_rows(row_audit_rows_from_dataframe(frame))


def report_rows_from_duplicate_site_resolution_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[PreprocessingReportRow, ...]:
    return report_rows_from_duplicate_site_resolution_rows(
        duplicate_site_resolution_rows_from_dataframe(frame)
    )


def report_rows_from_metadata_conflicts_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[PreprocessingReportRow, ...]:
    return report_rows_from_metadata_conflict_rows(
        metadata_conflict_rows_from_dataframe(frame)
    )


def report_rows_from_comparison_group_stats_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[PreprocessingReportRow, ...]:
    return report_rows_from_comparison_group_stats_rows(
        comparison_group_stats_rows_from_dataframe(frame)
    )


def report_rows_from_comparison_pair_stats_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[PreprocessingReportRow, ...]:
    return report_rows_from_comparison_pair_stats_rows(
        comparison_pair_stats_rows_from_dataframe(frame)
    )


__all__ = [
    "compose_stage_owned_report_tables",
    "report_rows_from_comparison_group_stats_dataframe",
    "report_rows_from_comparison_group_stats_rows",
    "report_rows_from_comparison_pair_stats_dataframe",
    "report_rows_from_comparison_pair_stats_rows",
    "report_rows_from_duplicate_site_resolution_dataframe",
    "report_rows_from_duplicate_site_resolution_rows",
    "report_rows_from_metadata_conflict_rows",
    "report_rows_from_metadata_conflicts_dataframe",
    "report_rows_from_row_audit_dataframe",
    "report_rows_from_row_audit_rows",
    "StageOwnedPreprocessingReportTables",
    "StageOwnedPreprocessingReportValue",
    "validate_preprocessing_report_row",
]
