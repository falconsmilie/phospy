"""Preprocessing-owned stage-only report schema and table registry."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import cast

import pandas as pd

import phospy.preparation_records.report_schema as _durable_schema

PREPROCESSING_REPORT_ROW_COUNTS_TABLE = "row_counts"
PREPROCESSING_REPORT_OPERATIONS_TABLE = "operations"
PREPROCESSING_REPORT_ROW_AUDIT_TABLE = "row_audit"
PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE = "duplicate_site_resolution"
PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE = "metadata_conflicts"
PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE = "comparison_group_stats"
PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE = "comparison_pair_stats"
PREPROCESSING_REPORT_SITE_SEQUENCE_RESOLUTION_TABLE = "site_sequence_resolution"


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionRow:
    row_id: str
    status: str
    reason: object
    accession: object
    site_token: object
    existing_site_sequence: object
    resolved_site_sequence: object
    action: str

    @classmethod
    def columns(cls) -> tuple[str, ...]:
        return tuple(
            field.name
            for field in fields(
                cast(type[object], cls)  # pyright: ignore[reportArgumentType]
            )
        )


SITE_SEQUENCE_RESOLUTION_COLUMNS = SiteSequenceResolutionRow.columns()

PREPROCESSING_REPORT_TABLE_COLUMNS = {
    PREPROCESSING_REPORT_ROW_COUNTS_TABLE: _durable_schema.ROW_COUNTS_COLUMNS,
    PREPROCESSING_REPORT_OPERATIONS_TABLE: _durable_schema.OPERATIONS_COLUMNS,
    PREPROCESSING_REPORT_ROW_AUDIT_TABLE: _durable_schema.ROW_AUDIT_COLUMNS,
    PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE: (
        _durable_schema.DUPLICATE_SITE_RESOLUTION_COLUMNS
    ),
    PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE: (
        _durable_schema.METADATA_CONFLICT_COLUMNS
    ),
    PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE: (
        _durable_schema.COMPARISON_GROUP_STATS_COLUMNS
    ),
    PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE: (
        _durable_schema.COMPARISON_PAIR_STATS_COLUMNS
    ),
    PREPROCESSING_REPORT_SITE_SEQUENCE_RESOLUTION_TABLE: (
        SITE_SEQUENCE_RESOLUTION_COLUMNS
    ),
}


def dataframe_from_site_sequence_resolution_rows(
    rows: tuple[SiteSequenceResolutionRow, ...],
) -> pd.DataFrame:
    return _durable_schema.dataframe_from_rows(
        rows,
        row_type=SiteSequenceResolutionRow,
    )


def site_sequence_resolution_rows_from_dataframe(
    frame: pd.DataFrame | None,
) -> tuple[SiteSequenceResolutionRow, ...]:
    return _durable_schema.rows_from_dataframe(
        frame,
        row_type=SiteSequenceResolutionRow,
    )


__all__ = [
    "PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_TABLE",
    "PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_TABLE",
    "PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_TABLE",
    "PREPROCESSING_REPORT_METADATA_CONFLICTS_TABLE",
    "PREPROCESSING_REPORT_OPERATIONS_TABLE",
    "PREPROCESSING_REPORT_ROW_AUDIT_TABLE",
    "PREPROCESSING_REPORT_ROW_COUNTS_TABLE",
    "PREPROCESSING_REPORT_SITE_SEQUENCE_RESOLUTION_TABLE",
    "PREPROCESSING_REPORT_TABLE_COLUMNS",
    "SITE_SEQUENCE_RESOLUTION_COLUMNS",
    "SiteSequenceResolutionRow",
    "dataframe_from_site_sequence_resolution_rows",
    "site_sequence_resolution_rows_from_dataframe",
]
