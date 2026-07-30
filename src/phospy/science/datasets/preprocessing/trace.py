"""Preprocessing state and row-audit trace models."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from enum import Enum

import pandas as pd

from phospy.science.datasets.preprocessing.batch_correction import BatchCorrectionReport
from phospy.science.datasets.preprocessing.batch_correction_metadata import (
    ResolvedBatchCorrectionMetadata,
)
from phospy.science.datasets.preprocessing.plan import PreprocessingPlan
from phospy.science.datasets.preprocessing.report_schema import (
    ROW_AUDIT_COLUMNS,
    PreprocessingRowAuditRow,
    dataframe_from_row_audit_rows,
    reorder_columns,
)
from phospy.science.datasets.preprocessing.results import PreprocessingReportRow


class PreprocessingStateTableKey(str, Enum):
    """Supported preprocessing state/report tables addressable in stage metadata."""

    DATASET_PHOSPHO = "dataset.phospho"
    DATASET_SITE_METADATA = "dataset.site_metadata"
    DATASET_SAMPLE_METADATA = "dataset.sample_metadata"
    DATASET_TOTAL = "dataset.total"
    DATASET_COMPARISONS = "dataset.comparisons"
    DATASET_IMPUTATION_OBSERVATION_MASK = "dataset.imputation_observation_mask"
    REPORT_COMPARISON_GROUP_STATS = "report.comparison_group_stats"
    REPORT_COMPARISON_PAIR_STATS = "report.comparison_pair_stats"
    REPORT_DUPLICATE_SITE_RESOLUTION = "report.duplicate_site_resolution"
    REPORT_METADATA_CONFLICTS = "report.metadata_conflicts"
    REPORT_ROW_AUDIT = "report.row_audit"


PREPROCESSING_STATE_TABLE_KEYS: tuple[PreprocessingStateTableKey, ...] = tuple(
    PreprocessingStateTableKey
)


@dataclass(frozen=True, slots=True)
class PreprocessingState:
    """Internal preprocessing state carried between ordered stages."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    plan: PreprocessingPlan
    comparisons: pd.DataFrame | None = None
    imputation_observation_mask: pd.DataFrame | None = None
    comparison_group_stats: pd.DataFrame | None = None
    comparison_pair_stats: pd.DataFrame | None = None
    duplicate_site_resolution: pd.DataFrame | None = None
    metadata_conflicts: pd.DataFrame | None = None
    row_audit: pd.DataFrame | None = None
    batch_correction_metadata: ResolvedBatchCorrectionMetadata | None = None
    batch_correction_report: BatchCorrectionReport | None = None
    report_rows: tuple[PreprocessingReportRow, ...] = ()


def empty_preprocessing_row_audit() -> pd.DataFrame:
    """Return an empty stable-schema preprocessing row-audit table."""

    return dataframe_from_row_audit_rows(())


def append_row_audit_records(
    state: PreprocessingState,
    records: Sequence[PreprocessingRowAuditRow],
) -> PreprocessingState:
    """Append row-audit records without mutating existing state frames."""

    if not records:
        return state

    existing = (
        empty_preprocessing_row_audit()
        if state.row_audit is None
        else reorder_columns(
            state.row_audit,
            expected_columns=ROW_AUDIT_COLUMNS,
        ).copy(deep=True)
    )
    appended = dataframe_from_row_audit_rows(records)
    combined = pd.concat([existing, appended], axis=0, ignore_index=True)
    return replace(state, row_audit=combined)
