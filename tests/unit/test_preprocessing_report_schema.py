from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.datasets.models import DatasetPreprocessingReport
from phospy.datasets.preprocessing.report_schema import (
    COMPARISON_GROUP_STATS_COLUMNS,
    COMPARISON_PAIR_STATS_COLUMNS,
    DUPLICATE_SITE_RESOLUTION_COLUMNS,
    METADATA_CONFLICT_COLUMNS,
    OPERATIONS_COLUMNS,
    ROW_AUDIT_COLUMNS,
    ROW_COUNTS_COLUMNS,
    PreprocessingOperationRow,
    PreprocessingRowAuditRow,
    PreprocessingRowCountRow,
    dataframe_from_comparison_group_stats_rows,
    dataframe_from_comparison_pair_stats_rows,
    dataframe_from_duplicate_site_resolution_rows,
    dataframe_from_metadata_conflict_rows,
    dataframe_from_operation_rows,
    dataframe_from_row_audit_rows,
    dataframe_from_row_count_rows,
)


def test_row_count_rows_convert_to_dataframe_with_schema_order() -> None:
    frame = dataframe_from_row_count_rows(
        (
            PreprocessingRowCountRow(
                stage="missing_data",
                input_rows=3,
                output_rows=2,
                dropped_rows=1,
            ),
        )
    )

    assert tuple(frame.columns) == ROW_COUNTS_COLUMNS
    assert frame.loc[0, "stage"] == "missing_data"
    assert int(frame.loc[0, "dropped_rows"]) == 1


def test_empty_row_count_table_has_schema_columns() -> None:
    frame = dataframe_from_row_count_rows(())

    assert frame.empty
    assert tuple(frame.columns) == ROW_COUNTS_COLUMNS


@pytest.mark.parametrize(
    ("factory", "expected_columns"),
    (
        (dataframe_from_row_count_rows, ROW_COUNTS_COLUMNS),
        (dataframe_from_operation_rows, OPERATIONS_COLUMNS),
        (dataframe_from_row_audit_rows, ROW_AUDIT_COLUMNS),
        (
            dataframe_from_duplicate_site_resolution_rows,
            DUPLICATE_SITE_RESOLUTION_COLUMNS,
        ),
        (dataframe_from_metadata_conflict_rows, METADATA_CONFLICT_COLUMNS),
        (dataframe_from_comparison_group_stats_rows, COMPARISON_GROUP_STATS_COLUMNS),
        (dataframe_from_comparison_pair_stats_rows, COMPARISON_PAIR_STATS_COLUMNS),
    ),
)
def test_empty_tables_preserve_schema_columns(
    factory: object,
    expected_columns: tuple[str, ...],
) -> None:
    frame = factory(())

    assert frame.empty
    assert tuple(frame.columns) == expected_columns


def test_preprocessing_report_from_rows_uses_central_schemas() -> None:
    report = DatasetPreprocessingReport.from_rows(
        row_count_rows=(
            PreprocessingRowCountRow(
                stage="missing_data",
                input_rows=2,
                output_rows=2,
                dropped_rows=0,
            ),
        ),
        operation_rows=(
            PreprocessingOperationRow(
                step_order=1,
                stage="missing_data",
                operation="forbid",
                parameters={"missing_data_policy": "forbid"},
                input_rows=2,
                output_rows=2,
                notes="stage executed",
            ),
        ),
        row_audit_rows=(
            PreprocessingRowAuditRow(
                stage="missing_data",
                action="retained",
                reason="stage passthrough",
                source_row_id="row_a",
                site_id="row_a",
                retained=True,
                retained_row_id="row_a",
                source_rows=("row_a",),
                retained_row="row_a",
                parameter_snapshot={"missing_data_policy": "forbid"},
            ),
        ),
        duplicate_site_resolution_rows=(),
        metadata_conflict_rows=(),
        comparison_group_stats_rows=(),
        comparison_pair_stats_rows=(),
    )

    assert tuple(report.row_counts.columns) == ROW_COUNTS_COLUMNS
    assert tuple(report.operations.columns) == OPERATIONS_COLUMNS
    assert tuple(report.row_audit.columns) == ROW_AUDIT_COLUMNS
    assert tuple(report.duplicate_site_resolution.columns) == (
        DUPLICATE_SITE_RESOLUTION_COLUMNS
    )
    assert tuple(report.metadata_conflicts.columns) == METADATA_CONFLICT_COLUMNS
    assert (
        tuple(report.comparison_group_stats.columns) == COMPARISON_GROUP_STATS_COLUMNS
    )
    assert tuple(report.comparison_pair_stats.columns) == COMPARISON_PAIR_STATS_COLUMNS


def test_preprocessing_report_integration_preserves_representative_values() -> None:
    phospho = pd.DataFrame(
        {
            "sample_1": [10.0, 8.0, 5.0],
            "sample_2": [12.0, 10.0, 5.0],
            "sample_3": [4.0, 2.0, 5.0],
            "sample_4": [6.0, 4.0, 5.0],
        },
        index=pd.Index(["row_a", "row_b", "row_c"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14", "AKT1"],
            "site": ["Y182", "Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {"comparison_group": ["g1", "g1", "g2", "g2"]},
        index=phospho.columns.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            preprocessing_config=DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_policy="first",
                ),
                comparisons=DatasetComparisonBuildingConfig(
                    policy="sample_metadata_pairs"
                ),
            ),
        )
    )

    assert built.preprocessing_report is not None
    report = built.preprocessing_report

    site_matrix_counts = report.row_counts.loc[
        report.row_counts.loc[:, "stage"] == "site_matrix"
    ]
    assert site_matrix_counts.shape[0] == 1
    assert int(site_matrix_counts.iloc[0]["input_rows"]) == 3
    assert int(site_matrix_counts.iloc[0]["output_rows"]) == 2
    assert int(site_matrix_counts.iloc[0]["dropped_rows"]) == 1

    collapsed = report.row_audit.loc[
        (report.row_audit.loc[:, "stage"] == "site_matrix")
        & (report.row_audit.loc[:, "action"] == "collapsed")
    ]
    assert collapsed.shape[0] == 1
    assert collapsed.iloc[0]["source_row_id"] == "row_b"

    duplicate_rows = report.duplicate_site_resolution
    assert not duplicate_rows.empty
    dropped_duplicate = duplicate_rows.loc[
        duplicate_rows.loc[:, "source_row_id"] == "row_b"
    ]
    assert dropped_duplicate.shape[0] == 1
    assert bool(dropped_duplicate.iloc[0]["retained"]) is False

    pair_stats = report.comparison_pair_stats.loc[
        (report.comparison_pair_stats.loc[:, "site_id"] == "MAPK14;Y182;")
        & (report.comparison_pair_stats.loc[:, "comparison"] == "p_g1_g2")
    ]
    assert pair_stats.shape[0] == 1
    assert float(pair_stats.iloc[0]["effect_size"]) == 6.0
