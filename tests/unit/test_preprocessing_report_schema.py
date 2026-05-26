from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    DatasetComparisonBuildingConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.errors.build import DatasetBuildError
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.models import PreprocessingReportRow
from phospy.science.datasets.preprocessing.report_rows import (
    validate_preprocessing_report_row,
)
from phospy.science.datasets.preprocessing.report_schema import (
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
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "localisation_confidence": [0.95, 0.9, 0.92],
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
            input_intensity_scale="linear",
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
    assert int(site_matrix_counts.iloc[0]["output_rows"]) == 3
    assert int(site_matrix_counts.iloc[0]["dropped_rows"]) == 0
    assert report.duplicate_site_resolution.empty

    pair_stats = report.comparison_pair_stats.loc[
        (report.comparison_pair_stats.loc[:, "site_id"] == "MAPK14;Y182;")
        & (report.comparison_pair_stats.loc[:, "comparison"] == "p_g1_g2")
    ]
    assert pair_stats.shape[0] == 1
    assert float(pair_stats.iloc[0]["effect_size"]) == 6.0


def test_report_schema_stable_with_missing_data_diagnostics_enabled() -> None:
    phospho = pd.DataFrame(
        {
            "sample_1": [10.0, 8.0],
            "sample_2": [12.0, float("nan")],
        },
        index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.9],
        },
        index=phospho.index.copy(),
    )

    built = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            input_intensity_scale="linear",
            preprocessing_config=DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                ),
            ),
        )
    )

    assert built.preprocessing_report is not None
    assert tuple(built.preprocessing_report.operations.columns) == OPERATIONS_COLUMNS
    assert built.processing_state.missing_data.diagnostics is not None
    diagnostics_payload = built.processing_state.missing_data.diagnostics.to_payload()
    assert isinstance(diagnostics_payload["missingness_mask_hash"], str)


def test_row_audit_parameter_snapshot_must_be_mapping() -> None:
    with pytest.raises(
        DatasetBuildError,
        match="row_audit.parameter_snapshot with invalid type",
    ):
        validate_preprocessing_report_row(
            PreprocessingReportRow(
                table="row_audit",
                values=PreprocessingRowAuditRow(
                    stage="missing_data",
                    action="imputed",
                    reason="invalid test payload",
                    source_row_id="row_a",
                    site_id="row_a",
                    retained=True,
                    retained_row_id="row_a",
                    source_rows=("row_a",),
                    retained_row="row_a",
                    parameter_snapshot="not-a-mapping",
                ),
            )
        )


def test_docs_state_row_median_missing_data_semantics() -> None:
    root = Path(__file__).resolve().parents[2]
    workflow_contracts = (root / "docs" / "workflow_contracts.md").read_text(
        encoding="utf-8"
    )
    validation = (root / "docs" / "validation.md").read_text(encoding="utf-8")

    combined = f"{workflow_contracts}\n{validation}"
    assert "missing-data handling runs before normalisation" in combined
    assert "row-median imputation is deterministic" in combined
    assert "row-median imputation is not left-censored imputation" in combined
