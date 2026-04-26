from __future__ import annotations

import pandas as pd
import pytest

from phospy.api.configs import (
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
)
from phospy.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.datasets.models import DatasetPreprocessingReport
from phospy.datasets.preprocessing.models import PreprocessingPlan
from phospy.errors.validation import DatasetValidationError

_ROW_AUDIT_COLUMNS = (
    "stage",
    "action",
    "reason",
    "source_row_id",
    "site_id",
    "retained",
    "retained_row_id",
    "source_rows",
    "retained_row",
    "parameter_snapshot",
)


def test_missing_data_stage_audits_rows_dropped_below_min_observed_values() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan"), float("nan")],
            "sample_b": [2.0, 2.0, float("nan")],
            "sample_c": [3.0, 4.0, 7.0],
        },
        index=pd.Index(["row_keep", "row_impute", "row_drop"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C"],
            "site": ["S1", "S2", "S3"],
            "site_sequence": ["SEQ_A", "SEQ_B", "SEQ_C"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=2,
                )
            )
        ),
    )

    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit["stage"] == "missing_data")
        & (preprocessed.row_audit["action"] == "dropped")
    ]
    assert dropped.shape[0] == 1
    assert dropped.iloc[0]["source_row_id"] == "row_drop"
    assert dropped.iloc[0]["site_id"] == "row_drop"
    assert dropped.iloc[0]["retained"] is False
    assert "below missing_data.min_observed_values" in str(dropped.iloc[0]["reason"])


def test_missing_data_stage_audits_row_median_imputation() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, float("nan")],
            "sample_b": [2.0, 2.0],
            "sample_c": [3.0, float("nan")],
        },
        index=pd.Index(["row_keep", "row_impute"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B"],
            "site": ["S1", "S2"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                missing_data=DatasetMissingDataConfig(
                    policy="impute_row_median",
                    min_observed_values=1,
                )
            )
        ),
    )

    imputed = preprocessed.row_audit.loc[
        (preprocessed.row_audit["stage"] == "missing_data")
        & (preprocessed.row_audit["action"] == "imputed")
    ]
    assert imputed.shape[0] == 1
    assert imputed.iloc[0]["source_row_id"] == "row_impute"
    snapshot = imputed.iloc[0]["parameter_snapshot"]
    assert isinstance(snapshot, dict)
    assert "imputed_columns" in snapshot
    assert set(snapshot["imputed_columns"]) == {"sample_a", "sample_c"}


def test_site_matrix_stage_audits_missing_sequence_drops() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            )
        ),
    )

    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit["stage"] == "site_matrix")
        & (preprocessed.row_audit["action"] == "dropped")
    ]
    assert dropped.shape[0] == 1
    assert dropped.iloc[0]["source_row_id"] == "SRC_ROW_1"
    assert dropped.iloc[0]["site_id"] == "MAPK14;Y182;"
    assert "site_sequence is missing or blank" in str(dropped.iloc[0]["reason"])


def test_site_matrix_stage_audits_incomplete_value_drops() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, float("nan")],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
            )
        ),
    )

    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit["stage"] == "site_matrix")
        & (preprocessed.row_audit["action"] == "dropped")
    ]
    assert dropped.shape[0] == 1
    assert dropped.iloc[0]["source_row_id"] == "SRC_ROW_2"
    assert dropped.iloc[0]["site_id"] == "AKT1;T308;"
    assert dropped.iloc[0]["reason"] == "dropped by site_matrix missing-data policy"


def test_site_matrix_stage_audits_duplicate_resolution_first_strategy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_strategy="first",
                )
            )
        ),
    )

    site_matrix_audit = preprocessed.row_audit.loc[
        preprocessed.row_audit["stage"] == "site_matrix"
    ].set_index("source_row_id")
    assert site_matrix_audit.loc["SRC_ROW_1", "action"] == "retained"
    assert site_matrix_audit.loc["SRC_ROW_1", "retained"] is True
    assert site_matrix_audit.loc["SRC_ROW_2", "action"] == "collapsed"
    assert site_matrix_audit.loc["SRC_ROW_2", "retained"] is False
    assert site_matrix_audit.loc["SRC_ROW_2", "retained_row_id"] == "SRC_ROW_1"


def test_site_matrix_stage_audits_duplicate_resolution_max_mean_signal_strategy() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 5.0],
            "sample_b": [2.0, 6.0],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_strategy="max_mean_signal",
                )
            )
        ),
    )

    site_matrix_audit = preprocessed.row_audit.loc[
        preprocessed.row_audit["stage"] == "site_matrix"
    ].set_index("source_row_id")
    assert site_matrix_audit.loc["SRC_ROW_2", "action"] == "retained"
    assert site_matrix_audit.loc["SRC_ROW_2", "retained"] is True
    assert site_matrix_audit.loc["SRC_ROW_1", "action"] == "collapsed"
    assert site_matrix_audit.loc["SRC_ROW_1", "retained"] is False
    assert site_matrix_audit.loc["SRC_ROW_1", "retained_row_id"] == "SRC_ROW_2"


@pytest.mark.parametrize(
    "duplicate_strategy",
    ("aggregate_mean", "aggregate_median"),
)
def test_site_matrix_stage_audits_aggregate_duplicate_contributors(
    duplicate_strategy: str,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "MAPK14"],
            "site": ["Y182", "Y182"],
            "site_sequence": ["SEQ_A", "SEQ_B"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.from_config(
            DatasetPreprocessingConfig(
                site_matrix=DatasetSiteMatrixConfig(
                    policy="build_from_metadata",
                    duplicate_site_strategy=duplicate_strategy,
                )
            )
        ),
    )

    aggregated = preprocessed.row_audit.loc[
        (preprocessed.row_audit["stage"] == "site_matrix")
        & (preprocessed.row_audit["action"] == "aggregated")
    ]
    assert aggregated.shape[0] == 2
    assert set(aggregated["source_row_id"]) == {"SRC_ROW_1", "SRC_ROW_2"}
    assert set(aggregated["retained_row_id"]) == {"MAPK14;Y182;"}
    assert {
        tuple(source_rows) for source_rows in aggregated["source_rows"].tolist()
    } == {("SRC_ROW_1", "SRC_ROW_2")}


def test_empty_row_audit_has_stable_schema() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["SEQ_A"],
        },
        index=phospho.index.copy(),
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=PreprocessingPlan.default(),
    )

    assert preprocessed.row_audit is not None
    assert tuple(preprocessed.row_audit.columns) == _ROW_AUDIT_COLUMNS
    assert preprocessed.row_audit.empty


def test_preprocessing_report_validation_requires_row_audit_columns() -> None:
    row_counts = pd.DataFrame.from_records(
        [
            {
                "stage": "missing_data",
                "input_rows": 1,
                "output_rows": 1,
                "dropped_rows": 0,
            }
        ]
    )
    operations = pd.DataFrame.from_records(
        [
            {
                "step_order": 1,
                "stage": "missing_data",
                "operation": "forbid",
                "parameters": {},
                "input_rows": 1,
                "output_rows": 1,
                "notes": "stage executed",
            }
        ]
    )
    malformed_row_audit = pd.DataFrame.from_records(
        [{"stage": "missing_data", "action": "retained"}]
    )

    with pytest.raises(
        DatasetValidationError,
        match="dataset.preprocessing_report.row_audit is missing required columns",
    ):
        DatasetPreprocessingReport(
            row_counts=row_counts,
            operations=operations,
            row_audit=malformed_row_audit,
        )
