from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from phospy.api.configs import (
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
)
from phospy.errors.validation import DatasetValidationError
from phospy.science.datasets.builders.preprocessing import DatasetPreprocessor
from phospy.science.datasets.models import DatasetPreprocessingReport
from phospy.science.datasets.preprocessing.models import PreprocessingPlan
from phospy.science.datasets.preprocessing.report_schema import ROW_AUDIT_COLUMNS
from tests.support.site_keys import site_key_index_from_display_ids


def _site_key(display_id: str) -> str:
    return str(
        site_key_index_from_display_ids(
            [display_id],
            protein_namespace="gene_symbol",
        )[0]
    )


def _with_site_identity(site_metadata: pd.DataFrame) -> pd.DataFrame:
    identified = site_metadata.copy()
    display_ids = [
        f"{gene_symbol};{site};"
        for gene_symbol, site in zip(
            identified.loc[:, "gene_symbol"].astype(str).tolist(),
            identified.loc[:, "site"].astype(str).tolist(),
            strict=True,
        )
    ]
    site_keys = site_key_index_from_display_ids(
        display_ids,
        protein_namespace="gene_symbol",
    )
    identified.insert(0, "display_id", display_ids)
    identified.insert(0, "site_key", site_keys.astype(str).tolist())
    return identified


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
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "localisation_confidence": [0.95, 0.9, 0.92],
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
                    input_scale="linear",
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
    assert bool(dropped.iloc[0]["retained"]) is False
    assert "below missing_data.min_observed_values" in str(dropped.iloc[0]["reason"])
    snapshot = dropped.iloc[0]["parameter_snapshot"]
    assert isinstance(snapshot, dict)
    assert snapshot["observed_values"] == 1
    assert snapshot["input_missing_cell_count"] == 3
    assert snapshot["output_missing_cell_count"] == 0
    assert snapshot["imputed_cell_count"] == 1
    assert snapshot["affected_row_count"] == 2
    assert snapshot["affected_column_count"] == 2
    assert isinstance(snapshot["missingness_mask_hash"], str)
    assert tuple(snapshot["stage_order"]) == (
        "localisation_confidence",
        "missing_data",
    )


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
            "site_sequence": ["SEQ_A", "SEQ_R"],
            "localisation_confidence": [0.95, 0.9],
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
                    input_scale="linear",
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
    assert snapshot["row_median"] == 2.0
    assert snapshot["imputed_cell_count"] == 2
    assert snapshot["input_missing_cell_count"] == 2
    assert snapshot["output_missing_cell_count"] == 0
    assert isinstance(snapshot["missingness_mask_hash"], str)
    assert tuple(snapshot["stage_order"]) == (
        "localisation_confidence",
        "missing_data",
    )


def test_missing_data_stage_audits_knn_column_mean_fallback_imputation() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [float("nan"), 2.0, 4.0],
            "sample_b": [float("nan"), 10.0, 20.0],
            "sample_c": [float("nan"), 30.0, 50.0],
        },
        index=pd.Index(["all_missing", "row_ref_1", "row_ref_2"], name="source_row"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["A", "B", "C"],
            "site": ["S1", "S2", "S3"],
            "site_sequence": ["SEQ_A", "SEQ_R", "SEQ_C"],
            "localisation_confidence": [0.95, 0.9, 0.92],
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
                    policy="impute_knn",
                    k=1,
                    distance="nan_euclidean",
                    max_missing_fraction_per_row=1.0,
                    input_scale="linear",
                )
            )
        ),
    )

    imputed = preprocessed.row_audit.loc[
        (preprocessed.row_audit["stage"] == "missing_data")
        & (preprocessed.row_audit["action"] == "imputed")
    ]
    assert imputed.shape[0] == 1
    assert imputed.iloc[0]["source_row_id"] == "all_missing"
    snapshot = imputed.iloc[0]["parameter_snapshot"]
    assert isinstance(snapshot, dict)
    assert snapshot["nearest_neighbour_imputed_columns"] == ()
    assert snapshot["column_mean_fallback_columns"] == (
        "sample_a",
        "sample_b",
        "sample_c",
    )
    assert snapshot["nearest_neighbour_imputed_cell_count"] == 0
    assert snapshot["column_mean_fallback_imputed_cell_count"] == 3
    assert snapshot["fully_column_mean_fallback_imputed"] is True
    assert snapshot["no_overlap_policy"] == "column_mean_with_caveat"
    assert snapshot["no_overlap_policy_version"] == 1


def test_site_matrix_stage_audits_missing_sequence_drops() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = _with_site_identity(
        pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["", "SEQ_R"],
                "localisation_confidence": [0.95, 0.9],
            },
            index=phospho.index.copy(),
        )
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=replace(
            PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                )
            ),
            stage_order=("site_matrix",),
        ),
    )

    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit["stage"] == "site_matrix")
        & (preprocessed.row_audit["action"] == "dropped")
    ]
    assert dropped.shape[0] == 1
    assert dropped.iloc[0]["source_row_id"] == "SRC_ROW_1"
    assert dropped.iloc[0]["site_id"] == _site_key("MAPK14;Y182;")
    assert "site_sequence is missing or blank" in str(dropped.iloc[0]["reason"])


def test_site_matrix_stage_audits_incomplete_value_drops() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, float("nan")],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = _with_site_identity(
        pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "AKT1"],
                "site": ["Y182", "T308"],
                "site_sequence": ["SEQ_A", "SEQ_R"],
                "localisation_confidence": [0.95, 0.9],
            },
            index=phospho.index.copy(),
        )
    )

    preprocessed = DatasetPreprocessor().run(
        phospho=phospho,
        site_metadata=site_metadata,
        sample_metadata=None,
        total=None,
        plan=replace(
            PreprocessingPlan.from_config(
                DatasetPreprocessingConfig(
                    site_matrix=DatasetSiteMatrixConfig(policy="build_from_metadata")
                )
            ),
            stage_order=("site_matrix",),
        ),
    )

    dropped = preprocessed.row_audit.loc[
        (preprocessed.row_audit["stage"] == "site_matrix")
        & (preprocessed.row_audit["action"] == "dropped")
    ]
    assert dropped.shape[0] == 1
    assert dropped.iloc[0]["source_row_id"] == "SRC_ROW_2"
    assert dropped.iloc[0]["site_id"] == _site_key("AKT1;T308;")
    assert dropped.iloc[0]["reason"] == "dropped by site_matrix missing-data policy"


def test_site_matrix_stage_audits_duplicate_resolution_first_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = _with_site_identity(
        pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "MAPK14"],
                "site": ["Y182", "Y182"],
                "site_sequence": ["SEQ_A", "SEQ_R"],
                "localisation_confidence": [0.95, 0.9],
            },
            index=phospho.index.copy(),
        )
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
                    duplicate_site_policy="first",
                )
            )
        ),
    )

    site_matrix_audit = preprocessed.row_audit.loc[
        preprocessed.row_audit["stage"] == "site_matrix"
    ].set_index("source_row_id")
    assert site_matrix_audit.loc["SRC_ROW_1", "action"] == "retained"
    assert bool(site_matrix_audit.loc["SRC_ROW_1", "retained"]) is True
    assert site_matrix_audit.loc["SRC_ROW_2", "action"] == "collapsed"
    assert bool(site_matrix_audit.loc["SRC_ROW_2", "retained"]) is False
    assert site_matrix_audit.loc["SRC_ROW_2", "retained_row_id"] == "SRC_ROW_1"
    snapshot = site_matrix_audit.loc["SRC_ROW_2", "parameter_snapshot"]
    assert isinstance(snapshot, dict)
    assert snapshot["duplicate_site_policy"] == "first"
    assert snapshot["site_matrix_duplicate_site_policy"] == "first"
    assert "duplicate_site_strategy" not in snapshot
    assert "site_matrix_duplicate_site_strategy" not in snapshot


def test_site_matrix_stage_audits_duplicate_resolution_max_mean_signal_policy() -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 5.0],
            "sample_b": [2.0, 6.0],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = _with_site_identity(
        pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "MAPK14"],
                "site": ["Y182", "Y182"],
                "site_sequence": ["SEQ_A", "SEQ_R"],
                "localisation_confidence": [0.95, 0.9],
            },
            index=phospho.index.copy(),
        )
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
                    duplicate_site_policy="max_mean_signal",
                )
            )
        ),
    )

    site_matrix_audit = preprocessed.row_audit.loc[
        preprocessed.row_audit["stage"] == "site_matrix"
    ].set_index("source_row_id")
    assert site_matrix_audit.loc["SRC_ROW_2", "action"] == "retained"
    assert bool(site_matrix_audit.loc["SRC_ROW_2", "retained"]) is True
    assert site_matrix_audit.loc["SRC_ROW_1", "action"] == "collapsed"
    assert bool(site_matrix_audit.loc["SRC_ROW_1", "retained"]) is False
    assert site_matrix_audit.loc["SRC_ROW_1", "retained_row_id"] == "SRC_ROW_2"


@pytest.mark.parametrize(
    "duplicate_policy",
    ("aggregate_mean", "aggregate_median"),
)
def test_site_matrix_stage_audits_aggregate_duplicate_contributors(
    duplicate_policy: str,
) -> None:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
        },
        index=pd.Index(["SRC_ROW_1", "SRC_ROW_2"], name="source_row"),
    )
    site_metadata = _with_site_identity(
        pd.DataFrame(
            {
                "gene_symbol": ["MAPK14", "MAPK14"],
                "site": ["Y182", "Y182"],
                "site_sequence": ["SEQ_A", "SEQ_R"],
                "localisation_confidence": [0.95, 0.9],
            },
            index=phospho.index.copy(),
        )
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
                    duplicate_site_policy=duplicate_policy,
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
    assert set(aggregated["retained_row_id"]) == {_site_key("MAPK14;Y182;")}
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
            "localisation_confidence": [0.95],
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
    assert tuple(preprocessed.row_audit.columns) == ROW_AUDIT_COLUMNS
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
