from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import (
    DatasetBuildRequest,
    DatasetGroupCoverageFilterConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    Organism,
)
from phospy.errors.input import PhosPyInputError


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "control_1": [1.0, 1.0, 1.0, float("nan")],
            "control_2": [2.0, 2.0, float("nan"), float("nan")],
            "control_3": [3.0, float("nan"), float("nan"), float("nan")],
            "treated_1": [float("nan"), float("nan"), 4.0, 7.0],
            "treated_2": [float("nan"), float("nan"), float("nan"), 8.0],
            "treated_3": [float("nan"), float("nan"), float("nan"), float("nan")],
        },
        index=pd.Index(
            ["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;", "PRKACA;S339;"],
            name="site_id",
        ),
    )


def _complete_phospho() -> pd.DataFrame:
    phospho = _phospho().copy(deep=True)
    values = [
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        [2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        [3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        [4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
    ]
    return pd.DataFrame(values, index=phospho.index.copy(), columns=phospho.columns)


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B", "PRKACA"],
            "protein_id": ["MAPK14", "AKT1", "GSK3B", "PRKACA"],
            "site": ["Y182", "T308", "S9", "S339"],
            "site_sequence": [
                ("A" * 15) + "Y" + ("A" * 15),
                ("A" * 15) + "T" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
                ("A" * 15) + "S" + ("A" * 15),
            ],
            "localisation_confidence": [0.95, 0.94, 0.93, 0.92],
        },
        index=index.copy(),
    )


def _sample_metadata(index: pd.Index | None = None) -> pd.DataFrame:
    sample_index = _phospho().columns.copy() if index is None else index
    conditions = {
        "control_1": "control",
        "control_2": "control",
        "control_3": "control",
        "treated_1": "treated",
        "treated_2": "treated",
        "treated_3": "treated",
    }
    return pd.DataFrame(
        {"condition": [conditions[str(sample)] for sample in sample_index]},
        index=sample_index.copy(),
    )


def _coverage_config(
    *,
    min_finite_observations_per_group: int | None = 2,
    min_finite_fraction_per_group: float | None = None,
    min_groups_passing_threshold: int = 1,
) -> DatasetPreprocessingConfig:
    return DatasetPreprocessingConfig(
        group_coverage_filter=DatasetGroupCoverageFilterConfig(
            enabled=True,
            group_column="condition",
            min_finite_observations_per_group=min_finite_observations_per_group,
            min_finite_fraction_per_group=min_finite_fraction_per_group,
            min_groups_passing_threshold=min_groups_passing_threshold,
        ),
        missing_data=DatasetMissingDataConfig(
            policy="impute_row_median",
            min_observed_values=1,
        ),
    )


def _build_dataset(
    *,
    phospho: pd.DataFrame | None = None,
    sample_metadata: pd.DataFrame | None = None,
    preprocessing_config: DatasetPreprocessingConfig | None = None,
):
    resolved_phospho = _phospho() if phospho is None else phospho
    return AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=resolved_phospho,
            site_metadata=_site_metadata(resolved_phospho.index),
            sample_metadata=(
                _sample_metadata(resolved_phospho.columns)
                if sample_metadata is None
                else sample_metadata
            ),
            organism=Organism.RAT,
            input_intensity_scale="linear",
            preprocessing_config=(
                _coverage_config()
                if preprocessing_config is None
                else preprocessing_config
            ),
        )
    )


def test_group_coverage_filter_retains_sites_passing_in_one_group() -> None:
    dataset = _build_dataset()

    display_ids = set(dataset.site_metadata.loc[:, "display_id"].astype(str))
    assert display_ids == {"MAPK14;Y182;", "AKT1;T308;", "PRKACA;S339;"}
    assert int(dataset.phospho.isna().to_numpy().sum()) == 0

    assert dataset.preprocessing_report is not None
    report = dataset.preprocessing_report
    coverage_counts = report.row_counts.loc[
        report.row_counts.loc[:, "stage"] == "group_coverage_filter"
    ].iloc[0]
    assert int(coverage_counts["input_rows"]) == 4
    assert int(coverage_counts["output_rows"]) == 3
    assert int(coverage_counts["dropped_rows"]) == 1

    operation = report.operations.loc[
        report.operations.loc[:, "stage"] == "group_coverage_filter"
    ].iloc[0]
    parameters = operation["parameters"]
    assert parameters["group_column"] == "condition"
    assert parameters["threshold_type"] == "count"
    assert parameters["min_finite_observations_per_group"] == 2
    assert parameters["min_groups_passing_threshold"] == 1
    summary = parameters["execution_summary"]["diagnostic_summary"]
    assert summary["input_feature_count"] == 4
    assert summary["retained_feature_count"] == 3
    assert summary["removed_feature_count"] == 1
    assert summary["group_column"] == "condition"

    dropped = report.row_audit.loc[
        report.row_audit.loc[:, "stage"] == "group_coverage_filter"
    ]
    assert dropped.shape[0] == 1
    assert str(dropped.iloc[0]["action"]) == "dropped"
    assert "insufficient finite coverage" in str(dropped.iloc[0]["reason"])


def test_group_coverage_filter_count_threshold_boundary_is_inclusive() -> None:
    dataset = _build_dataset()

    display_ids = set(dataset.site_metadata.loc[:, "display_id"].astype(str))
    assert "AKT1;T308;" in display_ids


def test_group_coverage_filter_fraction_threshold_boundary_is_inclusive() -> None:
    dataset = _build_dataset(
        preprocessing_config=_coverage_config(
            min_finite_observations_per_group=None,
            min_finite_fraction_per_group=2.0 / 3.0,
        )
    )

    display_ids = set(dataset.site_metadata.loc[:, "display_id"].astype(str))
    assert display_ids == {"MAPK14;Y182;", "AKT1;T308;", "PRKACA;S339;"}


def test_group_coverage_filter_fails_when_all_sites_are_dropped() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="group_coverage_filter removed all phosphosites/features",
    ):
        _build_dataset(
            preprocessing_config=_coverage_config(
                min_finite_observations_per_group=3,
                min_groups_passing_threshold=2,
            )
        )


def test_group_coverage_filter_rejects_missing_group_column() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="group_column references missing sample_metadata column 'condition'",
    ):
        _build_dataset(sample_metadata=pd.DataFrame(index=_phospho().columns.copy()))


def test_group_coverage_filter_rejects_sample_metadata_mismatch() -> None:
    mismatched_metadata = _sample_metadata(
        pd.Index(
            [
                "control_1",
                "control_2",
                "control_3",
                "treated_1",
                "treated_2",
            ],
            name="sample_id",
        )
    )

    with pytest.raises(
        PhosPyInputError,
        match="missing rows for coverage-filter samples.*'treated_3'",
    ):
        _build_dataset(sample_metadata=mismatched_metadata)


def test_disabled_group_coverage_filter_leaves_dataset_building_unchanged() -> None:
    phospho = _complete_phospho()
    dataset = _build_dataset(
        phospho=phospho,
        preprocessing_config=DatasetPreprocessingConfig(),
    )

    assert set(dataset.site_metadata.loc[:, "display_id"].astype(str)) == {
        "MAPK14;Y182;",
        "AKT1;T308;",
        "GSK3B;S9;",
        "PRKACA;S339;",
    }
    assert dataset.preprocessing_report is not None
    assert "group_coverage_filter" not in set(
        dataset.preprocessing_report.operations.loc[:, "stage"].astype(str)
    )
