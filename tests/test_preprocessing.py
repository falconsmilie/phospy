from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset, PhosphoDataset, SimpleKinaseWorkflow
from phospy.datasets import DatasetLoader, DatasetSchema
from phospy.errors import (
    InputCompatibilityError,
    PhospyValidationError,
    TableSchemaError,
)
from phospy.preprocessing import (
    CoverageFilterResult,
    DatasetPreprocessing,
    LocalizationFilterResult,
    PhosphoPreprocessor,
    ProteinCorrectionResult,
    ProteinCorrectionService,
    TotalPreprocessor,
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_localized_sites,
    filter_min_observed,
    filter_sites_by_coverage,
    replace_sentinel_with_nan,
)
from phospy.preprocessing.analysis_ready import build_analysis_ready_dataset


def test_filter_localized_sites_filters_inclusive_threshold_and_returns_summary() -> (
    None
):
    df = pd.DataFrame(
        {
            "gene_p_site": ["PRKACA_S339", "BTK_Y551", "LYN_Y397"],
            "localization_prob": [0.75, 0.7499, 1.0],
        }
    )

    result = filter_localized_sites(df, threshold=0.75, return_summary=True)

    assert isinstance(result, LocalizationFilterResult)
    assert result.filtered["gene_p_site"].tolist() == ["PRKACA_S339", "LYN_Y397"]
    assert result.summary.input_rows == 3
    assert result.summary.retained_rows == 2
    assert result.summary.removed_rows == 1
    assert result.summary.retention_fraction == pytest.approx(2 / 3)
    assert result.summary.threshold == 0.75
    assert result.summary.localization_col == "localization_prob"


def test_filter_localized_sites_accepts_threshold_boundary_of_one() -> None:
    df = pd.DataFrame(
        {
            "localization_prob": [1.0, 0.999],
            "gene_p_site": ["A", "B"],
        }
    )

    filtered = filter_localized_sites(df, threshold=1.0)

    assert filtered["gene_p_site"].tolist() == ["A"]


def test_filter_localized_sites_rejects_invalid_threshold() -> None:
    df = pd.DataFrame({"localization_prob": [0.95]})

    with pytest.raises(
        PhospyValidationError,
        match="threshold must be a finite numeric value between 0 and 1",
    ):
        filter_localized_sites(df, threshold=1.1)


def test_filter_localized_sites_reports_missing_required_column() -> None:
    df = pd.DataFrame({"score_for_localization": [0.95]})

    with pytest.raises(
        TableSchemaError,
        match=r"filter_localized_sites\(\) input is missing required columns: localization_prob",
    ):
        filter_localized_sites(df)


def test_filter_sites_by_coverage_filters_rows_and_returns_summary() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["A", "B", "C"],
            "p_group1": [1.0, 1.0, np.nan],
            "p_group2": [2.0, np.nan, np.nan],
            "p_group3": [3.0, np.nan, 5.0],
            "p_group4": [4.0, 4.0, np.nan],
        }
    )

    result = filter_sites_by_coverage(
        df,
        columns=["p_group1", "p_group2", "p_group3", "p_group4"],
        min_coverage=0.5,
        return_summary=True,
    )

    assert isinstance(result, CoverageFilterResult)
    assert result.filtered["gene_p_site"].tolist() == ["A", "B"]
    assert result.summary.input_rows == 3
    assert result.summary.retained_rows == 2
    assert result.summary.removed_rows == 1
    assert result.summary.retention_fraction == pytest.approx(2 / 3)
    assert result.summary.min_coverage == 0.5
    assert result.summary.required_observed_count == 2
    assert result.summary.value_columns == (
        "p_group1",
        "p_group2",
        "p_group3",
        "p_group4",
    )


def test_filter_sites_by_coverage_applies_inclusive_boundary_threshold() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["A", "B"],
            "p_group1": [1.0, 1.0],
            "p_group2": [2.0, np.nan],
            "p_group3": [np.nan, np.nan],
            "p_group4": [4.0, np.nan],
        }
    )

    filtered = filter_sites_by_coverage(
        df,
        columns=["p_group1", "p_group2", "p_group3", "p_group4"],
        min_coverage=0.75,
    )

    assert filtered["gene_p_site"].tolist() == ["A"]


def test_filter_sites_by_coverage_rejects_invalid_threshold() -> None:
    df = pd.DataFrame({"p_group1": [1.0], "p_group2": [2.0]})

    with pytest.raises(
        PhospyValidationError,
        match="min_coverage must be a finite numeric value between 0 and 1",
    ):
        filter_sites_by_coverage(
            df,
            columns=["p_group1", "p_group2"],
            min_coverage=-0.1,
        )


def test_filter_sites_by_coverage_reports_missing_required_column() -> None:
    df = pd.DataFrame({"p_group1": [1.0]})

    with pytest.raises(
        TableSchemaError,
        match=r"filter_sites_by_coverage\(\) input is missing required columns: p_group2",
    ):
        filter_sites_by_coverage(
            df,
            columns=["p_group1", "p_group2"],
            min_coverage=0.5,
        )


def test_filter_sites_by_coverage_returns_deterministic_empty_result() -> None:
    df = pd.DataFrame(
        {
            "gene_p_site": ["A", "B"],
            "p_group1": [1.0, np.nan],
            "p_group2": [np.nan, np.nan],
        }
    )

    result = filter_sites_by_coverage(
        df,
        columns=["p_group1", "p_group2"],
        min_coverage=1.0,
        return_summary=True,
    )

    assert result.filtered.empty
    assert result.filtered.columns.tolist() == df.columns.tolist()
    assert result.summary.input_rows == 2
    assert result.summary.retained_rows == 0
    assert result.summary.removed_rows == 2
    assert result.summary.retention_fraction == 0.0
    assert result.summary.required_observed_count == 2


def test_filter_sites_by_coverage_requires_at_least_one_column() -> None:
    df = pd.DataFrame({"gene_p_site": ["A"]})

    with pytest.raises(
        PhospyValidationError,
        match=(
            r"filter_sites_by_coverage\(\) requires at least one column name in 'columns'"
        ),
    ):
        filter_sites_by_coverage(df, columns=[])


def test_filter_min_observed_rejects_negative_threshold() -> None:
    df = pd.DataFrame({"gene": ["A"], "x1": [1.0]})

    with pytest.raises(
        PhospyValidationError,
        match="min_observed must be a non-negative integer",
    ):
        filter_min_observed(df, ["x1"], min_observed=-1)


def test_filter_min_observed_reports_missing_required_column() -> None:
    df = pd.DataFrame({"gene": ["A"], "x1": [1.0]})

    with pytest.raises(
        TableSchemaError,
        match=r"filter_min_observed\(\) input is missing required columns: x2",
    ):
        filter_min_observed(df, ["x1", "x2"], min_observed=1)


def test_replace_sentinel_with_nan_reports_missing_required_column() -> None:
    df = pd.DataFrame({"gene": ["A"], "x1": [1.0]})

    with pytest.raises(
        TableSchemaError,
        match=r"replace_sentinel_with_nan\(\) input is missing required columns: x2",
    ):
        replace_sentinel_with_nan(df, ["x1", "x2"], sentinel=12)


def test_replace_sentinel_with_nan_and_filter_min_observed() -> None:
    df = pd.DataFrame(
        {
            "gene": ["A", "B"],
            "x1": [1.0, 12.0],
            "x2": [2.0, np.nan],
            "x3": [12.0, 3.0],
        }
    )
    cleaned = replace_sentinel_with_nan(df, ["x1", "x2", "x3"], sentinel=12)
    assert cleaned["x1"].isna().sum() == 1
    filtered = filter_min_observed(cleaned, ["x1", "x2", "x3"], min_observed=2)
    assert filtered["gene"].tolist() == ["A"]


def test_replace_sentinel_with_nan_uses_block_numeric_conversion(monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "x1": [1.0, 12.0],
            "x2": [2.0, 12.0],
            "x3": [3.0, 4.0],
        }
    )

    dataframe_astype_calls = 0
    series_astype_calls = 0
    original_dataframe_astype = pd.DataFrame.astype
    original_series_astype = pd.Series.astype

    def counting_dataframe_astype(
        self, dtype=None, copy: bool = False, errors: str = "raise"
    ):
        nonlocal dataframe_astype_calls
        dataframe_astype_calls += 1
        return original_dataframe_astype(self, dtype=dtype, errors=errors)

    def counting_series_astype(
        self, dtype, copy: bool | None = None, errors: str = "raise"
    ):
        nonlocal series_astype_calls
        series_astype_calls += 1
        return original_series_astype(self, dtype=dtype, errors=errors)

    monkeypatch.setattr(pd.DataFrame, "astype", counting_dataframe_astype)
    monkeypatch.setattr(pd.Series, "astype", counting_series_astype)

    cleaned = replace_sentinel_with_nan(df, ["x1", "x2", "x3"], sentinel=12)

    assert dataframe_astype_calls == 1
    assert series_astype_calls == 0
    assert cleaned.isna().sum().to_dict() == {"x1": 1, "x2": 1, "x3": 0}


def test_collapse_duplicate_genes_keeps_highest_mean_signal() -> None:
    df = pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca", "Btk"],
            "group1": [1.0, 5.0, 2.0],
            "group2": [1.0, 5.0, 2.0],
            "group3": [1.0, 5.0, 2.0],
            "group4": [1.0, 5.0, 2.0],
            "group5": [1.0, 5.0, 2.0],
            "group6": [1.0, 5.0, 2.0],
        }
    )
    out = collapse_duplicate_genes(
        df=df,
        gene_col="genes",
        value_cols=["group1", "group2", "group3", "group4", "group5", "group6"],
    )
    prkaca = out.loc[out["genes"] == "PRKACA"].iloc[0]
    assert prkaca["group1"] == 5.0
    assert sorted(out["genes"].tolist()) == ["BTK", "PRKACA"]


def test_collapse_duplicate_genes_prefers_more_observed_values_before_mean() -> None:
    df = pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca"],
            "group1": [10.0, 6.0],
            "group2": [np.nan, 6.0],
            "group3": [np.nan, 6.0],
        }
    )

    out = collapse_duplicate_genes(
        df=df,
        gene_col="genes",
        value_cols=["group1", "group2", "group3"],
    )

    prkaca = out.loc[out["genes"] == "PRKACA"].iloc[0]
    assert prkaca["group1"] == 6.0
    assert prkaca["group2"] == 6.0
    assert prkaca["group3"] == 6.0


def test_collapse_duplicate_genes_collapses_mixed_case_duplicates_before_ranking() -> (
    None
):
    df = pd.DataFrame(
        {
            "genes": ["akt1", "AKT1", "Mapk1"],
            "group1": [1.0, 5.0, 3.0],
            "group2": [1.0, 5.0, 3.0],
        }
    )

    out = collapse_duplicate_genes(
        df=df,
        gene_col="genes",
        value_cols=["group1", "group2"],
    )

    assert out["genes"].tolist() == ["AKT1", "MAPK1"]
    akt1 = out.loc[out["genes"] == "AKT1"].iloc[0]
    assert akt1["group1"] == 5.0
    assert akt1["group2"] == 5.0


def test_collapse_duplicate_genes_uses_original_order_as_final_tiebreaker() -> None:
    df = pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca"],
            "row_id": ["first", "second"],
            "group1": [4.0, 4.0],
            "group2": [2.0, 2.0],
        }
    )

    out = collapse_duplicate_genes(
        df=df,
        gene_col="genes",
        value_cols=["group1", "group2"],
        uppercase=False,
    )

    prkaca = out.loc[out["genes"] == "Prkaca"].iloc[0]
    assert prkaca["row_id"] == "first"
    assert prkaca["group1"] == 4.0
    assert prkaca["group2"] == 2.0


def test_collapse_duplicate_genes_drops_all_nan_groups_before_deduping() -> None:
    df = pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca", "Btk"],
            "group1": [np.nan, np.nan, 2.0],
            "group2": [np.nan, np.nan, 3.0],
        }
    )

    out = collapse_duplicate_genes(
        df=df,
        gene_col="genes",
        value_cols=["group1", "group2"],
    )

    assert out["genes"].tolist() == ["BTK"]
    assert out[["group1", "group2"]].iloc[0].tolist() == [2.0, 3.0]


def test_collapse_duplicate_genes_reports_missing_required_columns() -> None:
    df = pd.DataFrame({"genes": ["Prkaca"], "group1": [1.0]})

    with pytest.raises(TableSchemaError, match="missing required columns: group2"):
        collapse_duplicate_genes(
            df=df,
            gene_col="genes",
            value_cols=["group1", "group2"],
        )


def test_correct_phospho_to_protein_and_pairwise_comparisons() -> None:
    phospho = pd.DataFrame(
        {
            "gene_names": ["PRKACA"],
            "p_group1": [8.0],
            "p_group2": [7.0],
            "p_group3": [6.0],
            "p_group4": [5.0],
            "p_group5": [4.0],
            "p_group6": [3.0],
        }
    )
    total = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )
    corrected = correct_phospho_to_protein(
        phospho,
        total,
        phospho_gene_col="gene_names",
        total_gene_col="genes",
        phospho_cols=[
            "p_group1",
            "p_group2",
            "p_group3",
            "p_group4",
            "p_group5",
            "p_group6",
        ],
        protein_cols=["group1", "group2", "group3", "group4", "group5", "group6"],
    )
    assert corrected["phospho_corrected_1"].iloc[0] == 7.0
    with_comparisons = add_pairwise_comparisons(
        corrected,
        comparisons=[("group1", "group4")],
    )
    assert with_comparisons["p_group1_group4"].iloc[0] == 3.0


def test_correct_phospho_to_protein_returns_summary_with_unmatched_rows() -> None:
    phospho = pd.DataFrame(
        {
            "gene_names": ["PRKACA", "MISSING"],
            "p_group1": [8.0, 6.0],
        }
    )
    total = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
        }
    )

    result = correct_phospho_to_protein(
        phospho,
        total,
        phospho_gene_col="gene_names",
        total_gene_col="genes",
        phospho_cols=["p_group1"],
        protein_cols=["group1"],
        return_summary=True,
    )

    assert isinstance(result, ProteinCorrectionResult)
    assert result.corrected["gene_names"].tolist() == ["PRKACA"]
    assert result.corrected["phospho_corrected_1"].tolist() == [7.0]
    assert result.summary.input_rows == 2
    assert result.summary.matched_rows == 1
    assert result.summary.unmatched_rows == 1
    assert result.summary.unmatched_fraction == pytest.approx(0.5)
    assert result.summary.phospho_gene_col == "gene_names"
    assert result.summary.total_gene_col == "genes"
    assert result.summary.unmatched_gene_preview == ("MISSING",)


def test_correct_phospho_to_protein_supports_strict_unmatched_row_mode() -> None:
    phospho = pd.DataFrame(
        {
            "gene_names": ["PRKACA", "MISSING"],
            "p_group1": [8.0, 6.0],
        }
    )
    total = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
        }
    )

    with pytest.raises(
        InputCompatibilityError,
        match=r"would drop 1 of 2 phosphosite rows \(50.0%\)",
    ):
        correct_phospho_to_protein(
            phospho,
            total,
            phospho_gene_col="gene_names",
            total_gene_col="genes",
            phospho_cols=["p_group1"],
            protein_cols=["group1"],
            max_unmatched_fraction=0.0,
        )


def test_correct_phospho_to_protein_reports_missing_required_value_columns() -> None:
    phospho = pd.DataFrame({"gene_names": ["PRKACA"]})
    total = pd.DataFrame({"genes": ["PRKACA"], "group1": [1.0]})

    with pytest.raises(
        TableSchemaError,
        match=(
            r"correct_phospho_to_protein\(\) phospho input is missing required columns: p_group1"
        ),
    ):
        correct_phospho_to_protein(
            phospho,
            total,
            phospho_gene_col="gene_names",
            total_gene_col="genes",
            phospho_cols=["p_group1"],
            protein_cols=["group1"],
        )


def test_correct_phospho_to_protein_rejects_duplicate_total_genes() -> None:
    phospho = pd.DataFrame({"gene_names": ["PRKACA"], "p_group1": [8.0]})
    total = pd.DataFrame({"genes": ["PRKACA", "PRKACA"], "group1": [1.0, 2.0]})

    with pytest.raises(
        InputCompatibilityError, match="must be unique before protein correction"
    ):
        correct_phospho_to_protein(
            phospho,
            total,
            phospho_gene_col="gene_names",
            total_gene_col="genes",
            phospho_cols=["p_group1"],
            protein_cols=["group1"],
        )


def test_correct_phospho_to_protein_matches_after_identifier_normalization() -> None:
    phospho = pd.DataFrame(
        {
            "gene_names": [" prkaca ", "BTK"],
            "p_group1": [8.0, 6.0],
        }
    )
    total = pd.DataFrame(
        {
            "genes": ["PRKACA", " btk "],
            "group1": [1.0, 2.0],
        }
    )

    corrected = correct_phospho_to_protein(
        phospho,
        total,
        phospho_gene_col="gene_names",
        total_gene_col="genes",
        phospho_cols=["p_group1"],
        protein_cols=["group1"],
    )

    assert corrected["gene_names"].tolist() == [" prkaca ", "BTK"]
    assert corrected["phospho_corrected_1"].tolist() == [7.0, 4.0]
    assert "genes" not in corrected.columns
    assert "__phospy_normalized_phospho_gene_key" not in corrected.columns
    assert "__phospy_normalized_total_gene_key" not in corrected.columns


def test_correct_phospho_to_protein_rejects_duplicate_normalized_total_genes() -> None:
    phospho = pd.DataFrame({"gene_names": ["PRKACA"], "p_group1": [8.0]})
    total = pd.DataFrame(
        {
            "genes": ["PRKACA", " prkaca "],
            "group1": [1.0, 2.0],
        }
    )

    with pytest.raises(
        InputCompatibilityError, match="must be unique before protein correction"
    ):
        correct_phospho_to_protein(
            phospho,
            total,
            phospho_gene_col="gene_names",
            total_gene_col="genes",
            phospho_cols=["p_group1"],
            protein_cols=["group1"],
        )


def test_total_preprocessor_matches_existing_total_preparation_flow() -> None:
    total_df = pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca", "Btk"],
            "group1": [1.0, 5.0, 9.0],
            "group2": [1.0, 5.0, 9.0],
            "group3": [1.0, 5.0, 9.0],
            "group4": [1.0, 5.0, 9.0],
            "group5": [1.0, 5.0, 9.0],
            "group6": [1.0, 5.0, 9.0],
        }
    )

    total_unique, total_filtered = TotalPreprocessor(schema=DatasetSchema()).prepare(
        total_df
    )

    assert total_unique["genes"].tolist() == ["BTK", "PRKACA"]
    assert total_filtered["genes"].tolist() == ["BTK", "PRKACA"]
    assert total_unique.loc[total_unique["genes"] == "PRKACA", "group1"].iloc[0] == 5.0


def test_total_preprocessor_does_not_mutate_input_dataframe() -> None:
    total_df = pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca", "Btk"],
            "group1": [1.0, 5.0, 12.0],
            "group2": [1.0, 5.0, 9.0],
            "group3": [1.0, 5.0, 9.0],
            "group4": [1.0, 5.0, 9.0],
            "group5": [1.0, 5.0, 9.0],
            "group6": [1.0, 5.0, 9.0],
        }
    )
    original = total_df.copy(deep=True)

    TotalPreprocessor(schema=DatasetSchema()).prepare(total_df, sentinel=12.0)

    pd.testing.assert_frame_equal(total_df, original)


def test_phospho_preprocessor_applies_localization_sentinel_and_coverage_rules() -> (
    None
):
    phospho_df = pd.DataFrame(
        {
            "gene_names": ["prkaca", "btk", "lyn"],
            "gene_p_site": ["prkaca_s339", "btk_y551", "lyn_y397"],
            "localization_prob": [0.95, 0.70, 0.95],
            "p_group1": [8.0, 6.0, 12.0],
            "p_group2": [7.0, 5.0, 12.0],
            "p_group3": [6.0, 4.0, 12.0],
            "p_group4": [5.0, 3.0, 12.0],
            "p_group5": [4.0, 2.0, 2.0],
            "p_group6": [3.0, 1.0, 12.0],
        }
    )

    filtered = PhosphoPreprocessor(schema=DatasetSchema()).prepare(
        phospho_df,
        sentinel=12.0,
        min_observed=4,
    )

    assert filtered["gene_names"].tolist() == ["PRKACA"]
    assert filtered["gene_p_site"].tolist() == ["prkaca_s339"]


def test_phospho_preprocessor_does_not_mutate_input_dataframe() -> None:
    phospho_df = pd.DataFrame(
        {
            "gene_names": ["prkaca", "btk"],
            "gene_p_site": ["prkaca_s339", "btk_y551"],
            "localization_prob": [0.95, 0.70],
            "p_group1": [8.0, 12.0],
            "p_group2": [7.0, 12.0],
            "p_group3": [6.0, 12.0],
            "p_group4": [5.0, 12.0],
            "p_group5": [4.0, 12.0],
            "p_group6": [3.0, 12.0],
        }
    )
    original = phospho_df.copy(deep=True)

    PhosphoPreprocessor(schema=DatasetSchema()).prepare(
        phospho_df,
        sentinel=12.0,
        min_observed=4,
    )

    pd.testing.assert_frame_equal(phospho_df, original)


def test_protein_correction_service_applies_correction_and_pairwise_augmentation() -> (
    None
):
    phospho_df = pd.DataFrame(
        {
            "gene_names": ["PRKACA"],
            "p_group1": [8.0],
            "p_group2": [7.0],
            "p_group3": [6.0],
            "p_group4": [5.0],
            "p_group5": [4.0],
            "p_group6": [3.0],
        }
    )
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )

    service = ProteinCorrectionService(
        schema=DatasetSchema(),
        comparisons=[("group1", "group4")],
    )
    corrected = service.correct(phospho_df, total_df)
    with_comparisons = service.add_pairwise_comparisons(corrected)

    assert corrected["phospho_corrected_1"].iloc[0] == 7.0
    assert with_comparisons["p_group1_group4"].iloc[0] == 3.0


def test_protein_correction_service_does_not_route_through_public_facade(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.preprocessing.services as preprocessing_services_module

    phospho_df = pd.DataFrame(
        {
            "gene_names": ["PRKACA"],
            "p_group1": [8.0],
            "p_group2": [7.0],
            "p_group3": [6.0],
            "p_group4": [5.0],
            "p_group5": [4.0],
            "p_group6": [3.0],
        }
    )
    total_df = pd.DataFrame(
        {
            "genes": ["PRKACA"],
            "group1": [1.0],
            "group2": [1.0],
            "group3": [1.0],
            "group4": [1.0],
            "group5": [1.0],
            "group6": [1.0],
        }
    )

    def _boom(*args: object, **kwargs: object) -> object:
        raise AssertionError("service should not call the public preprocessing facade")

    monkeypatch.setattr(
        preprocessing_services_module,
        "correct_phospho_to_protein",
        _boom,
        raising=False,
    )

    service = ProteinCorrectionService(schema=DatasetSchema())

    corrected = service.correct(phospho_df, total_df)

    assert corrected["phospho_corrected_1"].iloc[0] == 7.0


def test_add_pairwise_comparisons_uses_schema_group_names() -> None:
    corrected = pd.DataFrame(
        {
            "corrected_a": [7.0],
            "corrected_b": [4.0],
        }
    )
    schema = DatasetSchema(
        total_cols=("sample_a", "sample_b"),
        phospho_cols=("p_sample_a", "p_sample_b"),
        corrected_cols=("corrected_a", "corrected_b"),
    )

    with_comparisons = add_pairwise_comparisons(
        corrected,
        comparisons=[("sample_a", "sample_b")],
        schema=schema,
    )

    assert with_comparisons["p_sample_a_sample_b"].iloc[0] == 3.0


def test_add_pairwise_comparisons_rejects_self_comparisons_with_custom_mapping() -> (
    None
):
    corrected = pd.DataFrame(
        {
            "corrected_a": [7.0],
            "corrected_b": [4.0],
        }
    )

    with pytest.raises(InputCompatibilityError, match="Self comparison pair"):
        add_pairwise_comparisons(
            corrected,
            comparisons=[("sample_a", "sample_a")],
            group_to_corrected_col={
                "sample_a": "corrected_a",
                "sample_b": "corrected_b",
            },
        )


def test_add_pairwise_comparisons_rejects_reverse_duplicate_pairs_with_custom_mapping() -> (
    None
):
    corrected = pd.DataFrame(
        {
            "corrected_a": [7.0],
            "corrected_b": [4.0],
        }
    )

    with pytest.raises(
        InputCompatibilityError,
        match="Duplicate comparison pair regardless of direction",
    ):
        add_pairwise_comparisons(
            corrected,
            comparisons=[("sample_a", "sample_b"), ("sample_b", "sample_a")],
            group_to_corrected_col={
                "sample_a": "corrected_a",
                "sample_b": "corrected_b",
            },
        )


def test_dataset_preprocessing_run_rejects_mixed_config_styles() -> None:
    from phospy import PhosphoDataset
    from phospy.preprocessing import CorePreprocessingConfig

    dataset = PhosphoDataset(
        total_df=pd.DataFrame(
            {
                "genes": ["PRKACA"],
                "group1": [1.0],
                "group2": [1.0],
                "group3": [1.0],
                "group4": [1.0],
                "group5": [1.0],
                "group6": [1.0],
            }
        ),
        phospho_df=pd.DataFrame(
            {
                "uid": ["u1"],
                "gene_names": ["PRKACA"],
                "gene_p_site": ["PRKACA_S339"],
                "localization_prob": [0.95],
                "centralized_sequence": ["AAAAAA"],
                "p_group1": [1.0],
                "p_group2": [1.0],
                "p_group3": [1.0],
                "p_group4": [1.0],
                "p_group5": [1.0],
                "p_group6": [1.0],
            }
        ),
    )

    with pytest.raises(
        ValueError,
        match=(
            r"DatasetPreprocessing.run\(\): pass either config or scalar "
            r"preprocessing options, not both\."
        ),
    ):
        dataset.preprocessing.run(
            config=CorePreprocessingConfig(),
            min_observed=1,
        )


def test_dataset_preprocessing_facade_shares_live_workspace_state_with_dataset() -> (
    None
):
    from phospy import PhosphoDataset

    dataset = PhosphoDataset(
        total_df=pd.DataFrame(
            {
                "genes": ["PRKACA", "BTK"],
                "group1": [1.0, 2.0],
                "group2": [1.0, 2.0],
                "group3": [1.0, 2.0],
                "group4": [1.0, 2.0],
                "group5": [1.0, 2.0],
                "group6": [1.0, 2.0],
            }
        ),
        phospho_df=pd.DataFrame(
            {
                "uid": ["u1", "u2"],
                "gene_names": ["PRKACA", "BTK"],
                "gene_p_site": ["PRKACA_S339", "BTK_Y551"],
                "localization_prob": [0.95, 0.95],
                "centralized_sequence": ["AAAAAA", "BBBBBB"],
                "p_group1": [1.0, 2.0],
                "p_group2": [1.0, 2.0],
                "p_group3": [1.0, 2.0],
                "p_group4": [1.0, 2.0],
                "p_group5": [1.0, 2.0],
                "p_group6": [1.0, 2.0],
            }
        ),
    )
    preprocessing = dataset.preprocessing

    dataset.total_df_live.loc[0, "group1"] = 111.0
    dataset.phospho_df_live.loc[0, "p_group1"] = 222.0

    assert preprocessing.total_df.loc[0, "group1"] == 111.0
    assert preprocessing.phospho_df.loc[0, "p_group1"] == 222.0

    preprocessing.total_df.loc[1, "group2"] = 333.0
    preprocessing.phospho_df.loc[1, "p_group2"] = 444.0

    assert dataset.total_df_live.loc[1, "group2"] == 333.0
    assert dataset.phospho_df_live.loc[1, "p_group2"] == 444.0


def test_preprocessing_result_wrappers_with_pandas_state_are_not_frozen_dataclasses() -> (
    None
):
    assert LocalizationFilterResult.__dataclass_params__.frozen is False
    assert CoverageFilterResult.__dataclass_params__.frozen is False
    assert ProteinCorrectionResult.__dataclass_params__.frozen is False


def test_preprocessing_facades_and_services_use_plain_dataclass_construction() -> None:
    assert DatasetPreprocessing.__dataclass_params__.frozen is False
    assert ProteinCorrectionService.__dataclass_params__.frozen is False


def test_dataset_preprocessing_run_analysis_ready_uses_example_fixture_data() -> None:
    example_dir = Path(__file__).resolve().parents[1] / "examples" / "data"
    dataset = PhosphoDataset.from_files(
        total_path=example_dir / "total.tsv",
        phospho_path=example_dir / "phospho.tsv",
        phospho_encoding="utf-16le",
    )

    analysis_ready = dataset.preprocessing.run_analysis_ready(
        max_unmatched_fraction=0.1,
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)
    expected = AnalysisReadyPhosphoDataset.from_core_processing_result(
        core,
        schema=dataset.schema,
        comparisons=dataset.comparisons,
        source="dataset preprocessing",
    )

    pd.testing.assert_frame_equal(
        analysis_ready.phospho_matrix,
        expected.phospho_matrix,
    )
    pd.testing.assert_frame_equal(
        analysis_ready.site_metadata,
        expected.site_metadata,
    )
    pd.testing.assert_series_equal(
        analysis_ready.site_sequences,
        expected.site_sequences,
    )
    pd.testing.assert_frame_equal(
        analysis_ready.phospho_corrected,
        expected.phospho_corrected,
    )
    assert analysis_ready.provenance == expected.provenance
    assert analysis_ready.phospho_matrix.index.equals(
        analysis_ready.site_metadata.index
    )
    assert analysis_ready.phospho_matrix.index.equals(
        analysis_ready.site_sequences.index
    )


def test_simple_kinase_workflow_reuses_bound_analysis_ready_adapter_on_fixture_files() -> (
    None
):
    fixture_dir = (
        Path(__file__).resolve().parents[1] / "examples" / "data" / "simple_workflow"
    )
    dataset = PhosphoDataset.from_files(
        total_path=fixture_dir / "total.tsv",
        phospho_path=fixture_dir / "phospho.tsv",
    )
    expected = dataset.run_analysis_ready(
        source="simple kinase workflow",
    )

    with SimpleKinaseWorkflow(flank_size=7).run(
        total=fixture_dir / "total.tsv",
        phospho=fixture_dir / "phospho.tsv",
        species="rat",
        min_substrates=1,
        min_motif_size=1,
        ensemble_size=2,
        top=3,
        inclusion=2,
        n_iterations=2,
        random_state=7,
        kinase_activity_threshold=0.1,
        kinase_activity_min_substrates=1,
        kinase_activity_top_n_substrates=3,
    ) as result:
        pd.testing.assert_frame_equal(
            result.analysis_ready_dataset.phospho_matrix,
            expected.phospho_matrix,
        )
        pd.testing.assert_frame_equal(
            result.analysis_ready_dataset.site_metadata,
            expected.site_metadata,
        )
        pd.testing.assert_series_equal(
            result.analysis_ready_dataset.site_sequences,
            expected.site_sequences,
        )
        pd.testing.assert_frame_equal(
            result.analysis_ready_dataset.phospho_corrected,
            expected.phospho_corrected,
        )
        assert result.analysis_ready_dataset.provenance == expected.provenance


def test_dataset_loader_resolve_inputs_converges_in_memory_file_and_mixed_sources() -> (
    None
):
    fixture_dir = (
        Path(__file__).resolve().parents[1] / "examples" / "data" / "simple_workflow"
    )
    total_path = fixture_dir / "total.tsv"
    phospho_path = fixture_dir / "phospho.tsv"
    total_df = pd.read_csv(total_path, sep="\t")
    phospho_df = pd.read_csv(phospho_path, sep="\t")
    loader = DatasetLoader(schema=DatasetSchema())

    in_memory_inputs = loader.resolve_inputs(total=total_df, phospho=phospho_df)
    file_inputs = loader.resolve_inputs(total=total_path, phospho=phospho_path)
    mixed_inputs = loader.resolve_inputs(total=total_path, phospho=phospho_df)

    pd.testing.assert_frame_equal(in_memory_inputs.total_df, file_inputs.total_df)
    pd.testing.assert_frame_equal(in_memory_inputs.phospho_df, file_inputs.phospho_df)
    pd.testing.assert_frame_equal(mixed_inputs.total_df, file_inputs.total_df)
    pd.testing.assert_frame_equal(mixed_inputs.phospho_df, file_inputs.phospho_df)


def test_analysis_ready_builder_full_inputs_reuses_dataset_preprocessing_seam(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.preprocessing.modes as preprocessing_modes_module

    fixture_dir = (
        Path(__file__).resolve().parents[1] / "examples" / "data" / "simple_workflow"
    )
    total_path = fixture_dir / "total.tsv"
    phospho_path = fixture_dir / "phospho.tsv"
    calls: list[dict[str, object]] = []
    original_run_analysis_ready = (
        preprocessing_modes_module.DatasetPreprocessing.run_analysis_ready
    )

    def counting_run_analysis_ready(self, **kwargs: object):
        calls.append(
            {
                "schema": self.schema,
                "comparisons": self.comparisons,
                "kwargs": kwargs,
            }
        )
        return original_run_analysis_ready(self, **kwargs)

    monkeypatch.setattr(
        preprocessing_modes_module.DatasetPreprocessing,
        "run_analysis_ready",
        counting_run_analysis_ready,
    )

    result = build_analysis_ready_dataset(
        total=total_path,
        phospho=phospho_path,
        source="analysis ready dataset builder",
    )

    assert isinstance(result, AnalysisReadyPhosphoDataset)
    assert len(calls) == 1
    assert calls[0]["schema"] == DatasetSchema()
    assert calls[0]["comparisons"] is None
    assert calls[0]["kwargs"]["source"] == "analysis ready dataset builder"


def test_analysis_ready_builder_phospho_only_reuses_core_processor_seam(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.preprocessing.modes as preprocessing_modes_module

    fixture_dir = (
        Path(__file__).resolve().parents[1] / "examples" / "data" / "simple_workflow"
    )
    phospho_path = fixture_dir / "phospho.tsv"
    calls: list[dict[str, object]] = []
    original_process_phospho_only = (
        preprocessing_modes_module.CoreProcessor.process_phospho_only
    )

    def counting_process_phospho_only(self, phospho_df: pd.DataFrame, *, config=None):
        calls.append(
            {
                "schema": self.schema,
                "comparisons": self.comparisons,
                "config": config,
                "rows": len(phospho_df),
            }
        )
        return original_process_phospho_only(self, phospho_df, config=config)

    monkeypatch.setattr(
        preprocessing_modes_module.CoreProcessor,
        "process_phospho_only",
        counting_process_phospho_only,
    )

    result = build_analysis_ready_dataset(
        phospho=phospho_path,
        phospho_only_source="analysis ready dataset builder (phospho only)",
    )

    assert isinstance(result, AnalysisReadyPhosphoDataset)
    assert len(calls) == 1
    assert calls[0]["schema"] == DatasetSchema()
    assert calls[0]["comparisons"] is None
    assert calls[0]["config"] is not None
    assert calls[0]["rows"] > 0
    assert result.provenance.source == "analysis ready dataset builder (phospho only)"


def test_build_analysis_ready_dataset_accepts_mixed_input_sources() -> None:
    fixture_dir = (
        Path(__file__).resolve().parents[1] / "examples" / "data" / "simple_workflow"
    )
    total_path = fixture_dir / "total.tsv"
    phospho_path = fixture_dir / "phospho.tsv"
    phospho_df = pd.read_csv(phospho_path, sep="\t")

    expected = build_analysis_ready_dataset(
        total=total_path,
        phospho=phospho_path,
        source="analysis ready dataset builder",
    )
    result = build_analysis_ready_dataset(
        total=total_path,
        phospho=phospho_df,
        source="analysis ready dataset builder",
    )

    pd.testing.assert_frame_equal(result.phospho_matrix, expected.phospho_matrix)
    pd.testing.assert_frame_equal(result.site_metadata, expected.site_metadata)
    pd.testing.assert_series_equal(result.site_sequences, expected.site_sequences)
    pd.testing.assert_frame_equal(
        result.phospho_corrected,
        expected.phospho_corrected,
    )
    assert result.provenance == expected.provenance
