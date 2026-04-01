from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.dataset_schema import DatasetSchema
from phospy.preprocessing import (
    CoverageFilterResult,
    LocalizationFilterResult,
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_localized_sites,
    filter_min_observed,
    filter_sites_by_coverage,
    replace_sentinel_with_nan,
)
from phospy.preprocessing_services import (
    PhosphoPreprocessor,
    ProteinCorrectionService,
    TotalPreprocessor,
)
from phospy.validation.errors import (
    InputCompatibilityError,
    PhospyValidationError,
    TableSchemaError,
)


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
