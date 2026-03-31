from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy.preprocessing import (
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_min_observed,
    replace_sentinel_with_nan,
)
from phospy.validation.errors import InputCompatibilityError, TableSchemaError


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
