from __future__ import annotations

import numpy as np
import pandas as pd

from phosrpy.preprocessing import (
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_min_observed,
    replace_sentinel_with_nan,
)


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
