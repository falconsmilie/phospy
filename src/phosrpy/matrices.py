from __future__ import annotations

from typing import Sequence

import pandas as pd


def build_site_matrix(
    df: pd.DataFrame,
    gene_p_site_col: str,
    sequence_col: str,
    value_cols: Sequence[str],
    gene_col_name: str = "gene",
    p_site_col_name: str = "p_site",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    work = df.copy()

    split_cols = work[gene_p_site_col].astype("string").str.split("_", n=1, expand=True)
    work[gene_col_name] = split_cols[0].astype("string")
    work[p_site_col_name] = split_cols[1].astype("string")
    work["site_id"] = (
        work[gene_col_name].str.upper() + ";" + work[p_site_col_name].str.upper() + ";"
    )

    base_cols = [
        col for col in ["gene_names", gene_col_name, p_site_col_name, "uid", sequence_col, "site_id"]
        if col in work.columns
    ]
    keep_cols: list[str] = [*base_cols, *list(value_cols)]
    phosr_input = work.loc[:, keep_cols].copy()

    phosr_input = phosr_input[phosr_input[sequence_col].notna()].copy()
    phosr_input = phosr_input[phosr_input.loc[:, list(value_cols)].notna().all(axis=1)].copy()

    phosr_input["__mean_signal"] = phosr_input.loc[:, list(value_cols)].mean(axis=1, skipna=True)
    idx = phosr_input.groupby("site_id")["__mean_signal"].idxmax()
    phosr_input = phosr_input.loc[idx].drop(columns="__mean_signal").copy().reset_index(drop=True)

    matrix = phosr_input.loc[:, ["site_id", *value_cols]].set_index("site_id")
    sequences = phosr_input.set_index("site_id")[sequence_col].copy()
    return phosr_input, matrix, sequences
