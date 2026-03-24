from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd

from .constants import ComparisonSpec
from .validation.errors import InputCompatibilityError, TableSchemaError


def replace_sentinel_with_nan(
    df: pd.DataFrame,
    columns: Iterable[str],
    sentinel: float | int,
) -> pd.DataFrame:
    result = df.copy()
    cols = list(columns)
    for col in cols:
        result[col] = result[col].astype(float).replace(sentinel, np.nan)
    return result


def filter_min_observed(
    df: pd.DataFrame,
    columns: Sequence[str],
    min_observed: int,
) -> pd.DataFrame:
    mask = df.loc[:, list(columns)].notna().sum(axis=1) >= min_observed
    return df.loc[mask].copy()


def collapse_duplicate_genes(
    df: pd.DataFrame,
    gene_col: str,
    value_cols: Sequence[str],
    uppercase: bool = True,
) -> pd.DataFrame:
    if gene_col not in df.columns:
        msg = f"Missing gene column: {gene_col}"
        raise TableSchemaError(msg)

    work = df.copy()
    work[gene_col] = work[gene_col].astype("string")
    work["__mean_signal"] = work.loc[:, list(value_cols)].mean(axis=1, skipna=True)
    idx = work.groupby(gene_col)["__mean_signal"].idxmax()
    result = work.loc[idx].drop(columns="__mean_signal").copy()

    if uppercase:
        result[gene_col] = result[gene_col].str.upper()

    return result.reset_index(drop=True)


def correct_phospho_to_protein(
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    output_prefix: str = "phospho_corrected_",
) -> pd.DataFrame:
    if len(phospho_cols) != len(protein_cols):
        msg = "phospho_cols and protein_cols must have the same length"
        raise InputCompatibilityError(msg)

    merged = df_phospho.merge(
        df_total[[total_gene_col, *protein_cols]],
        left_on=phospho_gene_col,
        right_on=total_gene_col,
        how="inner",
    )

    if total_gene_col != phospho_gene_col and total_gene_col in merged.columns:
        merged = merged.drop(columns=[total_gene_col])

    for idx, (p_col, t_col) in enumerate(
        zip(phospho_cols, protein_cols, strict=True),
        start=1,
    ):
        merged[f"{output_prefix}{idx}"] = merged[p_col] - merged[t_col]

    return merged


def add_pairwise_comparisons(
    df: pd.DataFrame,
    comparisons: Sequence[ComparisonSpec],
    group_to_corrected_col: dict[str, str] | None = None,
    output_prefix: str = "p_",
) -> pd.DataFrame:
    result = df.copy()
    if group_to_corrected_col is None:
        group_to_corrected_col = {
            f"group{i}": f"phospho_corrected_{i}" for i in range(1, 7)
        }

    for left, right in comparisons:
        if left not in group_to_corrected_col or right not in group_to_corrected_col:
            msg = f"Missing group mapping for comparison: {(left, right)}"
            raise InputCompatibilityError(msg)
        result[f"{output_prefix}{left}_{right}"] = (
            result[group_to_corrected_col[left]] - result[group_to_corrected_col[right]]
        )

    return result
