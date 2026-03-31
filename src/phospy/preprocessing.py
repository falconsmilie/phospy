from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd

from .constants import ComparisonSpec
from .validation.errors import InputCompatibilityError, TableSchemaError
from .validation.normalization import normalize_identifier_series


def _require_columns(
    df: pd.DataFrame,
    *,
    required_columns: Sequence[str],
    context: str,
) -> None:
    missing_columns = [
        column for column in required_columns if column not in df.columns
    ]
    if missing_columns:
        joined_columns = ", ".join(missing_columns)
        msg = f"{context} is missing required columns: {joined_columns}"
        raise TableSchemaError(msg)


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
    """Collapse duplicate gene rows using an explicit ranking policy.

    Gene identifiers are converted to pandas ``string`` values before ranking.
    When ``uppercase`` is ``True``, identifiers are normalised to uppercase
    before duplicate grouping so mixed-case duplicates collapse together.

    Rows are ranked within each gene group by:
    1. highest observed-value count across ``value_cols``
    2. highest mean signal across ``value_cols``
    3. earliest original row order as a stable tie-breaker

    Rows with zero observed values across ``value_cols`` are dropped before the
    winning row is selected so all-missing groups do not survive deduplication.
    The top-ranked remaining row for each gene is retained.
    """
    _require_columns(
        df,
        required_columns=[gene_col, *value_cols],
        context="collapse_duplicate_genes() input",
    )

    work = df.copy()
    work[gene_col] = work[gene_col].astype("string")
    if uppercase:
        work[gene_col] = work[gene_col].str.upper()
    ranked_cols = list(value_cols)
    work["__observed_count"] = work.loc[:, ranked_cols].notna().sum(axis=1)
    work["__mean_signal"] = work.loc[:, ranked_cols].mean(axis=1, skipna=True)
    work["__original_order"] = np.arange(len(work), dtype=int)

    ranked = work.sort_values(
        by=[gene_col, "__observed_count", "__mean_signal", "__original_order"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    ranked = ranked.loc[ranked["__observed_count"] > 0]
    result = ranked.drop_duplicates(subset=[gene_col], keep="first").drop(
        columns=["__observed_count", "__mean_signal", "__original_order"]
    )

    return result.reset_index(drop=True)


def correct_phospho_to_protein(
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    output_prefix: str = "phospho_corrected_",
) -> pd.DataFrame:
    if len(phospho_cols) != len(protein_cols):
        raise ValueError("phospho_cols and protein_cols must have the same length")

    resolved_corrected_cols = (
        list(corrected_cols)
        if corrected_cols is not None
        else [f"{output_prefix}{idx}" for idx in range(1, len(phospho_cols) + 1)]
    )
    if len(resolved_corrected_cols) != len(phospho_cols):
        raise ValueError(
            "corrected_cols must have the same length as phospho_cols and protein_cols"
        )

    phospho_join_col = "__phospy_normalized_phospho_gene_key"
    total_join_col = "__phospy_normalized_total_gene_key"

    phospho_work = df_phospho.copy()
    total_work = df_total.copy()
    phospho_work[phospho_join_col] = normalize_identifier_series(
        phospho_work[phospho_gene_col]
    )
    total_work[total_join_col] = normalize_identifier_series(total_work[total_gene_col])

    if total_work[total_join_col].duplicated().any():
        msg = (
            f"{total_gene_col} must be unique before protein correction to avoid "
            "duplicating phosphosite rows during the merge"
        )
        raise InputCompatibilityError(msg)

    merged = phospho_work.merge(
        total_work[[total_join_col, total_gene_col, *protein_cols]],
        left_on=phospho_join_col,
        right_on=total_join_col,
        how="inner",
    )

    drop_columns: list[str] = [phospho_join_col, total_join_col]
    if total_gene_col != phospho_gene_col and total_gene_col in merged.columns:
        drop_columns.append(total_gene_col)
    merged = merged.drop(columns=drop_columns, errors="ignore")

    for corrected_col, p_col, t_col in zip(
        resolved_corrected_cols,
        phospho_cols,
        protein_cols,
        strict=True,
    ):
        merged[corrected_col] = merged[p_col] - merged[t_col]

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
            raise KeyError(f"Missing group mapping for comparison: {(left, right)}")
        result[f"{output_prefix}{left}_{right}"] = (
            result[group_to_corrected_col[left]] - result[group_to_corrected_col[right]]
        )

    return result
