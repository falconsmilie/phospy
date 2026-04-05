from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd

from .constants import ComparisonSpec
from .dataset_schema import DatasetSchema
from .validation.errors import TableSchemaError

"""Internal preprocessing transformation primitives.

These helpers implement the shared mutation-oriented building blocks used by
both the public preprocessing façade and the internal preprocessing service
layer. They are intentionally kept out of the public helper module so the
service layer does not depend on public-module internals.
"""


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


def _require_numeric_series(
    values: pd.Series,
    *,
    column: str,
    context: str,
) -> pd.Series:
    try:
        numeric_values = pd.to_numeric(values, errors="raise")
    except (TypeError, ValueError) as exc:
        raise TableSchemaError(
            f"{context} requires numeric values in column '{column}'"
        ) from exc

    return pd.Series(numeric_values, index=values.index, copy=False)


def _filter_localized_sites_without_copy(
    df: pd.DataFrame,
    *,
    localization_col: str,
    threshold: float,
    localization_values: pd.Series | None = None,
) -> pd.DataFrame:
    values = localization_values
    if values is None:
        values = _require_numeric_series(
            df[localization_col],
            column=localization_col,
            context="filter_localized_sites()",
        )
    return df.loc[values >= threshold]


def _replace_sentinel_with_nan_in_place(
    df: pd.DataFrame,
    columns: Iterable[str],
    sentinel: float | int,
) -> pd.DataFrame:
    cols = list(columns)
    if not cols:
        return df

    numeric_block = df.loc[:, cols].astype(float)
    cleaned_block = numeric_block.mask(numeric_block == float(sentinel), np.nan)
    df[cols] = cleaned_block
    return df


def _filter_min_observed_without_copy(
    df: pd.DataFrame,
    columns: Sequence[str],
    min_observed: int,
) -> pd.DataFrame:
    mask = df.loc[:, list(columns)].notna().sum(axis=1) >= min_observed
    return df.loc[mask]


def _collapse_duplicate_genes_owned(
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

    df[gene_col] = df[gene_col].astype("string")
    if uppercase:
        df[gene_col] = df[gene_col].str.upper()
    ranked_cols = list(value_cols)
    df["__observed_count"] = df.loc[:, ranked_cols].notna().sum(axis=1)
    df["__mean_signal"] = df.loc[:, ranked_cols].mean(axis=1, skipna=True)
    df["__original_order"] = np.arange(len(df), dtype=int)

    ranked = df.sort_values(
        by=[gene_col, "__observed_count", "__mean_signal", "__original_order"],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    ranked = ranked.loc[ranked["__observed_count"] > 0]
    result = ranked.drop_duplicates(subset=[gene_col], keep="first").drop(
        columns=["__observed_count", "__mean_signal", "__original_order"]
    )

    return result.reset_index(drop=True)


def _add_pairwise_comparisons_in_place(
    df: pd.DataFrame,
    comparisons: Sequence[ComparisonSpec],
    group_to_corrected_col: dict[str, str] | None = None,
    output_prefix: str = "p_",
    schema: DatasetSchema | None = None,
) -> pd.DataFrame:
    if group_to_corrected_col is None:
        resolved_schema = schema or DatasetSchema()
        resolved_schema.validate_comparisons(
            comparisons,
            context="Pairwise comparison configuration",
        )
        group_to_corrected_col = dict(resolved_schema.group_to_corrected_col)
    else:
        resolved = tuple(comparisons)
        valid_groups = frozenset(group_to_corrected_col)
        seen: set[tuple[str, str]] = set()
        for left, right in resolved:
            if left not in valid_groups or right not in valid_groups:
                raise KeyError(f"Missing group mapping for comparison: {(left, right)}")
            pair = (left, right)
            if pair in seen:
                raise ValueError(f"Duplicate comparison pair: {left!r}, {right!r}")
            seen.add(pair)

    for left, right in comparisons:
        if left not in group_to_corrected_col or right not in group_to_corrected_col:
            raise KeyError(f"Missing group mapping for comparison: {(left, right)}")
        df[f"{output_prefix}{left}_{right}"] = (
            df[group_to_corrected_col[left]] - df[group_to_corrected_col[right]]
        )

    return df
