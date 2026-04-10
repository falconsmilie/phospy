from __future__ import annotations

import pandas as pd

from ...constants import GENE_P_SITE_COLUMN
from ..errors import TableSchemaError


def normalize_identifier_series(series: pd.Series) -> pd.Series:
    """Normalize identifier values for case/whitespace-insensitive joins."""

    return series.astype("string").str.strip().str.upper()


def require_splitable_gene_p_site(
    series: pd.Series,
    *,
    context: str,
    column_name: str = GENE_P_SITE_COLUMN,
) -> None:
    """Validate gene-site identifiers that must split on a single underscore."""

    normalized = series.astype("string")
    split_columns = normalized.str.split("_", n=1, expand=True)
    underscore_count = normalized.str.count("_")
    if split_columns.shape[1] < 2:
        invalid_mask = pd.Series(True, index=series.index)
    else:
        invalid_mask = (
            normalized.isna()
            | (underscore_count != 1)
            | split_columns[0].isna()
            | split_columns[1].isna()
            | (split_columns[0].str.strip().str.len() == 0)
            | (split_columns[1].str.strip().str.len() == 0)
        )
    if invalid_mask.any():
        sample_values = series.loc[invalid_mask].astype(str).unique()[:3]
        sample_preview = ", ".join(str(value) for value in sample_values)
        msg = (
            f"{context} contains malformed {column_name} values that cannot be split "
            f"into non-empty gene and site parts using a single underscore: "
            f"{sample_preview}"
        )
        raise TableSchemaError(msg)


__all__ = ["normalize_identifier_series", "require_splitable_gene_p_site"]
