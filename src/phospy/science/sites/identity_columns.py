"""Science-owned phosphosite identity column guards."""

from __future__ import annotations

from typing import TypeVar

import pandas as pd

from phospy.science.sites.identifiers import canonicalize_site_series
from phospy.science.sites.site_keys import decode_site_key, encode_site_key

ErrorType = TypeVar("ErrorType", bound=Exception)


def enforce_site_key_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "site_key",
) -> pd.Series:
    """Require one present site_key column with decodable encoded key values."""

    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    values = pd.Series(site_metadata[column_name], dtype="object")
    normalized: list[str] = []
    for row_id, raw_value in values.items():
        key = decode_site_key(
            raw_value,
            field_name=f"{field_name}.{column_name}[{row_id!r}]",
            error_type=error_type,
        )
        normalized.append(encode_site_key(key))
    return pd.Series(
        normalized,
        index=pd.Index(site_metadata.index),
        name=column_name,
        dtype="object",
    )


def enforce_display_id_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "display_id",
) -> pd.Series:
    """Require one present display_id column with recommended site identifiers."""

    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    return canonicalize_site_series(
        pd.Series(site_metadata[column_name], dtype="object"),
        field_name=f"{field_name}.{column_name}",
        error_type=error_type,
    )


__all__ = ["enforce_display_id_column", "enforce_site_key_column"]
