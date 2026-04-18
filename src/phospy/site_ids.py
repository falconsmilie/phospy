"""Shared site-identifier canonicalization helpers."""

from __future__ import annotations

from typing import TypeVar

import pandas as pd

ErrorType = TypeVar("ErrorType", bound=Exception)


def canonicalize_site_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[ErrorType],
    require_unique: bool = True,
    index_name: str | None = None,
) -> pd.Index:
    """Canonicalize one site-identifier index to stripped string labels."""

    canonical = pd.Index(
        [
            canonicalize_site_identifier(
                value,
                field_name=field_name,
                error_type=error_type,
            )
            for value in index.tolist()
        ],
        name=index.name if index_name is None else index_name,
    )
    if require_unique and not canonical.is_unique:
        duplicates = list(dict.fromkeys(canonical[canonical.duplicated(keep=False)]))
        preview = ", ".join(duplicates[:5])
        suffix = "" if len(duplicates) <= 5 else " ..."
        raise error_type(
            f"{field_name} contains duplicate site identifiers after "
            f"canonicalization: {preview}{suffix}"
        )
    return canonical


def canonicalize_site_series(
    series: pd.Series,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> pd.Series:
    """Canonicalize one site-identifier series to stripped string labels."""

    return pd.Series(
        [
            canonicalize_site_identifier(
                value,
                field_name=field_name,
                error_type=error_type,
            )
            for value in series.tolist()
        ],
        index=series.index.copy(),
        name=series.name,
        dtype="object",
    )


def canonicalize_site_identifier(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    """Canonicalize one site identifier to a non-empty stripped string label."""

    try:
        is_missing = bool(pd.isna(value))
    except (TypeError, ValueError):
        is_missing = False
    if is_missing:
        raise error_type(f"{field_name} must not contain missing site identifiers")
    canonical = str(value).strip()
    if canonical == "":
        raise error_type(f"{field_name} must contain non-empty site identifiers")
    return canonical
