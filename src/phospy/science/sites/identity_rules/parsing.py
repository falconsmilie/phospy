"""Shared parsing helpers for phosphosite identity validation."""

from __future__ import annotations

import math
from typing import TypeVar, cast

import numpy as np
import pandas as pd

from phospy.science.sites.identifiers import (
    ParsedSiteToken,
    canonicalize_site_series,
    try_parse_site_token,
)
from phospy.science.sites.site_keys import require_positive_integer_position

ErrorType = TypeVar("ErrorType", bound=Exception)
_SITE_POSITION_CANDIDATE_COLUMNS = ("position", "site_position")


def required_text_value(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    column_name: str,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    """Require one row metadata value to be a non-empty string."""

    if column_name not in site_metadata.columns:
        raise error_type(f"{field_name} is missing required columns: {column_name}")
    value = site_metadata.at[row_id, column_name]
    if not isinstance(value, str):
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a non-empty string"
        )
    token = value.strip()
    if token == "":
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a non-empty string"
        )
    return token


def optional_text_value(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    column_name: str,
    field_name: str,
    error_type: type[ErrorType],
) -> str | None:
    """Resolve one optional row metadata value as stripped text."""

    if column_name not in site_metadata.columns:
        return None
    value = site_metadata.at[row_id, column_name]
    if is_missing(value):
        return None
    if not isinstance(value, str):
        raise error_type(
            f"{field_name}[{row_id!r}].{column_name} must be a string when provided"
        )
    token = value.strip()
    if token == "":
        return None
    return token


def resolve_row_residue(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> str:
    """Resolve one row's phosphosite residue from explicit or site-token metadata."""

    explicit_residue = (
        optional_text_value(
            site_metadata=site_metadata,
            row_id=row_id,
            column_name="residue",
            field_name=field_name,
            error_type=error_type,
        )
        if "residue" in site_metadata.columns
        else None
    )
    parsed_site = parse_row_site_token(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    if explicit_residue is None:
        if parsed_site is None:
            raise error_type(
                f"{field_name}[{row_id!r}] requires residue metadata or strict site "
                "token parsing to derive ProteinScopedPhosphositeKey"
            )
        return parsed_site.residue
    token = explicit_residue.upper()
    if len(token) != 1 or token not in {"S", "T", "Y"}:
        raise error_type(
            f"{field_name}[{row_id!r}].residue must be one of 'S', 'T', or 'Y'"
        )
    if parsed_site is not None and parsed_site.residue != token:
        raise error_type(
            f"{field_name}[{row_id!r}] has inconsistent residue metadata; "
            f"site_token={parsed_site.residue!r}, residue_column={token!r}"
        )
    return token


def resolve_row_position(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> int:
    """Resolve one row's phosphosite position from explicit or site-token metadata."""

    explicit_position = resolve_explicit_position(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    parsed_site = parse_row_site_token(
        site_metadata=site_metadata,
        row_id=row_id,
        field_name=field_name,
        error_type=error_type,
    )
    if explicit_position is None:
        if parsed_site is None:
            raise error_type(
                f"{field_name}[{row_id!r}] requires position metadata or strict site "
                "token parsing to derive ProteinScopedPhosphositeKey"
            )
        return int(parsed_site.position)
    if parsed_site is not None and explicit_position != int(parsed_site.position):
        raise error_type(
            f"{field_name}[{row_id!r}] has inconsistent position metadata; "
            f"site_token={parsed_site.position}, position_column={explicit_position}"
        )
    return explicit_position


def resolve_explicit_position(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> int | None:
    """Resolve the first supported explicit site-position column."""

    for column_name in _SITE_POSITION_CANDIDATE_COLUMNS:
        if column_name not in site_metadata.columns:
            continue
        raw_value = site_metadata.at[row_id, column_name]
        return require_positive_integer_position(
            raw_value,
            field_name=f"{field_name}[{row_id!r}].{column_name}",
            error_type=error_type,
        )
    return None


def parse_row_site_token(
    *,
    site_metadata: pd.DataFrame,
    row_id: object,
    field_name: str,
    error_type: type[ErrorType],
) -> ParsedSiteToken | None:
    """Parse one row's strict ``site`` token when present."""

    if "site" not in site_metadata.columns:
        return None
    parsed = try_parse_site_token(site_metadata.at[row_id, "site"])
    if parsed is not None:
        return parsed
    raw_value = site_metadata.at[row_id, "site"]
    if is_missing(raw_value):
        return None
    if isinstance(raw_value, str) and raw_value.strip() == "":
        return None
    raise error_type(
        f"{field_name}[{row_id!r}].site must use strict 'S/T/Y<position>' tokens"
    )


def parse_site_token(
    value: object,
    *,
    field_name: str,
    error_type: type[ErrorType],
) -> ParsedSiteToken:
    """Require one strict ``S/T/Y<position>`` site token."""

    parsed = try_parse_site_token(value)
    if parsed is None:
        raise error_type(f"{field_name} must use strict 'S/T/Y<position>' tokens")
    return parsed


def is_missing(value: object) -> bool:
    """Return whether a scalar value should be treated as missing."""

    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value = cast(object, value)
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value = cast(object, value)
        return str(temporal_value) == "NaT"
    return False


def looks_like_display_site_index(index: pd.Index) -> bool:
    """Return whether an index looks like display ``GENE;SITE;`` labels."""

    values = index.tolist()
    if not values:
        return False
    if not all(isinstance(value, str) for value in values):
        return False
    if not any(";" in value for value in values):
        return False
    try:
        canonicalize_site_series(
            pd.Series(values, dtype="object"),
            field_name="analysis-ready row identity",
            error_type=ValueError,
        )
    except ValueError:
        return False
    return True


__all__ = [
    "is_missing",
    "looks_like_display_site_index",
    "optional_text_value",
    "parse_row_site_token",
    "parse_site_token",
    "required_text_value",
    "resolve_explicit_position",
    "resolve_row_position",
    "resolve_row_residue",
]
