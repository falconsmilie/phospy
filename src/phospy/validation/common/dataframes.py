"""Shared DataFrame-level validation helpers."""

from __future__ import annotations

import re
from collections.abc import Iterable

import pandas as pd

from phospy.errors.validation import PhosPyValidationError

ValidationErrorType = type[PhosPyValidationError]
_SITE_IDENTITY_PATTERN = re.compile(
    r"^\s*(?P<gene_symbol>[^;]+?)\s*;\s*(?P<site>[^;]+?)\s*;\s*$"
)


def require_dataframe(
    value: object,
    *,
    field_name: str,
    allow_empty: bool,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require a pandas DataFrame and optionally reject empty frames."""

    if not isinstance(value, pd.DataFrame):
        raise error_type(f"{field_name} must be a pandas DataFrame")
    if not allow_empty and value.empty:
        raise error_type(f"{field_name} must be non-empty")
    return value


def require_numeric_dataframe(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require all columns in a DataFrame to be numeric and non-boolean."""

    boolean_columns = [
        str(column)
        for column in value.columns
        if pd.api.types.is_bool_dtype(value[column])
    ]
    if boolean_columns:
        joined_columns = ", ".join(boolean_columns)
        raise error_type(
            f"{field_name} must contain only scientific numeric columns; "
            f"boolean columns are invalid: {joined_columns}"
        )

    non_numeric_columns = [
        str(column)
        for column in value.columns
        if not pd.api.types.is_numeric_dtype(value[column])
    ]
    if non_numeric_columns:
        joined_columns = ", ".join(non_numeric_columns)
        raise error_type(
            f"{field_name} must contain only numeric columns; non-numeric columns: "
            f"{joined_columns}"
        )
    return value


def require_unique_index(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require unique index labels in a DataFrame."""

    if not value.index.is_unique:
        raise error_type(f"{field_name}.index must be unique")
    return value


def require_unique_columns(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require unique column labels in a DataFrame."""

    if not value.columns.is_unique:
        raise error_type(f"{field_name}.columns must be unique")
    return value


def require_columns(
    value: pd.DataFrame,
    *,
    field_name: str,
    required_columns: Iterable[str],
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require a DataFrame to include the given columns."""

    missing = [column for column in required_columns if column not in value.columns]
    if missing:
        joined = ", ".join(missing)
        raise error_type(f"{field_name} is missing required columns: {joined}")
    return value


def require_exact_index_match(
    *,
    left: pd.Index,
    right: pd.Index,
    left_name: str,
    right_name: str,
    error_type: ValidationErrorType,
) -> None:
    """Require two indexes to be exactly equal (labels and order)."""

    if not left.equals(right):
        raise error_type(f"{left_name} must exactly match {right_name}")


def require_non_empty_string_column(
    value: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require all values in a DataFrame column to be non-empty strings."""

    column_values = value[column_name]
    if column_values.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    non_string_or_blank = [
        idx
        for idx, raw_value in column_values.items()
        if not isinstance(raw_value, str) or not raw_value.strip()
    ]
    if non_string_or_blank:
        raise error_type(
            f"{field_name}.{column_name} must contain non-empty string values"
        )
    return value


def require_canonical_string_column(
    value: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require one string column to be stripped, non-empty, and non-missing."""

    column_values = value[column_name]
    if column_values.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    invalid = [
        idx
        for idx, raw_value in column_values.items()
        if not isinstance(raw_value, str)
        or raw_value == ""
        or raw_value != raw_value.strip()
    ]
    if invalid:
        raise error_type(
            f"{field_name}.{column_name} must contain canonical non-empty string values"
        )
    return value


def require_canonical_site_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.Index:
    """Require one site index to already be strict/canonical."""

    _require_strict_site_identifiers(
        index.tolist(),
        field_name=field_name,
        error_type=error_type,
    )
    return index


def require_canonical_site_series(
    series: pd.Series,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.Series:
    """Require one site-id series to already be strict/canonical."""

    _require_strict_site_identifiers(
        series.tolist(),
        field_name=field_name,
        error_type=error_type,
    )
    return series


def require_site_identity_coherence(
    *,
    site_index: pd.Index,
    site_metadata: pd.DataFrame,
    site_index_field_name: str,
    site_metadata_field_name: str,
    gene_symbol_column: str = "gene_symbol",
    site_column: str = "site",
    error_type: ValidationErrorType,
    error_preview_limit: int = 5,
) -> None:
    """Require canonical site IDs to agree with metadata gene/site columns."""

    unparseable_site_ids: list[str] = []
    mismatched_rows: list[str] = []

    for site_id in site_index:
        parsed = _parse_site_identity(site_id)
        if parsed is None:
            unparseable_site_ids.append(str(site_id))
            continue

        expected_gene_symbol, expected_site = parsed
        observed_gene_symbol = site_metadata.at[site_id, gene_symbol_column]
        observed_site = site_metadata.at[site_id, site_column]
        if (
            observed_gene_symbol != expected_gene_symbol
            or observed_site != expected_site
        ):
            mismatched_rows.append(
                f"{site_id} expected(gene_symbol={expected_gene_symbol!r}, "
                f"site={expected_site!r}) observed(gene_symbol={observed_gene_symbol!r}, "
                f"site={observed_site!r})"
            )

    if not unparseable_site_ids and not mismatched_rows:
        return

    details: list[str] = []
    if unparseable_site_ids:
        preview = ", ".join(
            repr(site_id) for site_id in unparseable_site_ids[:error_preview_limit]
        )
        suffix = "" if len(unparseable_site_ids) <= error_preview_limit else " ..."
        details.append(
            f"unparseable site IDs for '<gene_symbol>;<site>;': {preview}{suffix}"
        )
    if mismatched_rows:
        preview = "; ".join(mismatched_rows[:error_preview_limit])
        suffix = "" if len(mismatched_rows) <= error_preview_limit else " ..."
        details.append(f"mismatched rows: {preview}{suffix}")

    joined_details = "; ".join(details)
    raise error_type(
        "dataset site-identity coherence failed: "
        f"{site_index_field_name} canonical site IDs must agree with "
        f"{site_metadata_field_name}.{gene_symbol_column} and "
        f"{site_metadata_field_name}.{site_column}; {joined_details}"
    )


def _parse_site_identity(site_id: object) -> tuple[str, str] | None:
    if not isinstance(site_id, str):
        return None

    match = _SITE_IDENTITY_PATTERN.fullmatch(site_id)
    if match is None:
        return None

    gene_symbol = match.group("gene_symbol").strip()
    site = match.group("site").strip()
    if gene_symbol == "" or site == "":
        return None
    return gene_symbol, site


def _require_strict_site_identifiers(
    values: list[object],
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> None:
    if any(_is_missing_site_identifier(value) for value in values):
        raise error_type(f"{field_name} must not contain missing site identifiers")
    if not all(isinstance(value, str) for value in values):
        raise error_type(
            f"{field_name} must contain canonical site identifiers (non-empty stripped strings)"
        )
    raw_values = [value for value in values if isinstance(value, str)]
    stripped_values = [value.strip() for value in raw_values]
    if any(value == "" for value in stripped_values):
        raise error_type(f"{field_name} must contain non-empty site identifiers")

    collisions: dict[str, set[str]] = {}
    for raw_value, stripped_value in zip(raw_values, stripped_values, strict=False):
        collisions.setdefault(stripped_value, set()).add(raw_value)
    colliding = [value for value, raw_set in collisions.items() if len(raw_set) > 1]
    if colliding:
        preview = ", ".join(colliding[:5])
        suffix = "" if len(colliding) <= 5 else " ..."
        raise error_type(
            f"{field_name} contains colliding site identifiers when stripped: "
            f"{preview}{suffix}"
        )
    if any(
        raw_value != stripped_value
        for raw_value, stripped_value in zip(
            raw_values,
            stripped_values,
            strict=False,
        )
    ):
        raise error_type(
            f"{field_name} must contain canonical site identifiers (non-empty stripped strings)"
        )


def _is_missing_site_identifier(value: object) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def require_unique_row_pairs(
    value: pd.DataFrame,
    *,
    field_name: str,
    column_names: tuple[str, str],
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require unique row pairs for one two-column key."""

    duplicated = value.duplicated(subset=list(column_names), keep=False)
    if not bool(duplicated.any()):
        return value
    duplicate_pairs = (
        value.loc[duplicated, list(column_names)]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    preview_pairs = list(duplicate_pairs)
    preview = ", ".join(repr(pair) for pair in preview_pairs[:5])
    suffix = "" if len(preview_pairs) <= 5 else " ..."
    left, right = column_names
    raise error_type(
        f"{field_name} contains duplicate ({left}, {right}) pairs: {preview}{suffix}"
    )
