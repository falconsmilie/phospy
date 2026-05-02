"""Shared DataFrame-level validation helpers."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.site_ids import (
    canonicalize_site_index,
    canonicalize_site_series,
    parse_canonical_site_identifier,
)

ValidationErrorType = type[PhosPyValidationError]


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
    if not allow_empty:
        require_non_empty_dataframe(
            value,
            field_name=field_name,
            error_type=error_type,
        )
    return value


def require_non_empty_dataframe(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require a DataFrame to contain at least one row and one column."""

    if value.empty:
        raise error_type(
            f"{field_name} must be non-empty; "
            f"rows={int(value.shape[0])}, columns={int(value.shape[1])}"
        )
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
        column_count = int(len(non_numeric_columns))
        raise error_type(
            f"{field_name} must contain only numeric columns; non-numeric columns: "
            f"{joined_columns}; non_numeric_column_count={column_count}"
        )
    return value


def require_finite_numeric_dataframe(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
    allow_missing: bool = False,
) -> pd.DataFrame:
    """Require one numeric DataFrame to satisfy finite-value constraints."""

    missing_mask = value.isna()
    if not allow_missing and missing_mask.to_numpy().any():
        missing_count = int(missing_mask.to_numpy().sum())
        raise error_type(
            f"{field_name} must not contain missing values; "
            f"{_invalid_location_preview(missing_mask)}; missing_entries={missing_count}"
        )
    infinite_mask = value.isin([float("inf"), float("-inf")])
    if infinite_mask.to_numpy().any():
        infinite_count = int(infinite_mask.to_numpy().sum())
        raise error_type(
            f"{field_name} must contain finite numeric values; "
            f"{_invalid_location_preview(infinite_mask)}; infinite_entries={infinite_count}"
        )
    return value


def require_unique_index(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require unique index labels in a DataFrame."""

    require_no_duplicate_labels(
        value.index,
        field_name=f"{field_name}.index",
        error_type=error_type,
    )
    return value


def require_unique_columns(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require unique column labels in a DataFrame."""

    require_no_duplicate_labels(
        value.columns,
        field_name=f"{field_name}.columns",
        error_type=error_type,
    )
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
        raise error_type(
            f"{field_name} is missing required columns: {joined}; "
            f"missing_column_count={int(len(missing))}, "
            f"column_count={int(value.shape[1])}"
        )
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
        raise error_type(
            f"{left_name} must exactly match {right_name}; "
            f"{left_name}_count={int(left.size)}, {right_name}_count={int(right.size)}"
        )


def require_non_empty_index_intersection(
    *,
    left: pd.Index,
    right: pd.Index,
    left_name: str,
    right_name: str,
    error_type: ValidationErrorType,
) -> pd.Index:
    """Require two indexes to share at least one label and return the overlap."""

    shared = left.intersection(right)
    if shared.size > 0:
        return shared
    raise error_type(
        f"{left_name} and {right_name} must share at least one label; "
        f"{left_name}_count={int(left.size)}, {right_name}_count={int(right.size)}, "
        "shared_count=0"
    )


def require_no_duplicate_labels(
    labels: pd.Index,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.Index:
    """Require one index-like label sequence to be unique."""

    duplicated = labels[labels.duplicated(keep=False)]
    if duplicated.size == 0:
        return labels
    duplicate_labels = duplicated.unique()
    preview = ", ".join(repr(label) for label in duplicate_labels.tolist()[:5])
    suffix = "" if duplicate_labels.size <= 5 else " ..."
    raise error_type(
        f"{field_name} must be unique; duplicate_count={int(duplicated.size)}, "
        f"duplicate_labels={preview}{suffix}"
    )


def require_string_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.Index:
    """Require one index to contain only string labels and no missing labels."""

    values = index.tolist()
    missing = [value for value in values if _is_missing_site_identifier(value)]
    if missing:
        raise error_type(
            f"{field_name} must not contain missing labels; "
            f"missing_label_count={int(len(missing))}"
        )
    non_strings = [value for value in values if not isinstance(value, str)]
    if non_strings:
        preview = ", ".join(repr(value) for value in non_strings[:5])
        suffix = "" if len(non_strings) <= 5 else " ..."
        raise error_type(
            f"{field_name} must contain only string labels; "
            f"non_string_label_count={int(len(non_strings))}, "
            f"examples={preview}{suffix}"
        )
    return index


def require_aligned_dataframe_shape(
    *,
    left: pd.DataFrame,
    right: pd.DataFrame,
    left_name: str,
    right_name: str,
    error_type: ValidationErrorType,
) -> None:
    """Require two DataFrames to have identical row/column counts."""

    if left.shape == right.shape:
        return
    left_rows, left_columns = left.shape
    right_rows, right_columns = right.shape
    raise error_type(
        f"{left_name} shape must align with {right_name}; "
        f"{left_name}_rows={int(left_rows)}, {left_name}_columns={int(left_columns)}, "
        f"{right_name}_rows={int(right_rows)}, {right_name}_columns={int(right_columns)}"
    )


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
    strict_supported_format: bool = True,
) -> pd.Index:
    """Require one site index to already be strict/canonical."""

    if not strict_supported_format:
        _require_stripped_site_identifiers(
            index.tolist(),
            field_name=field_name,
            error_type=error_type,
        )
        return index
    canonical = canonicalize_site_index(
        index,
        field_name=field_name,
        error_type=error_type,
        require_unique=True,
        index_name=(str(index.name) if index.name is not None else None),
    )
    if canonical.tolist() != index.tolist():
        raise error_type(
            f"{field_name} must contain canonical site identifiers in "
            "'GENE;SITE;' format"
        )
    if not index.is_unique:
        require_no_duplicate_labels(
            index,
            field_name=field_name,
            error_type=error_type,
        )
    return index


def require_canonical_site_series(
    series: pd.Series,
    *,
    field_name: str,
    error_type: ValidationErrorType,
    strict_supported_format: bool = True,
) -> pd.Series:
    """Require one site-id series to already be strict/canonical."""

    if not strict_supported_format:
        _require_stripped_site_identifiers(
            series.tolist(),
            field_name=field_name,
            error_type=error_type,
        )
        return series
    canonical = canonicalize_site_series(
        series,
        field_name=field_name,
        error_type=error_type,
    )
    if canonical.tolist() != series.tolist():
        raise error_type(
            f"{field_name} must contain canonical site identifiers in "
            "'GENE;SITE;' format"
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
        parsed = _parse_site_identity(
            site_id,
            field_name=site_index_field_name,
            error_type=error_type,
        )
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


def _parse_site_identity(
    site_id: object,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> tuple[str, str] | None:
    if not isinstance(site_id, str):
        return None
    try:
        return parse_canonical_site_identifier(
            site_id,
            field_name=field_name,
            error_type=error_type,
        )
    except error_type:
        return None


def _require_stripped_site_identifiers(
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
        for raw_value, stripped_value in zip(raw_values, stripped_values, strict=False)
    ):
        raise error_type(
            f"{field_name} must contain canonical site identifiers (non-empty stripped strings)"
        )


def _is_missing_site_identifier(value: object) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _invalid_location_preview(mask: pd.DataFrame, *, max_items: int = 3) -> str:
    locations = np.argwhere(mask.to_numpy())
    count = int(locations.shape[0])
    preview = [
        f"({mask.index[row_idx]!r}, {mask.columns[col_idx]!r})"
        for row_idx, col_idx in locations[:max_items]
    ]
    suffix = "" if count <= max_items else f", +{count - max_items} more"
    return f"found at {', '.join(preview)}{suffix}"


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
