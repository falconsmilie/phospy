"""Shared DataFrame-level validation helpers."""

from __future__ import annotations

import numbers
from collections.abc import Iterable

import numpy as np
import pandas as pd

from phospy.errors.base import PhosPyError

ValidationErrorType = type[PhosPyError]
_ALIGNMENT_EXAMPLE_LIMIT = 5


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

    cell_values = value.to_numpy(dtype="object")
    missing_array = np.zeros(cell_values.shape, dtype=bool)
    infinite_array = np.zeros(cell_values.shape, dtype=bool)
    row_count = int(cell_values.shape[0])
    column_count = int(cell_values.shape[1])
    for row_idx in range(row_count):
        for col_idx in range(column_count):
            cell_value = cell_values[row_idx, col_idx]
            if _is_missing_label(cell_value):
                missing_array[row_idx, col_idx] = True
                continue
            if isinstance(cell_value, numbers.Real) and not np.isfinite(
                float(cell_value)
            ):
                infinite_array[row_idx, col_idx] = True

    missing_mask = pd.DataFrame(
        missing_array,
        index=value.index,
        columns=value.columns,
    )
    if not allow_missing and bool(missing_array.any()):
        missing_count = int(missing_array.sum())
        raise error_type(
            f"{field_name} must not contain missing values; "
            f"{_invalid_location_preview(missing_mask)}; missing_entries={missing_count}"
        )

    infinite_mask = pd.DataFrame(
        infinite_array,
        index=value.index,
        columns=value.columns,
    )
    if bool(infinite_array.any()):
        infinite_count = int(infinite_array.sum())
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
        summariser = (
            summarise_column_mismatch
            if left_name.endswith(".columns") and right_name.endswith(".columns")
            else summarise_index_mismatch
        )
        detail = summariser(
            left=left,
            right=right,
            left_name=left_name,
            right_name=right_name,
            max_examples=_ALIGNMENT_EXAMPLE_LIMIT,
        )
        raise error_type(
            f"{left_name} must exactly match {right_name}; "
            f"expected_length={int(right.size)}, actual_length={int(left.size)}; "
            f"{left_name}_count={int(left.size)}, {right_name}_count={int(right.size)}; "
            f"{detail}"
        )


def summarise_index_mismatch(
    *,
    left: pd.Index,
    right: pd.Index,
    left_name: str,
    right_name: str,
    max_examples: int = _ALIGNMENT_EXAMPLE_LIMIT,
) -> str:
    """Summarise one index alignment mismatch without dumping full label sets."""

    left_only, right_only = _labels_only_in_each_side(left=left, right=right)
    first_mismatch = _first_positional_mismatch(left=left, right=right)
    labels_match_set_order_differs = (
        not left_only and not right_only and not left.equals(right)
    )
    mismatch_text = (
        f"position {first_mismatch[0]}, "
        f"{left_name}={_format_label(first_mismatch[1])}, "
        f"{right_name}={_format_label(first_mismatch[2])}"
        if first_mismatch is not None
        else "none"
    )
    return "; ".join(
        (
            f"Only in {left_name}: {format_label_examples(left_only, max_examples=max_examples)}",
            f"Only in {right_name}: {format_label_examples(right_only, max_examples=max_examples)}",
            f"First positional mismatch: {mismatch_text}",
            (
                "Labels match as a set but order differs: "
                f"{str(labels_match_set_order_differs).lower()}"
            ),
        )
    )


def summarise_column_mismatch(
    *,
    left: pd.Index,
    right: pd.Index,
    left_name: str,
    right_name: str,
    max_examples: int = _ALIGNMENT_EXAMPLE_LIMIT,
) -> str:
    """Summarise one column-label mismatch between two aligned tables."""

    return summarise_index_mismatch(
        left=left,
        right=right,
        left_name=left_name,
        right_name=right_name,
        max_examples=max_examples,
    )


def format_label_examples(
    labels: Iterable[object],
    *,
    max_examples: int = _ALIGNMENT_EXAMPLE_LIMIT,
) -> str:
    """Format capped label examples for mismatch diagnostics."""

    values = list(labels)
    if not values:
        return "(none)"
    preview_count = max(int(max_examples), 1)
    preview = ", ".join(_format_label(value) for value in values[:preview_count])
    remaining = len(values) - preview_count
    if remaining > 0:
        return f"{preview}, +{remaining} more"
    return preview


def require_non_empty_index_intersection(
    *,
    left: pd.Index,
    right: pd.Index,
    left_name: str,
    right_name: str,
    error_type: ValidationErrorType,
) -> pd.Index:
    """Require two indexes to share at least one label and return the overlap."""

    left_labels = left.tolist()
    right_labels = right.tolist()
    try:
        right_lookup = set(right_labels)
        shared_labels = [label for label in left_labels if label in right_lookup]
    except TypeError:
        shared_labels = [label for label in left_labels if label in right_labels]
    shared = pd.Index(shared_labels)
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

    label_values = labels.tolist()
    duplicate_counts: dict[object, int] = {}
    for label in label_values:
        duplicate_counts[label] = duplicate_counts.get(label, 0) + 1
    duplicated_values = [label for label in label_values if duplicate_counts[label] > 1]
    if not duplicated_values:
        return labels
    duplicate_labels = list(dict.fromkeys(duplicated_values))
    preview = ", ".join(repr(label) for label in duplicate_labels[:5])
    suffix = "" if len(duplicate_labels) <= 5 else " ..."
    raise error_type(
        f"{field_name} must be unique; duplicate_count={int(len(duplicated_values))}, "
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
    missing = [value for value in values if _is_missing_label(value)]
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


def _is_missing_label(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _labels_only_in_each_side(
    *, left: pd.Index, right: pd.Index
) -> tuple[list[object], list[object]]:
    left_values = left.tolist()
    right_values = right.tolist()
    try:
        right_lookup = set(right_values)
        left_lookup = set(left_values)
        left_only = [label for label in left_values if label not in right_lookup]
        right_only = [label for label in right_values if label not in left_lookup]
    except TypeError:
        left_only = [label for label in left_values if label not in right_values]
        right_only = [label for label in right_values if label not in left_values]
    return left_only, right_only


def _first_positional_mismatch(
    *,
    left: pd.Index,
    right: pd.Index,
) -> tuple[int, object, object] | None:
    shared_length = min(int(left.size), int(right.size))
    for position in range(shared_length):
        left_label = left[position]
        right_label = right[position]
        if left_label != right_label:
            return position, left_label, right_label
    if left.size == right.size:
        return None
    missing = "<missing>"
    if left.size > right.size:
        return shared_length, left[shared_length], missing
    return shared_length, missing, right[shared_length]


def _format_label(label: object) -> str:
    return repr(label)


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

    left_column, right_column = column_names
    left_values = value[left_column].tolist()
    right_values = value[right_column].tolist()
    pair_counts: dict[tuple[object, object], int] = {}
    for pair in zip(left_values, right_values, strict=True):
        pair_counts[pair] = pair_counts.get(pair, 0) + 1
    duplicate_pairs = [pair for pair, count in pair_counts.items() if count > 1]
    if not duplicate_pairs:
        return value
    preview_pairs = duplicate_pairs
    preview = ", ".join(repr(pair) for pair in preview_pairs[:5])
    suffix = "" if len(preview_pairs) <= 5 else " ..."
    left, right = column_names
    raise error_type(
        f"{field_name} contains duplicate ({left}, {right}) pairs: {preview}{suffix}"
    )
