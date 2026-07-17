"""Compatibility imports for neutral DataFrame validation primitives."""

from __future__ import annotations

from phospy.frames.validation import (
    ValidationErrorType,
    format_label_examples,
    require_aligned_dataframe_shape,
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_no_duplicate_labels,
    require_non_empty_dataframe,
    require_non_empty_index_intersection,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_string_index,
    require_unique_columns,
    require_unique_index,
    require_unique_row_pairs,
    summarise_column_mismatch,
    summarise_index_mismatch,
)

__all__ = [
    "ValidationErrorType",
    "format_label_examples",
    "require_aligned_dataframe_shape",
    "require_canonical_string_column",
    "require_columns",
    "require_dataframe",
    "require_exact_index_match",
    "require_finite_numeric_dataframe",
    "require_no_duplicate_labels",
    "require_non_empty_dataframe",
    "require_non_empty_index_intersection",
    "require_non_empty_string_column",
    "require_numeric_dataframe",
    "require_string_index",
    "require_unique_columns",
    "require_unique_index",
    "require_unique_row_pairs",
    "summarise_column_mismatch",
    "summarise_index_mismatch",
]
