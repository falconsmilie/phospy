"""Compatibility re-exports for differential table validation helpers."""

from phospy.science.tables._differential_validation import (
    require_boolean,
    require_column_name,
    require_differential_result_columns,
    require_na_position,
    require_non_negative_threshold,
    require_numeric_result_column,
    require_probability_threshold,
)

__all__ = [
    "require_boolean",
    "require_column_name",
    "require_differential_result_columns",
    "require_na_position",
    "require_non_negative_threshold",
    "require_numeric_result_column",
    "require_probability_threshold",
]
