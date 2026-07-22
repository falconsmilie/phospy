"""Compatibility route for generic table schema infrastructure."""

from phospy.frames.table_schema import (
    TableSchema,
    ValidationErrorType,
    require_canonical_label_index,
)

__all__ = [
    "TableSchema",
    "ValidationErrorType",
    "require_canonical_label_index",
]
