"""Compatibility re-exports for activity table schemas."""

from phospy.science.tables.activity import (
    ActivityCountMatrix,
    ActivityCountSeries,
    ActivityMatrix,
    ActivityStatisticsTable,
    ActivityTargetTable,
    SeriesSchema,
)
from phospy.science.tables.activity import (
    _column_series as _column_series,
)
from phospy.science.tables.activity import (
    _require_integer_compatible_column as _require_integer_compatible_column,
)
from phospy.science.tables.activity import (
    _require_numeric_column as _require_numeric_column,
)
from phospy.science.tables.activity import (
    _require_numeric_column_allowing_missing as _require_numeric_column_allowing_missing,
)

__all__ = [
    "ActivityCountMatrix",
    "ActivityCountSeries",
    "ActivityMatrix",
    "ActivityStatisticsTable",
    "ActivityTargetTable",
    "SeriesSchema",
]
