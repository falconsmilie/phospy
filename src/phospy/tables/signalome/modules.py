from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.frames.validation import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.tables.base import TableSchema, require_canonical_label_index
from phospy.tables.signalome.common import (
    _require_integer_compatible_index,
    _require_non_negative_integer_index,
)

_ROW_TOTAL_ATOL = 0.05
_VALUE_BOUNDS_ATOL = 1e-6


@dataclass(frozen=True, slots=True)
class SignalomeModulesTable(TableSchema):
    """Schema wrapper for ``signalome_result.signalome_modules.table``."""

    _field_name = "signalome_result.signalome_modules.table"

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=True,
            error_type=self._error_type,
        )
        require_unique_index(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        _require_integer_compatible_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        _require_non_negative_integer_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        require_unique_columns(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_canonical_label_index(
            frame.columns,
            field_name=f"{self._field_name}.columns",
            error_type=self._error_type,
        )
        require_numeric_dataframe(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_finite_numeric_dataframe(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
            allow_missing=False,
        )
        if frame.empty:
            return frame
        values = frame.to_numpy(dtype=float, copy=False)
        if (
            (values < -_VALUE_BOUNDS_ATOL) | (values > 100.0 + _VALUE_BOUNDS_ATOL)
        ).any():
            raise self._error_type(
                f"{self._field_name} values must be between 0.0 and 100.0"
            )
        row_totals = frame.sum(axis=1).to_numpy(dtype=float, copy=False)
        valid_totals = np.isclose(row_totals, 0.0, atol=_ROW_TOTAL_ATOL) | np.isclose(
            row_totals,
            100.0,
            atol=_ROW_TOTAL_ATOL,
        )
        if not valid_totals.all():
            raise self._error_type(
                f"{self._field_name} row totals must be approximately 0.0 or 100.0"
            )
        return frame
