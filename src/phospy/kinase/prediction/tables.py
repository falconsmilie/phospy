"""Passive table contracts for Kinase prediction."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.frames.ownership import own_dataframe
from phospy.frames.table_schema import TableSchema
from phospy.kinase._matrix_validation import _validate_kinase_score_like_matrix


@dataclass(frozen=True, slots=True, eq=False)
class KinasePredictionMatrix(TableSchema):
    """Schema wrapper for ``prediction_result.pred_mat``."""

    field_name: str = field(
        default="prediction_result.pred_mat",
        repr=False,
        compare=False,
    )
    allow_missing: bool = field(default=True, repr=False, compare=False)
    enforce_unit_interval: bool = field(default=True, repr=False, compare=False)

    _field_name = "prediction_result.pred_mat"
    _error_type = PhosPyValidationError

    def _own_frame(self, assume_owned: bool) -> pd.DataFrame:
        return own_dataframe(
            self.frame,
            field_name=self.field_name,
            error_type=self._error_type,
            assume_owned=assume_owned,
        )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        _validate_kinase_score_like_matrix(
            frame=frame,
            field_name=self.field_name,
            error_type=self._error_type,
            allow_missing=self.allow_missing,
            enforce_unit_interval=self.enforce_unit_interval,
            require_site_index=True,
        )
        return frame


__all__ = ["KinasePredictionMatrix"]
