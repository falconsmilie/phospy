"""Kinase/prediction scientific table wrappers."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from phospy._frame_ownership import own_dataframe
from phospy.errors.validation import PhosPyValidationError
from phospy.tables.base import TableSchema, require_canonical_label_index
from phospy.validation.common.dataframes import (
    require_canonical_site_index,
    require_dataframe,
    require_finite_numeric_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.common.numeric_frames import require_numeric_unit_interval


@dataclass(frozen=True, slots=True)
class KinaseScoreMatrix(TableSchema):
    """Schema wrapper for kinase scoring matrices."""

    field_name: str = field(
        default="scoring_result.profile_scores",
        repr=False,
        compare=False,
    )
    allow_missing: bool = field(default=True, repr=False, compare=False)
    enforce_unit_interval: bool = field(default=False, repr=False, compare=False)

    _field_name = "scoring_result.profile_scores"
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
        )
        return frame


@dataclass(frozen=True, slots=True)
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
        )
        return frame


def _validate_kinase_score_like_matrix(
    *,
    frame: pd.DataFrame,
    field_name: str,
    error_type: type[PhosPyValidationError],
    allow_missing: bool,
    enforce_unit_interval: bool,
) -> None:
    require_dataframe(
        frame,
        field_name=field_name,
        allow_empty=False,
        error_type=error_type,
    )
    require_unique_columns(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    require_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    require_finite_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=error_type,
        allow_missing=allow_missing,
    )
    require_unique_index(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    require_canonical_site_index(
        frame.index,
        field_name=f"{field_name}.index",
        error_type=error_type,
    )
    require_canonical_label_index(
        frame.columns,
        field_name=f"{field_name}.columns",
        error_type=error_type,
    )
    if enforce_unit_interval:
        require_numeric_unit_interval(
            frame,
            field_name=field_name,
            error_type=error_type,
        )
