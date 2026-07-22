"""Activity-stage scientific table/series wrappers."""

from __future__ import annotations

from dataclasses import InitVar, dataclass, field
from typing import ClassVar

import numpy as np
import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.frames.ownership import own_dataframe, own_series
from phospy.frames.table_schema import (
    TableSchema,
    ValidationErrorType,
    require_canonical_label_index,
)
from phospy.frames.validation import (
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_finite_numeric_dataframe,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.science.sites.validation import (
    require_canonical_site_series,
    require_site_key_series,
)


@dataclass(frozen=True, slots=True)
class ActivityMatrix(TableSchema):
    """Schema wrapper for activity score matrices."""

    field_name: str = field(
        default="activity_result.weighted_activity",
        repr=False,
        compare=False,
    )
    allow_missing: bool = field(default=True, repr=False, compare=False)

    _field_name = "activity_result.weighted_activity"
    _error_type = PhosPyValidationError

    def _own_frame(self, assume_owned: bool) -> pd.DataFrame:
        return own_dataframe(
            self.frame,
            field_name=self.field_name,
            error_type=self._error_type,
            assume_owned=assume_owned,
        )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self.field_name,
            allow_empty=True,
            error_type=self._error_type,
        )
        require_numeric_dataframe(
            frame,
            field_name=self.field_name,
            error_type=self._error_type,
        )
        require_finite_numeric_dataframe(
            frame,
            field_name=self.field_name,
            error_type=self._error_type,
            allow_missing=self.allow_missing,
        )
        require_unique_index(
            frame,
            field_name=self.field_name,
            error_type=self._error_type,
        )
        require_unique_columns(
            frame,
            field_name=self.field_name,
            error_type=self._error_type,
        )
        require_canonical_label_index(
            frame.index,
            field_name=f"{self.field_name}.index",
            error_type=self._error_type,
        )
        require_canonical_label_index(
            frame.columns,
            field_name=f"{self.field_name}.columns",
            error_type=self._error_type,
        )
        return frame


@dataclass(frozen=True, slots=True)
class ActivityCountMatrix(TableSchema):
    """Schema wrapper for condition-specific activity count matrices."""

    field_name: str = field(
        default="activity_result.activity_substrate_counts",
        repr=False,
        compare=False,
    )

    _field_name = "activity_result.activity_substrate_counts"
    _error_type = PhosPyValidationError

    def _own_frame(self, assume_owned: bool) -> pd.DataFrame:
        return own_dataframe(
            self.frame,
            field_name=self.field_name,
            error_type=self._error_type,
            assume_owned=assume_owned,
        )

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self.field_name,
            allow_empty=True,
            error_type=self._error_type,
        )
        require_numeric_dataframe(
            frame,
            field_name=self.field_name,
            error_type=self._error_type,
        )
        require_unique_index(
            frame,
            field_name=self.field_name,
            error_type=self._error_type,
        )
        require_unique_columns(
            frame,
            field_name=self.field_name,
            error_type=self._error_type,
        )
        require_canonical_label_index(
            frame.index,
            field_name=f"{self.field_name}.index",
            error_type=self._error_type,
        )
        require_canonical_label_index(
            frame.columns,
            field_name=f"{self.field_name}.columns",
            error_type=self._error_type,
        )
        if frame.empty:
            return frame.astype("int64")

        values = frame.to_numpy(dtype="float64", copy=False)
        if not np.isfinite(values).all():
            raise self._error_type(f"{self.field_name} must contain finite counts")
        if (values < 0.0).any():
            raise self._error_type(
                f"{self.field_name} must contain non-negative counts"
            )
        if not np.isclose(values, np.round(values)).all():
            raise self._error_type(
                f"{self.field_name} must contain integer-compatible counts"
            )
        return pd.DataFrame(
            np.round(values).astype("int64"),
            index=frame.index.copy(),
            columns=frame.columns.copy(),
            dtype="int64",
        )


@dataclass(frozen=True, slots=True)
class ActivityTargetTable(TableSchema):
    """Schema wrapper for ``activity_result.target_table``."""

    _field_name = "activity_result.target_table"
    _error_type = PhosPyValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=True,
            error_type=self._error_type,
        )
        require_columns(
            frame,
            field_name=self._field_name,
            required_columns=("site_id", "kinase", "score"),
            error_type=self._error_type,
        )
        if frame.empty:
            return frame
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="site_id",
            error_type=self._error_type,
        )
        site_series = frame.loc[:, "site_id"]
        try:
            require_canonical_site_series(
                site_series,
                field_name=f"{self._field_name}.site_id",
                error_type=self._error_type,
            )
        except self._error_type:
            require_site_key_series(
                site_series,
                field_name=f"{self._field_name}.site_id",
                error_type=self._error_type,
            )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name="kinase",
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name="kinase",
            error_type=self._error_type,
        )
        _require_numeric_column(
            frame,
            field_name=self._field_name,
            column_name="score",
            error_type=self._error_type,
        )
        score_values = pd.to_numeric(
            _column_series(frame, "score"),
            errors="coerce",
        ).to_numpy(
            dtype="float64",
            copy=False,
        )
        invalid_mask = (score_values < 0.0) | (score_values > 1.0)
        if invalid_mask.any():
            raise self._error_type(
                f"{self._field_name}.score must be between 0.0 and 1.0"
            )
        return frame


@dataclass(frozen=True, slots=True)
class ActivityStatisticsTable(TableSchema):
    """Schema wrapper for method-specific activity statistics tables."""

    _field_name = "activity_result.statistics_table"
    _error_type = PhosPyValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        require_dataframe(
            frame,
            field_name=self._field_name,
            allow_empty=True,
            error_type=self._error_type,
        )
        require_unique_columns(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_columns(
            frame,
            field_name=self._field_name,
            required_columns=(
                "kinase",
                "condition",
                "z_score",
                "p_value",
                "q_value",
                "n_substrates",
                "n_background_sites",
                "evidence_threshold",
                "evidence_threshold_operator",
                "evidence_threshold_description",
                "min_substrates",
                "computability_status",
                "reason",
            ),
            error_type=self._error_type,
        )
        if frame.empty:
            return frame

        for column_name in (
            "kinase",
            "condition",
            "evidence_threshold_operator",
            "evidence_threshold_description",
            "computability_status",
            "reason",
        ):
            values = _column_series(frame, column_name)
            if values.isna().any():
                raise self._error_type(
                    f"{self._field_name}.{column_name} must not contain missing values"
                )
            if not values.astype(str).str.len().ge(0).all():
                raise self._error_type(
                    f"{self._field_name}.{column_name} must contain string values"
                )

        _require_numeric_column_allowing_missing(
            frame,
            field_name=self._field_name,
            column_name="z_score",
            error_type=self._error_type,
        )
        _require_numeric_column_allowing_missing(
            frame,
            field_name=self._field_name,
            column_name="p_value",
            error_type=self._error_type,
        )
        _require_numeric_column_allowing_missing(
            frame,
            field_name=self._field_name,
            column_name="q_value",
            error_type=self._error_type,
        )
        for count_column in ("n_substrates", "n_background_sites", "min_substrates"):
            numeric = _require_integer_compatible_column(
                frame,
                field_name=self._field_name,
                column_name=count_column,
                error_type=self._error_type,
            )
            if (numeric < 0.0).any():
                raise self._error_type(
                    f"{self._field_name}.{count_column} must be non-negative"
                )
        evidence_threshold = _require_numeric_column_allowing_missing(
            frame,
            field_name=self._field_name,
            column_name="evidence_threshold",
            error_type=self._error_type,
        )
        if (evidence_threshold < 0.0).any() or (evidence_threshold > 1.0).any():
            raise self._error_type(
                f"{self._field_name}.evidence_threshold must be between 0.0 and 1.0"
            )
        p_values = _require_numeric_column_allowing_missing(
            frame,
            field_name=self._field_name,
            column_name="p_value",
            error_type=self._error_type,
        )
        if (p_values < 0.0).any() or (p_values > 1.0).any():
            raise self._error_type(
                f"{self._field_name}.p_value must be between 0.0 and 1.0 when present"
            )
        q_values = _require_numeric_column_allowing_missing(
            frame,
            field_name=self._field_name,
            column_name="q_value",
            error_type=self._error_type,
        )
        if (q_values < 0.0).any() or (q_values > 1.0).any():
            raise self._error_type(
                f"{self._field_name}.q_value must be between 0.0 and 1.0 when present"
            )
        return frame


@dataclass(frozen=True, slots=True)
class SeriesSchema:
    """Base wrapper for one owned, validated Series contract."""

    series: pd.Series
    _assume_owned: InitVar[bool] = False

    _field_name: ClassVar[str] = "series"
    _error_type: ClassVar[ValidationErrorType] = PhosPyValidationError

    def __post_init__(self, _assume_owned: bool) -> None:
        series = own_series(
            self.series,
            field_name=self._field_name,
            error_type=self._error_type,
            assume_owned=_assume_owned,
        )
        validated = self._validate_series(series)
        object.__setattr__(self, "series", validated)

    def _validate_series(self, series: pd.Series) -> pd.Series:
        return series

    @classmethod
    def _from_owned(cls, *, series: pd.Series) -> SeriesSchema:
        return cls(series=series, _assume_owned=True)


@dataclass(frozen=True, slots=True)
class ActivityCountSeries(SeriesSchema):
    """Schema wrapper for activity count series."""

    field_name: str = field(
        default="activity_result.thresholded_substrate_counts",
        repr=False,
        compare=False,
    )
    allow_empty: bool = field(default=True, repr=False, compare=False)

    _field_name = "activity_result.thresholded_substrate_counts"
    _error_type = PhosPyValidationError

    def __post_init__(self, _assume_owned: bool) -> None:
        series = own_series(
            self.series,
            field_name=self.field_name,
            error_type=self._error_type,
            assume_owned=_assume_owned,
        )
        validated = self._validate_series(series)
        object.__setattr__(self, "series", validated)

    def _validate_series(self, series: pd.Series) -> pd.Series:
        if not self.allow_empty and series.empty:
            raise self._error_type(f"{self.field_name} must be non-empty")
        if not series.index.is_unique:
            raise self._error_type(f"{self.field_name}.index must be unique")
        require_canonical_label_index(
            series.index,
            field_name=f"{self.field_name}.index",
            error_type=self._error_type,
        )
        coerced = pd.to_numeric(series, errors="coerce")
        if coerced.isna().any():
            raise self._error_type(
                f"{self.field_name} must contain integer-compatible counts"
            )
        values = coerced.to_numpy(dtype="float64", copy=False)
        if not np.isfinite(values).all():
            raise self._error_type(f"{self.field_name} must contain finite counts")
        if (values < 0.0).any():
            raise self._error_type(
                f"{self.field_name} must contain non-negative counts"
            )
        if not np.isclose(values, np.round(values)).all():
            raise self._error_type(
                f"{self.field_name} must contain integer-compatible counts"
            )
        return pd.Series(
            np.round(values).astype("int64"),
            index=series.index.copy(),
            name=series.name,
            dtype="int64",
        )


def _require_numeric_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[PhosPyValidationError],
) -> None:
    values = pd.to_numeric(_column_series(frame, column_name), errors="coerce")
    if values.isna().any():
        raise error_type(f"{field_name}.{column_name} must contain numeric values")
    if not np.isfinite(values.to_numpy(dtype="float64", copy=False)).all():
        raise error_type(
            f"{field_name}.{column_name} must contain finite numeric values"
        )


def _require_numeric_column_allowing_missing(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[PhosPyValidationError],
) -> np.ndarray:
    raw_values = _column_series(frame, column_name)
    values = pd.to_numeric(raw_values, errors="coerce")
    coerced_missing = values.isna() & raw_values.notna()
    if coerced_missing.any():
        raise error_type(
            f"{field_name}.{column_name} must contain numeric values when present"
        )
    array = values.to_numpy(dtype="float64", copy=False)
    finite_mask = np.isfinite(array)
    nan_mask = np.isnan(array)
    if ~(finite_mask | nan_mask).all():
        raise error_type(
            f"{field_name}.{column_name} must contain finite numeric values when present"
        )
    return array[finite_mask]


def _require_integer_compatible_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[PhosPyValidationError],
) -> np.ndarray:
    values = pd.to_numeric(_column_series(frame, column_name), errors="coerce")
    if values.isna().any():
        raise error_type(f"{field_name}.{column_name} must be integer-compatible")
    array = values.to_numpy(dtype="float64", copy=False)
    if not np.isfinite(array).all():
        raise error_type(f"{field_name}.{column_name} must be finite")
    if not np.isclose(array, np.round(array)).all():
        raise error_type(f"{field_name}.{column_name} must be integer-compatible")
    return array


def _column_series(frame: pd.DataFrame, column_name: str) -> pd.Series:
    values = frame.loc[:, column_name]
    if not isinstance(values, pd.Series):
        raise PhosPyValidationError(
            f"Expected '{column_name}' column lookup to return a Series"
        )
    return values
