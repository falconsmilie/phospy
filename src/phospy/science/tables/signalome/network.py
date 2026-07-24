from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.frames.table_schema import (
    TableSchema,
    ValidationErrorType,
    require_canonical_label_index,
)
from phospy.frames.validation import (
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_unique_index,
)
from phospy.science.signalomes.constants import (
    CORRELATION_COLUMN,
    CORRELATION_REASON_COLUMN,
    CORRELATION_STATUS_COLUMN,
    DEGREE_COLUMN,
    N_SUBSTRATES_COLUMN,
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
    VALID_OBSERVATIONS_COLUMN,
)
from phospy.science.tables.signalome.common import (
    _column_series,
    _numeric_series,
    _require_non_negative_integer_column,
    _require_numeric_bounds,
    _require_numeric_column,
)

_KINASE_NETWORK_EDGES_REQUIRED_COLUMNS = (
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
    CORRELATION_COLUMN,
    VALID_OBSERVATIONS_COLUMN,
)
_KINASE_NETWORK_NODES_REQUIRED_COLUMNS = (
    DEGREE_COLUMN,
    N_SUBSTRATES_COLUMN,
)
_KINASE_NETWORK_CANDIDATE_CORRELATIONS_REQUIRED_COLUMNS = (
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
    CORRELATION_COLUMN,
    CORRELATION_STATUS_COLUMN,
    VALID_OBSERVATIONS_COLUMN,
    CORRELATION_REASON_COLUMN,
)
_ALLOWED_CORRELATION_STATUSES = frozenset(
    {
        "finite",
        "constant_profile",
        "insufficient_observations",
        "missing_values",
        "non_finite_values",
        "undefined",
    }
)


@dataclass(frozen=True, slots=True)
class KinaseNetworkEdgesTable(TableSchema):
    """Schema wrapper for ``signalome_result.kinase_network.edges``."""

    _field_name = "signalome_result.kinase_network.edges"

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
            required_columns=_KINASE_NETWORK_EDGES_REQUIRED_COLUMNS,
            error_type=self._error_type,
        )
        if frame.empty:
            return frame
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=SOURCE_KINASE_COLUMN,
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name=SOURCE_KINASE_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=TARGET_KINASE_COLUMN,
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name=TARGET_KINASE_COLUMN,
            error_type=self._error_type,
        )
        _require_numeric_column(
            frame,
            field_name=self._field_name,
            column_name=CORRELATION_COLUMN,
            error_type=self._error_type,
            allow_missing=False,
        )
        _require_numeric_bounds(
            frame,
            field_name=self._field_name,
            column_name=CORRELATION_COLUMN,
            error_type=self._error_type,
            minimum=-1.0,
            maximum=1.0,
            allow_missing=False,
        )
        _require_non_negative_integer_column(
            frame,
            field_name=self._field_name,
            column_name=VALID_OBSERVATIONS_COLUMN,
            error_type=self._error_type,
        )
        return frame


@dataclass(frozen=True, slots=True)
class KinaseNetworkNodesTable(TableSchema):
    """Schema wrapper for ``signalome_result.kinase_network.nodes``."""

    _field_name = "signalome_result.kinase_network.nodes"

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
            required_columns=_KINASE_NETWORK_NODES_REQUIRED_COLUMNS,
            error_type=self._error_type,
        )
        require_unique_index(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_canonical_label_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        if frame.empty:
            return frame
        _require_non_negative_integer_column(
            frame,
            field_name=self._field_name,
            column_name=DEGREE_COLUMN,
            error_type=self._error_type,
        )
        _require_non_negative_integer_column(
            frame,
            field_name=self._field_name,
            column_name=N_SUBSTRATES_COLUMN,
            error_type=self._error_type,
        )
        return frame


@dataclass(frozen=True, slots=True)
class KinaseNetworkCandidateCorrelationsTable(TableSchema):
    """Schema wrapper for ``signalome_result.kinase_network.candidate_correlations``."""

    _field_name = "signalome_result.kinase_network.candidate_correlations"

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
            required_columns=_KINASE_NETWORK_CANDIDATE_CORRELATIONS_REQUIRED_COLUMNS,
            error_type=self._error_type,
        )
        if frame.empty:
            return frame
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=SOURCE_KINASE_COLUMN,
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name=SOURCE_KINASE_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=TARGET_KINASE_COLUMN,
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name=TARGET_KINASE_COLUMN,
            error_type=self._error_type,
        )
        _require_numeric_column(
            frame,
            field_name=self._field_name,
            column_name=CORRELATION_COLUMN,
            error_type=self._error_type,
            allow_missing=True,
        )
        _require_non_negative_integer_column(
            frame,
            field_name=self._field_name,
            column_name=VALID_OBSERVATIONS_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=CORRELATION_STATUS_COLUMN,
            error_type=self._error_type,
        )
        require_canonical_string_column(
            frame,
            field_name=self._field_name,
            column_name=CORRELATION_STATUS_COLUMN,
            error_type=self._error_type,
        )
        statuses = frame.loc[:, CORRELATION_STATUS_COLUMN].astype(str)
        unknown_statuses = sorted(
            status
            for status in set(statuses.tolist())
            if status not in _ALLOWED_CORRELATION_STATUSES
        )
        if unknown_statuses:
            preview = ", ".join(unknown_statuses[:3])
            suffix = "..." if len(unknown_statuses) > 3 else ""
            raise self._error_type(
                f"{self._field_name}.{CORRELATION_STATUS_COLUMN} contains unsupported "
                f"values: {preview}{suffix}"
            )
        _require_correlation_reason_column(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        _require_candidate_correlation_value_semantics(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        return frame


def _require_correlation_reason_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> None:
    values = _column_series(frame, CORRELATION_REASON_COLUMN)
    for value in values.tolist():
        if pd.isna(value):
            continue
        if not isinstance(value, str) or value.strip() == "":
            raise error_type(
                f"{field_name}.{CORRELATION_REASON_COLUMN} must contain strings or missing values"
            )


def _require_candidate_correlation_value_semantics(
    frame: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> None:
    statuses = _column_series(frame, CORRELATION_STATUS_COLUMN).astype(str)
    numeric = _numeric_series(frame, CORRELATION_COLUMN)
    finite_mask = statuses.eq("finite")
    finite_values = numeric.loc[finite_mask]
    if finite_values.isna().any():
        raise error_type(
            f"{field_name}.{CORRELATION_COLUMN} must be present when "
            f"{CORRELATION_STATUS_COLUMN}='finite'"
        )
    finite_array = finite_values.to_numpy(dtype="float64", copy=False)
    if not np.isfinite(finite_array).all():
        raise error_type(
            f"{field_name}.{CORRELATION_COLUMN} must be finite when "
            f"{CORRELATION_STATUS_COLUMN}='finite'"
        )
    if ((finite_array < -1.0) | (finite_array > 1.0)).any():
        raise error_type(
            f"{field_name}.{CORRELATION_COLUMN} must be between -1.0 and 1.0"
        )
    undefined_mask = ~finite_mask
    if numeric.loc[undefined_mask].notna().any():
        raise error_type(
            f"{field_name}.{CORRELATION_COLUMN} must be missing when "
            f"{CORRELATION_STATUS_COLUMN} is not 'finite'"
        )
