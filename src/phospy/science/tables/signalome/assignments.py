from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.frames.validation import (
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_string_index,
    require_unique_index,
)
from phospy.science.signalomes.constants import (
    DISPLAY_ID_COLUMN,
    GENE_SYMBOL_COLUMN,
    ISOFORM_ID_COLUMN,
    MODULE_ID_COLUMN,
    MODULE_TOP_KINASE_CANDIDATES_COLUMN,
    MODULE_TOP_KINASE_COLUMN,
    MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
    MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
    PROTEIN_ACCESSION_COLUMN,
    PROTEIN_COLUMN,
    SITE_COLUMN,
    SITE_KEY_COLUMN,
    TOP_KINASE_CANDIDATES_COLUMN,
    TOP_KINASE_COLUMN,
    TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    TOP_KINASE_SELECTION_POLICY_COLUMN,
    TOP_KINASE_TIE_COUNT_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_SCORE_COLUMN,
    UNSUPPORTED_KINASE,
)
from phospy.science.sites.identity_contracts import (
    enforce_analysis_ready_site_key_index,
    enforce_site_key_column_matches_index,
)
from phospy.science.tables.base import TableSchema, ValidationErrorType
from phospy.science.tables.signalome.common import (
    _column_series,
    _numeric_series,
    _require_boolean_column,
    _require_integer_compatible_column,
    _require_non_negative_integer_column,
    _require_string_column,
)

_SIGNALOME_ASSIGNMENTS_REQUIRED_COLUMNS = (
    SITE_KEY_COLUMN,
    DISPLAY_ID_COLUMN,
    GENE_SYMBOL_COLUMN,
    SITE_COLUMN,
    PROTEIN_COLUMN,
    PROTEIN_ACCESSION_COLUMN,
    ISOFORM_ID_COLUMN,
    MODULE_ID_COLUMN,
    TOP_KINASE_COLUMN,
    TOP_SCORE_COLUMN,
    TOP_KINASE_CANDIDATES_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_KINASE_TIE_COUNT_COLUMN,
    TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    TOP_KINASE_SELECTION_POLICY_COLUMN,
    MODULE_TOP_KINASE_COLUMN,
    MODULE_TOP_KINASE_CANDIDATES_COLUMN,
    MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
    MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
)


@dataclass(frozen=True, slots=True)
class SignalomeAssignmentsTable(TableSchema):
    """Schema wrapper for ``signalome_result.module_assignments.table``."""

    _field_name = "signalome_result.module_assignments.table"

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
            required_columns=_SIGNALOME_ASSIGNMENTS_REQUIRED_COLUMNS,
            error_type=self._error_type,
        )
        require_unique_index(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        require_string_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        enforce_analysis_ready_site_key_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        if frame.empty:
            return frame
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_KEY_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=DISPLAY_ID_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=GENE_SYMBOL_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_ACCESSION_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=ISOFORM_ID_COLUMN,
            error_type=self._error_type,
        )
        enforce_site_key_column_matches_index(
            site_metadata=frame,
            field_name=self._field_name,
            error_type=self._error_type,
            site_key_column=SITE_KEY_COLUMN,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=TOP_KINASE_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=TOP_KINASE_SELECTION_POLICY_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=MODULE_TOP_KINASE_COLUMN,
            error_type=self._error_type,
        )
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
            error_type=self._error_type,
        )
        _require_integer_compatible_column(
            frame,
            field_name=self._field_name,
            column_name=MODULE_ID_COLUMN,
            error_type=self._error_type,
            allow_missing=False,
        )
        _require_non_negative_integer_column(
            frame,
            field_name=self._field_name,
            column_name=MODULE_ID_COLUMN,
            error_type=self._error_type,
        )
        _require_non_negative_integer_column(
            frame,
            field_name=self._field_name,
            column_name=TOP_KINASE_TIE_COUNT_COLUMN,
            error_type=self._error_type,
        )
        _require_non_negative_integer_column(
            frame,
            field_name=self._field_name,
            column_name=MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
            error_type=self._error_type,
        )
        _require_boolean_column(
            frame,
            field_name=self._field_name,
            column_name=TOP_KINASE_IS_AMBIGUOUS_COLUMN,
            error_type=self._error_type,
        )
        _require_boolean_column(
            frame,
            field_name=self._field_name,
            column_name=MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
            error_type=self._error_type,
        )
        _require_assignment_top_score_column(
            frame,
            field_name=self._field_name,
            error_type=self._error_type,
        )
        _require_kinase_candidates_sequence_column(
            frame,
            field_name=self._field_name,
            column_name=TOP_KINASE_CANDIDATES_COLUMN,
            error_type=self._error_type,
        )
        _require_kinase_candidates_sequence_column(
            frame,
            field_name=self._field_name,
            column_name=MODULE_TOP_KINASE_CANDIDATES_COLUMN,
            error_type=self._error_type,
        )
        _require_kinase_weights_sequence_column(
            frame,
            field_name=self._field_name,
            column_name=TOP_KINASE_WEIGHTS_COLUMN,
            error_type=self._error_type,
        )
        return frame


def _require_assignment_top_score_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> None:
    numeric = _numeric_series(frame, TOP_SCORE_COLUMN)
    raw_scores = _column_series(frame, TOP_SCORE_COLUMN)
    missing_mask = raw_scores.isna()
    if (~missing_mask).any():
        finite_values = numeric.loc[~missing_mask].to_numpy(
            dtype="float64",
            copy=False,
        )
        if not np.isfinite(finite_values).all():
            raise error_type(
                f"{field_name}.{TOP_SCORE_COLUMN} must contain finite numeric values"
            )
    if not missing_mask.any():
        return
    top_kinases = _column_series(frame, TOP_KINASE_COLUMN).astype(str)
    tie_counts = _numeric_series(frame, TOP_KINASE_TIE_COUNT_COLUMN).fillna(-1.0)
    allowed_missing_mask = top_kinases.eq(UNSUPPORTED_KINASE) & tie_counts.eq(0.0)
    invalid_missing = missing_mask & ~allowed_missing_mask
    if bool(invalid_missing.any()):
        raise error_type(
            f"{field_name}.{TOP_SCORE_COLUMN} must be finite for supported top kinases; "
            f"missing values are only allowed when {TOP_KINASE_COLUMN}='{UNSUPPORTED_KINASE}' "
            f"and {TOP_KINASE_TIE_COUNT_COLUMN}=0"
        )


def _require_kinase_candidates_sequence_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> None:
    for row_label, value in frame.loc[:, column_name].items():
        if not isinstance(value, (tuple, list)):
            raise error_type(
                f"{field_name}.{column_name} must contain tuple/list values; "
                f"invalid value at index '{row_label}'"
            )
        for kinase in value:
            if not isinstance(kinase, str) or kinase.strip() == "":
                raise error_type(
                    f"{field_name}.{column_name} must contain non-empty string kinase IDs"
                )


def _require_kinase_weights_sequence_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: ValidationErrorType,
) -> None:
    for row_label, value in frame.loc[:, column_name].items():
        if not isinstance(value, (tuple, list)):
            raise error_type(
                f"{field_name}.{column_name} must contain tuple/list values; "
                f"invalid value at index '{row_label}'"
            )
        for pair in value:
            if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                raise error_type(
                    f"{field_name}.{column_name} entries must be (kinase, weight) pairs"
                )
            kinase_id, weight = pair
            if not isinstance(kinase_id, str) or kinase_id.strip() == "":
                raise error_type(
                    f"{field_name}.{column_name} must contain non-empty string kinase IDs"
                )
            try:
                numeric_weight = float(weight)
            except (TypeError, ValueError) as exc:
                raise error_type(
                    f"{field_name}.{column_name} weights must be numeric"
                ) from exc
            if not np.isfinite(numeric_weight):
                raise error_type(
                    f"{field_name}.{column_name} weights must be finite numeric values"
                )
            if numeric_weight < 0.0:
                raise error_type(
                    f"{field_name}.{column_name} weights must be non-negative"
                )
