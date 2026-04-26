"""Signalome sidecar table schema wrappers."""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.signalomes.constants import (
    PROTEIN_COLUMN,
    SITE_CLUSTER_COLUMN,
    SITE_ID_COLUMN,
    TOP_KINASE_COLUMN,
)
from phospy.signalomes.context import (
    PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN,
    PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN,
    PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN,
    PROTEIN_SITE_CONTEXT_N_SITES_COLUMN,
    PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN,
    PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN,
    SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN,
    SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN,
    SITE_MEMBERSHIP_INCLUDED_COLUMN,
    SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN,
    SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN,
    SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN,
    SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN,
)
from phospy.tables.base import TableSchema
from phospy.validation.common.dataframes import (
    require_canonical_site_series,
    require_columns,
    require_dataframe,
)

_SITE_MEMBERSHIP_REQUIRED_COLUMNS = (
    SITE_ID_COLUMN,
    PROTEIN_COLUMN,
    SITE_CLUSTER_COLUMN,
    SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN,
    SITE_MEMBERSHIP_INCLUDED_COLUMN,
    SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN,
    SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN,
    TOP_KINASE_COLUMN,
    SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN,
    SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN,
    SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN,
)
_PROTEIN_SITE_CONTEXT_REQUIRED_COLUMNS = (
    PROTEIN_COLUMN,
    PROTEIN_SITE_CONTEXT_N_SITES_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN,
    PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN,
    PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN,
    PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN,
    PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN,
    PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN,
)


@dataclass(frozen=True, slots=True)
class SignalomeSiteContext(TableSchema):
    """Schema wrapper for ``signalome_result.site_membership``."""

    _field_name = "signalome_result.site_membership"
    _error_type = WorkflowValidationError

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
            required_columns=_SITE_MEMBERSHIP_REQUIRED_COLUMNS,
            error_type=self._error_type,
        )
        if frame.empty:
            return frame
        require_canonical_site_series(
            frame.loc[:, SITE_ID_COLUMN],
            field_name=f"{self._field_name}.{SITE_ID_COLUMN}",
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=TOP_KINASE_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN,
            error_type=self._error_type,
        )
        _require_integer_compatible_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_CLUSTER_COLUMN,
            error_type=self._error_type,
            allow_missing=True,
        )
        _require_integer_compatible_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN,
            error_type=self._error_type,
            allow_missing=False,
        )
        _require_integer_compatible_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN,
            error_type=self._error_type,
            allow_missing=False,
        )
        _require_numeric_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN,
            error_type=self._error_type,
            allow_missing=True,
        )
        _require_numeric_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN,
            error_type=self._error_type,
            allow_missing=True,
        )
        _require_boolean_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_MEMBERSHIP_INCLUDED_COLUMN,
            error_type=self._error_type,
        )
        return frame


@dataclass(frozen=True, slots=True)
class SignalomeProteinSiteContext(TableSchema):
    """Schema wrapper for ``signalome_result.protein_site_context``."""

    _field_name = "signalome_result.protein_site_context"
    _error_type = WorkflowValidationError

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
            required_columns=_PROTEIN_SITE_CONTEXT_REQUIRED_COLUMNS,
            error_type=self._error_type,
        )
        if frame.empty:
            return frame
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN,
            error_type=self._error_type,
        )
        _require_integer_compatible_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_N_SITES_COLUMN,
            error_type=self._error_type,
            allow_missing=False,
        )
        _require_integer_compatible_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN,
            error_type=self._error_type,
            allow_missing=False,
        )
        _require_integer_compatible_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN,
            error_type=self._error_type,
            allow_missing=False,
        )
        _require_boolean_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN,
            error_type=self._error_type,
        )
        _require_boolean_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN,
            error_type=self._error_type,
        )
        _require_json_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN,
            error_type=self._error_type,
        )
        _require_json_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN,
            error_type=self._error_type,
        )
        _require_json_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN,
            error_type=self._error_type,
        )
        _require_json_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN,
            error_type=self._error_type,
        )
        return frame


def _require_string_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[WorkflowValidationError],
) -> None:
    values = frame.loc[:, column_name]
    if values.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    if not all(isinstance(value, str) for value in values.tolist()):
        raise error_type(f"{field_name}.{column_name} must contain string values")


def _require_boolean_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[WorkflowValidationError],
) -> None:
    values = frame.loc[:, column_name]
    if values.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    invalid = [
        value for value in values.tolist() if not isinstance(value, (bool, np.bool_))
    ]
    if invalid:
        raise error_type(f"{field_name}.{column_name} must contain boolean values")


def _require_integer_compatible_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[WorkflowValidationError],
    allow_missing: bool,
) -> None:
    numeric = pd.to_numeric(frame.loc[:, column_name], errors="coerce")
    if not allow_missing and numeric.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    finite_values = numeric.dropna().to_numpy(dtype=float, copy=False)
    if not np.isfinite(finite_values).all():
        raise error_type(
            f"{field_name}.{column_name} must contain finite integer-compatible values"
        )
    if not np.isclose(finite_values, np.round(finite_values)).all():
        raise error_type(
            f"{field_name}.{column_name} must contain integer-compatible values"
        )


def _require_numeric_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[WorkflowValidationError],
    allow_missing: bool,
) -> None:
    numeric = pd.to_numeric(frame.loc[:, column_name], errors="coerce")
    if not allow_missing and numeric.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    finite_values = numeric.dropna().to_numpy(dtype=float, copy=False)
    if not np.isfinite(finite_values).all():
        raise error_type(
            f"{field_name}.{column_name} must contain finite numeric values"
        )


def _require_json_string_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[WorkflowValidationError],
) -> None:
    for value in frame.loc[:, column_name].tolist():
        if not isinstance(value, str):
            raise error_type(
                f"{field_name}.{column_name} must contain JSON-encoded strings"
            )
        try:
            json.loads(value)
        except json.JSONDecodeError as exc:
            raise error_type(
                f"{field_name}.{column_name} must contain parseable JSON strings"
            ) from exc
