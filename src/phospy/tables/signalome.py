"""Signalome sidecar table schema wrappers."""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.signalomes.constants import (
    CORRELATION_COLUMN,
    CORRELATION_REASON_COLUMN,
    CORRELATION_STATUS_COLUMN,
    DEGREE_COLUMN,
    MODULE_ID_COLUMN,
    MODULE_TOP_KINASE_CANDIDATES_COLUMN,
    MODULE_TOP_KINASE_COLUMN,
    MODULE_TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    MODULE_TOP_KINASE_SELECTION_POLICY_COLUMN,
    MODULE_TOP_KINASE_TIE_COUNT_COLUMN,
    N_SUBSTRATES_COLUMN,
    PROTEIN_COLUMN,
    SITE_CLUSTER_COLUMN,
    SITE_ID_COLUMN,
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
    TOP_KINASE_CANDIDATES_COLUMN,
    TOP_KINASE_COLUMN,
    TOP_KINASE_IS_AMBIGUOUS_COLUMN,
    TOP_KINASE_SELECTION_POLICY_COLUMN,
    TOP_KINASE_TIE_COUNT_COLUMN,
    TOP_KINASE_WEIGHTS_COLUMN,
    TOP_SCORE_COLUMN,
    UNSUPPORTED_KINASE,
    VALID_OBSERVATIONS_COLUMN,
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
from phospy.tables.base import TableSchema, require_canonical_label_index
from phospy.validation.common.dataframes import (
    require_canonical_site_index,
    require_canonical_site_series,
    require_canonical_string_column,
    require_columns,
    require_dataframe,
    require_non_empty_string_column,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
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
_SIGNALOME_ASSIGNMENTS_REQUIRED_COLUMNS = (
    PROTEIN_COLUMN,
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
_KINASE_NETWORK_EDGES_REQUIRED_COLUMNS = (
    SOURCE_KINASE_COLUMN,
    TARGET_KINASE_COLUMN,
    CORRELATION_COLUMN,
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
_ROW_TOTAL_ATOL = 0.05
_VALUE_BOUNDS_ATOL = 1e-6


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
        require_canonical_site_index(
            frame.index,
            field_name=f"{self._field_name}.index",
            error_type=self._error_type,
        )
        if frame.empty:
            return frame
        require_non_empty_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_COLUMN,
            error_type=self._error_type,
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
        if frame.empty:
            return frame
        values = frame.to_numpy(dtype=float, copy=False)
        if not np.isfinite(values).all():
            raise self._error_type(
                f"{self._field_name} must contain finite numeric values"
            )
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


def _require_non_negative_integer_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[WorkflowValidationError],
) -> None:
    _require_integer_compatible_column(
        frame,
        field_name=field_name,
        column_name=column_name,
        error_type=error_type,
        allow_missing=False,
    )
    numeric = pd.to_numeric(frame.loc[:, column_name], errors="coerce")
    values = numeric.to_numpy(dtype=float, copy=False)
    if (values < 0.0).any():
        raise error_type(
            f"{field_name}.{column_name} must contain non-negative integer values"
        )


def _require_integer_compatible_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[WorkflowValidationError],
) -> None:
    numeric = pd.to_numeric(index.to_series(index=index), errors="coerce")
    if numeric.isna().any():
        raise error_type(f"{field_name} must contain integer-compatible labels")
    values = numeric.to_numpy(dtype=float, copy=False)
    if not np.isfinite(values).all():
        raise error_type(f"{field_name} must contain finite integer-compatible labels")
    if not np.isclose(values, np.round(values)).all():
        raise error_type(f"{field_name} must contain integer-compatible labels")


def _require_non_negative_integer_index(
    index: pd.Index,
    *,
    field_name: str,
    error_type: type[WorkflowValidationError],
) -> None:
    numeric = pd.to_numeric(index.to_series(index=index), errors="coerce")
    values = numeric.to_numpy(dtype=float, copy=False)
    if (values < 0.0).any():
        raise error_type(f"{field_name} must contain non-negative integer labels")


def _require_assignment_top_score_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    error_type: type[WorkflowValidationError],
) -> None:
    numeric = pd.to_numeric(frame.loc[:, TOP_SCORE_COLUMN], errors="coerce")
    raw_scores = frame.loc[:, TOP_SCORE_COLUMN]
    missing_mask = raw_scores.isna()
    if (~missing_mask).any():
        finite_values = numeric.loc[~missing_mask].to_numpy(dtype=float, copy=False)
        if not np.isfinite(finite_values).all():
            raise error_type(
                f"{field_name}.{TOP_SCORE_COLUMN} must contain finite numeric values"
            )
    if not missing_mask.any():
        return
    top_kinases = frame.loc[:, TOP_KINASE_COLUMN].astype(str)
    tie_counts = pd.to_numeric(
        frame.loc[:, TOP_KINASE_TIE_COUNT_COLUMN], errors="coerce"
    ).fillna(-1.0)
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
    error_type: type[WorkflowValidationError],
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
    error_type: type[WorkflowValidationError],
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


def _require_numeric_bounds(
    frame: pd.DataFrame,
    *,
    field_name: str,
    column_name: str,
    error_type: type[WorkflowValidationError],
    minimum: float,
    maximum: float,
    allow_missing: bool,
) -> None:
    numeric = pd.to_numeric(frame.loc[:, column_name], errors="coerce")
    if not allow_missing and numeric.isna().any():
        raise error_type(f"{field_name}.{column_name} must not contain missing values")
    values = numeric.dropna().to_numpy(dtype=float, copy=False)
    if ((values < float(minimum)) | (values > float(maximum))).any():
        raise error_type(
            f"{field_name}.{column_name} must be between {minimum:.1f} and {maximum:.1f}"
        )


def _require_correlation_reason_column(
    frame: pd.DataFrame,
    *,
    field_name: str,
    error_type: type[WorkflowValidationError],
) -> None:
    values = frame.loc[:, CORRELATION_REASON_COLUMN]
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
    error_type: type[WorkflowValidationError],
) -> None:
    statuses = frame.loc[:, CORRELATION_STATUS_COLUMN].astype(str)
    numeric = pd.to_numeric(frame.loc[:, CORRELATION_COLUMN], errors="coerce")
    finite_mask = statuses.eq("finite")
    finite_values = numeric.loc[finite_mask]
    if finite_values.isna().any():
        raise error_type(
            f"{field_name}.{CORRELATION_COLUMN} must be present when "
            f"{CORRELATION_STATUS_COLUMN}='finite'"
        )
    finite_array = finite_values.to_numpy(dtype=float, copy=False)
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
