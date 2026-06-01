from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.validation import WorkflowValidationError
from phospy.science.signalomes.constants import (
    DISPLAY_ID_COLUMN,
    ISOFORM_ID_COLUMN,
    PROTEIN_ACCESSION_COLUMN,
    PROTEIN_COLUMN,
    SITE_CLUSTER_COLUMN,
    SITE_COLUMN,
    SITE_ID_COLUMN,
    SITE_KEY_COLUMN,
    TOP_KINASE_COLUMN,
)
from phospy.science.signalomes.context import (
    PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN,
    PROTEIN_SITE_CONTEXT_DISPLAY_IDS_COLUMN,
    PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN,
    PROTEIN_SITE_CONTEXT_ISOFORM_ID_COLUMN,
    PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN,
    PROTEIN_SITE_CONTEXT_N_SITES_COLUMN,
    PROTEIN_SITE_CONTEXT_PROTEIN_ACCESSION_COLUMN,
    PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_KEY_TO_DISPLAY_ID_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_KEYS_COLUMN,
    PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN,
    SITE_MEMBERSHIP_EXCLUDED_REASON_COLUMN,
    SITE_MEMBERSHIP_GENE_SYMBOL_COLUMN,
    SITE_MEMBERSHIP_INCLUDED_COLUMN,
    SITE_MEMBERSHIP_N_SUPPORTED_KINASES_COLUMN,
    SITE_MEMBERSHIP_PROTEIN_MODULE_COLUMN,
    SITE_MEMBERSHIP_TOP_KINASE_SCORE_COLUMN,
    SITE_MEMBERSHIP_TOP_KINASE_WEIGHT_COLUMN,
)
from phospy.science.sites.validation import require_canonical_site_series
from phospy.tables.base import TableSchema
from phospy.tables.signalome.common import (
    _require_boolean_column,
    _require_integer_compatible_column,
    _require_json_string_column,
    _require_numeric_column,
    _require_string_column,
)
from phospy.validation.common.dataframes import (
    require_columns,
    require_dataframe,
)

_SITE_MEMBERSHIP_REQUIRED_COLUMNS = (
    SITE_KEY_COLUMN,
    DISPLAY_ID_COLUMN,
    SITE_ID_COLUMN,
    SITE_COLUMN,
    PROTEIN_COLUMN,
    PROTEIN_ACCESSION_COLUMN,
    ISOFORM_ID_COLUMN,
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
    PROTEIN_SITE_CONTEXT_SITE_KEYS_COLUMN,
    PROTEIN_SITE_CONTEXT_DISPLAY_IDS_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN,
    PROTEIN_SITE_CONTEXT_N_DISTINCT_SITE_CLUSTERS_COLUMN,
    PROTEIN_SITE_CONTEXT_PROTEIN_MODULE_COLUMN,
    PROTEIN_SITE_CONTEXT_MULTI_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_AMBIGUOUS_CONTEXT_COLUMN,
    PROTEIN_SITE_CONTEXT_GENE_SYMBOL_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_PROTEIN_ACCESSION_COLUMN,
    PROTEIN_SITE_CONTEXT_ISOFORM_ID_COLUMN,
    PROTEIN_SITE_CONTEXT_TOP_KINASES_BY_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_MODULE_IDS_BY_SITE_COLUMN,
    PROTEIN_SITE_CONTEXT_SITE_KEY_TO_DISPLAY_ID_COLUMN,
)


@dataclass(frozen=True, slots=True)
class SignalomeSiteContext(TableSchema):
    """Schema wrapper for ``signalome_result.site_membership``."""

    _field_name = "signalome_result.site_membership"
    _error_type = WorkflowValidationError

    def _validate_frame(self, frame: pd.DataFrame) -> pd.DataFrame:
        frame = _normalize_legacy_site_membership_frame(frame)
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
            frame.loc[:, DISPLAY_ID_COLUMN],
            field_name=f"{self._field_name}.{DISPLAY_ID_COLUMN}",
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_KEY_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=DISPLAY_ID_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_ID_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=SITE_COLUMN,
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
            column_name=PROTEIN_ACCESSION_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=ISOFORM_ID_COLUMN,
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
        frame = _normalize_legacy_protein_site_context_frame(frame)
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
            column_name=PROTEIN_SITE_CONTEXT_SITE_KEYS_COLUMN,
            error_type=self._error_type,
        )
        _require_json_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_DISPLAY_IDS_COLUMN,
            error_type=self._error_type,
        )
        _require_json_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_SITE_CLUSTERS_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_SITE_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_PROTEIN_ACCESSION_COLUMN,
            error_type=self._error_type,
        )
        _require_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_ISOFORM_ID_COLUMN,
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
        _require_json_string_column(
            frame,
            field_name=self._field_name,
            column_name=PROTEIN_SITE_CONTEXT_SITE_KEY_TO_DISPLAY_ID_COLUMN,
            error_type=self._error_type,
        )
        return frame


def _normalize_legacy_site_membership_frame(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    if (
        SITE_KEY_COLUMN not in normalized.columns
        and SITE_ID_COLUMN in normalized.columns
    ):
        normalized.loc[:, SITE_KEY_COLUMN] = (
            normalized.loc[:, SITE_ID_COLUMN].fillna("").astype(str)
        )
    if (
        DISPLAY_ID_COLUMN not in normalized.columns
        and SITE_ID_COLUMN in normalized.columns
    ):
        normalized.loc[:, DISPLAY_ID_COLUMN] = (
            normalized.loc[:, SITE_ID_COLUMN].fillna("").astype(str)
        )
    if SITE_COLUMN not in normalized.columns:
        normalized.loc[:, SITE_COLUMN] = ""
    if PROTEIN_ACCESSION_COLUMN not in normalized.columns:
        normalized.loc[:, PROTEIN_ACCESSION_COLUMN] = ""
    if ISOFORM_ID_COLUMN not in normalized.columns:
        normalized.loc[:, ISOFORM_ID_COLUMN] = ""
    return normalized


def _normalize_legacy_protein_site_context_frame(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    if (
        PROTEIN_SITE_CONTEXT_SITE_KEYS_COLUMN not in normalized.columns
        and PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN in normalized.columns
    ):
        normalized.loc[:, PROTEIN_SITE_CONTEXT_SITE_KEYS_COLUMN] = (
            normalized.loc[:, PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN]
            .fillna("[]")
            .astype(str)
        )
    if (
        PROTEIN_SITE_CONTEXT_DISPLAY_IDS_COLUMN not in normalized.columns
        and PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN in normalized.columns
    ):
        normalized.loc[:, PROTEIN_SITE_CONTEXT_DISPLAY_IDS_COLUMN] = (
            normalized.loc[:, PROTEIN_SITE_CONTEXT_SITE_IDS_COLUMN]
            .fillna("[]")
            .astype(str)
        )
    if PROTEIN_SITE_CONTEXT_SITE_COLUMN not in normalized.columns:
        normalized.loc[:, PROTEIN_SITE_CONTEXT_SITE_COLUMN] = ""
    if PROTEIN_SITE_CONTEXT_PROTEIN_ACCESSION_COLUMN not in normalized.columns:
        normalized.loc[:, PROTEIN_SITE_CONTEXT_PROTEIN_ACCESSION_COLUMN] = ""
    if PROTEIN_SITE_CONTEXT_ISOFORM_ID_COLUMN not in normalized.columns:
        normalized.loc[:, PROTEIN_SITE_CONTEXT_ISOFORM_ID_COLUMN] = ""
    if PROTEIN_SITE_CONTEXT_SITE_KEY_TO_DISPLAY_ID_COLUMN not in normalized.columns:
        normalized.loc[:, PROTEIN_SITE_CONTEXT_SITE_KEY_TO_DISPLAY_ID_COLUMN] = "{}"
    return normalized
