from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import pandas as pd

from ..constants import (
    DEFAULT_PHOSPHO_COLS,
    DEFAULT_TOTAL_COLS,
    GENE_P_SITE_COLUMN,
    LOCALIZATION_PROB_COLUMN,
    PHOSPHO_GENE_COLUMN,
    PHOSPHO_REQUIRED_METADATA_COLUMNS,
    PHOSPHO_UID_COLUMN,
    TOTAL_GENE_COLUMN,
)
from .errors import TableSchemaError
from .frames import (
    coerce_numeric_columns,
    require_columns,
    require_dataframe,
    require_finite_numeric_values,
    require_non_null_column_names,
    require_non_null_index,
    require_non_null_values,
    require_unique_columns,
    require_unique_index,
    require_value_range,
)
from .identifiers import require_splitable_gene_p_site

if TYPE_CHECKING:
    from ..prediction.models import PredMatResult


class TotalInputSchema:
    """Validate raw total-proteome input tables."""

    @classmethod
    def validate(
        cls,
        frame: pd.DataFrame,
        *,
        total_cols: Sequence[str] | None = None,
        context: str = "total input table",
    ) -> pd.DataFrame:
        total_columns = list(total_cols or DEFAULT_TOTAL_COLS)
        validated = require_dataframe(frame, context=context)
        require_unique_columns(validated.columns, context=context)
        require_columns(
            validated,
            required_columns=[TOTAL_GENE_COLUMN, *total_columns],
            context=context,
        )
        require_non_null_values(
            validated,
            columns=[TOTAL_GENE_COLUMN],
            context=context,
        )
        validated = coerce_numeric_columns(
            validated,
            columns=total_columns,
            context=context,
        )
        return validated


class PhosphoInputSchema:
    """Validate raw phosphoproteome input tables."""

    @classmethod
    def validate(
        cls,
        frame: pd.DataFrame,
        *,
        phospho_cols: Sequence[str] | None = None,
        context: str = "phospho input table",
    ) -> pd.DataFrame:
        phospho_columns = list(phospho_cols or DEFAULT_PHOSPHO_COLS)
        validated = require_dataframe(frame, context=context)
        require_unique_columns(validated.columns, context=context)
        require_columns(
            validated,
            required_columns=[
                *PHOSPHO_REQUIRED_METADATA_COLUMNS,
                *phospho_columns,
            ],
            context=context,
        )
        require_non_null_values(
            validated,
            columns=[PHOSPHO_UID_COLUMN, PHOSPHO_GENE_COLUMN, GENE_P_SITE_COLUMN],
            context=context,
        )
        validated = coerce_numeric_columns(
            validated,
            columns=[LOCALIZATION_PROB_COLUMN, *phospho_columns],
            context=context,
        )
        require_value_range(
            validated,
            columns=[LOCALIZATION_PROB_COLUMN],
            minimum=0.0,
            maximum=1.0,
            context=context,
        )
        require_splitable_gene_p_site(
            validated[GENE_P_SITE_COLUMN],
            context=context,
            column_name=GENE_P_SITE_COLUMN,
        )
        return validated


class PredMatSchema:
    """Validate prediction matrices used for kinase activity analysis."""

    @classmethod
    def validate(
        cls,
        frame: pd.DataFrame,
        *,
        context: str = "pred_mat",
    ) -> pd.DataFrame:
        validated = require_dataframe(frame, context=context)
        require_unique_columns(validated.columns, context=context)
        if validated.shape[0] == 0:
            msg = f"{context} must contain at least one row"
            raise TableSchemaError(msg)
        if validated.shape[1] == 0:
            msg = f"{context} must contain at least one kinase column"
            raise TableSchemaError(msg)
        validated = coerce_numeric_columns(
            validated,
            columns=list(validated.columns),
            context=context,
        )
        require_unique_index(validated, context=context)
        require_non_null_index(validated, context=context)
        require_value_range(
            validated,
            columns=list(validated.columns),
            minimum=0.0,
            maximum=1.0,
            context=context,
        )
        return validated


class PredictionScoreMatrixSchema:
    """Validate combined score matrices used by kinase prediction workflows."""

    @classmethod
    def validate(
        cls,
        frame: pd.DataFrame,
        *,
        context: str = "combined_scores",
    ) -> pd.DataFrame:
        validated = require_dataframe(frame, context=context)
        require_unique_columns(validated.columns, context=context)
        require_non_null_column_names(validated.columns, context=context)
        if validated.shape[0] == 0:
            msg = f"{context} must contain at least one phosphosite row"
            raise TableSchemaError(msg)
        if validated.shape[1] == 0:
            msg = f"{context} must contain at least one kinase column"
            raise TableSchemaError(msg)
        validated = coerce_numeric_columns(
            validated,
            columns=list(validated.columns),
            context=context,
        )
        require_unique_index(validated, context=context)
        require_non_null_index(validated, context=context)
        require_finite_numeric_values(
            validated,
            columns=list(validated.columns),
            context=context,
        )
        return validated


class SiteMatrixSourceSchema:
    """Validate corrected phosphosite rows before site-matrix construction."""

    @classmethod
    def validate(
        cls,
        frame: pd.DataFrame,
        *,
        gene_p_site_col: str,
        sequence_col: str,
        value_cols: Sequence[str],
        context: str = "site-matrix source table",
    ) -> pd.DataFrame:
        validated = require_dataframe(frame, context=context)
        require_unique_columns(validated.columns, context=context)
        requested_value_cols = list(value_cols)
        if not requested_value_cols:
            msg = f"{context} must declare at least one numeric value column"
            raise TableSchemaError(msg)
        require_columns(
            validated,
            required_columns=[gene_p_site_col, sequence_col, *requested_value_cols],
            context=context,
        )
        require_non_null_values(
            validated,
            columns=[gene_p_site_col],
            context=context,
        )
        validated = coerce_numeric_columns(
            validated,
            columns=requested_value_cols,
            context=context,
        )
        require_splitable_gene_p_site(
            validated[gene_p_site_col],
            context=context,
            column_name=gene_p_site_col,
        )
        return validated


class SiteMatrixSchema:
    """Validate phosphosite matrices used by downstream analysis and workflows."""

    @classmethod
    def validate(
        cls,
        frame: pd.DataFrame,
        *,
        context: str = "site matrix",
    ) -> pd.DataFrame:
        validated = require_dataframe(frame, context=context)
        require_unique_columns(validated.columns, context=context)
        if validated.shape[0] == 0:
            msg = f"{context} must contain at least one phosphosite row"
            raise TableSchemaError(msg)
        if validated.shape[1] == 0:
            msg = f"{context} must contain at least one numeric value column"
            raise TableSchemaError(msg)
        validated = coerce_numeric_columns(
            validated,
            columns=list(validated.columns),
            context=context,
        )
        require_unique_index(validated, context=context)
        require_non_null_index(validated, context=context)
        require_finite_numeric_values(
            validated,
            columns=list(validated.columns),
            context=context,
        )
        return validated


def normalize_pred_mat_input(
    pred_mat: pd.DataFrame | PredMatResult | None,
) -> pd.DataFrame | None:
    """Normalize public predMat inputs to the internal DataFrame contract."""

    from ..prediction.models import PredMatResult

    if isinstance(pred_mat, PredMatResult):
        return pred_mat.to_frame(copy=False)
    return pred_mat


__all__ = [
    "PhosphoInputSchema",
    "PredMatSchema",
    "PredictionScoreMatrixSchema",
    "SiteMatrixSchema",
    "SiteMatrixSourceSchema",
    "TotalInputSchema",
    "normalize_pred_mat_input",
]
