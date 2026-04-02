from __future__ import annotations

from collections.abc import Iterable, Sequence

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
        validated = _ensure_dataframe(frame, context=context)
        _ensure_unique_columns(validated.columns, context=context)
        _ensure_required_columns(
            validated, [TOTAL_GENE_COLUMN, *total_columns], context=context
        )
        _ensure_non_null(validated, [TOTAL_GENE_COLUMN], context=context)
        validated = _coerce_numeric_columns(
            validated,
            total_columns,
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
        validated = _ensure_dataframe(frame, context=context)
        _ensure_unique_columns(validated.columns, context=context)
        _ensure_required_columns(
            validated,
            [
                *PHOSPHO_REQUIRED_METADATA_COLUMNS,
                *phospho_columns,
            ],
            context=context,
        )
        _ensure_non_null(
            validated,
            [PHOSPHO_UID_COLUMN, PHOSPHO_GENE_COLUMN, GENE_P_SITE_COLUMN],
            context=context,
        )
        validated = _coerce_numeric_columns(
            validated,
            [LOCALIZATION_PROB_COLUMN, *phospho_columns],
            context=context,
        )
        _ensure_value_range(
            validated,
            [LOCALIZATION_PROB_COLUMN],
            minimum=0.0,
            maximum=1.0,
            context=context,
        )
        _ensure_splitable_gene_p_site(
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
        validated = _ensure_dataframe(frame, context=context)
        _ensure_unique_columns(validated.columns, context=context)
        if validated.empty:
            msg = f"{context} must contain at least one row"
            raise TableSchemaError(msg)
        if validated.shape[1] == 0:
            msg = f"{context} must contain at least one kinase column"
            raise TableSchemaError(msg)
        validated = _coerce_numeric_columns(
            validated,
            list(validated.columns),
            context=context,
        )
        _ensure_unique_index(validated, context=context)
        _ensure_non_null_index(validated, context=context)
        _ensure_value_range(
            validated,
            list(validated.columns),
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
        validated = _ensure_dataframe(frame, context=context)
        _ensure_unique_columns(validated.columns, context=context)
        _ensure_non_null_columns(validated.columns, context=context)
        if validated.empty:
            msg = f"{context} must contain at least one phosphosite row"
            raise TableSchemaError(msg)
        if validated.shape[1] == 0:
            msg = f"{context} must contain at least one kinase column"
            raise TableSchemaError(msg)
        validated = _coerce_numeric_columns(
            validated,
            list(validated.columns),
            context=context,
        )
        _ensure_unique_index(validated, context=context)
        _ensure_non_null_index(validated, context=context)
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
        validated = _ensure_dataframe(frame, context=context)
        _ensure_unique_columns(validated.columns, context=context)
        requested_value_cols = list(value_cols)
        if not requested_value_cols:
            msg = f"{context} must declare at least one numeric value column"
            raise TableSchemaError(msg)
        _ensure_required_columns(
            validated,
            [gene_p_site_col, sequence_col, *requested_value_cols],
            context=context,
        )
        _ensure_non_null(validated, [gene_p_site_col], context=context)
        validated = _coerce_numeric_columns(
            validated,
            requested_value_cols,
            context=context,
        )
        _ensure_splitable_gene_p_site(
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
        validated = _ensure_dataframe(frame, context=context)
        _ensure_unique_columns(validated.columns, context=context)
        if validated.empty:
            msg = f"{context} must contain at least one phosphosite row"
            raise TableSchemaError(msg)
        if validated.shape[1] == 0:
            msg = f"{context} must contain at least one numeric value column"
            raise TableSchemaError(msg)
        validated = _coerce_numeric_columns(
            validated,
            list(validated.columns),
            context=context,
        )
        _ensure_unique_index(validated, context=context)
        _ensure_non_null_index(validated, context=context)
        return validated


def _ensure_dataframe(frame: pd.DataFrame, *, context: str) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        msg = f"{context} must be a pandas DataFrame"
        raise TableSchemaError(msg)
    return frame.copy()


def _ensure_required_columns(
    frame: pd.DataFrame,
    required_columns: Sequence[str],
    *,
    context: str,
) -> None:
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        missing_str = ", ".join(missing)
        msg = f"{context} is missing required columns: {missing_str}"
        raise TableSchemaError(msg)


def _ensure_unique_columns(columns: Iterable[object], *, context: str) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    for column in columns:
        value = str(column)
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    if duplicates:
        duplicates_str = ", ".join(duplicates)
        msg = f"{context} contains duplicate column names: {duplicates_str}"
        raise TableSchemaError(msg)


def _ensure_non_null_columns(columns: Iterable[object], *, context: str) -> None:
    null_like: list[str] = []
    for column in columns:
        if pd.isna(column):
            null_like.append("<null>")
    if null_like:
        msg = f"{context} contains null column names"
        raise TableSchemaError(msg)


def _ensure_non_null(
    frame: pd.DataFrame,
    columns: Sequence[str],
    *,
    context: str,
) -> None:
    failures: list[str] = []
    for column in columns:
        if frame[column].isna().any():
            failures.append(column)
    if failures:
        failures_str = ", ".join(failures)
        msg = f"{context} contains null values in required columns: {failures_str}"
        raise TableSchemaError(msg)


def _coerce_numeric_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
    *,
    context: str,
) -> pd.DataFrame:
    validated = frame.copy()
    failures: list[str] = []
    for column in columns:
        converted = pd.to_numeric(validated[column], errors="coerce")
        invalid_mask = validated[column].notna() & converted.isna()
        if invalid_mask.any():
            sample_values = validated.loc[invalid_mask, column].astype(str).unique()[:3]
            sample_preview = ", ".join(sample_values)
            failures.append(f"{column} ({sample_preview})")
        validated[column] = converted
    if failures:
        failures_str = "; ".join(failures)
        msg = (
            f"{context} contains non-numeric values in numeric columns: {failures_str}"
        )
        raise TableSchemaError(msg)
    return validated


def _ensure_value_range(
    frame: pd.DataFrame,
    columns: Sequence[str],
    *,
    minimum: float,
    maximum: float,
    context: str,
) -> None:
    failures: list[str] = []
    for column in columns:
        series = frame[column]
        out_of_range = series.notna() & ((series < minimum) | (series > maximum))
        if out_of_range.any():
            failures.append(column)
    if failures:
        failures_str = ", ".join(failures)
        msg = (
            f"{context} contains values outside the allowed range "
            f"[{minimum}, {maximum}] in columns: {failures_str}"
        )
        raise TableSchemaError(msg)


def _ensure_splitable_gene_p_site(
    series: pd.Series,
    *,
    context: str,
    column_name: str = GENE_P_SITE_COLUMN,
) -> None:
    normalized = series.astype("string")
    split_columns = normalized.str.split("_", n=1, expand=True)
    underscore_count = normalized.str.count("_")
    if split_columns.shape[1] < 2:
        invalid_mask = pd.Series(True, index=series.index)
    else:
        invalid_mask = (
            normalized.isna()
            | (underscore_count != 1)
            | split_columns[0].isna()
            | split_columns[1].isna()
            | (split_columns[0].str.strip().str.len() == 0)
            | (split_columns[1].str.strip().str.len() == 0)
        )
    if invalid_mask.any():
        sample_values = series.loc[invalid_mask].astype(str).unique()[:3]
        sample_preview = ", ".join(sample_values)
        msg = (
            f"{context} contains malformed {column_name} values that cannot be split "
            f"into non-empty gene and site parts using a single underscore: "
            f"{sample_preview}"
        )
        raise TableSchemaError(msg)


def _ensure_unique_index(frame: pd.DataFrame, *, context: str) -> None:
    duplicated = frame.index.duplicated(keep=False)
    if duplicated.any():
        duplicates = frame.index[duplicated].astype(str).unique().tolist()[:5]
        duplicates_str = ", ".join(duplicates)
        msg = f"{context} contains duplicate index entries: {duplicates_str}"
        raise TableSchemaError(msg)


def _ensure_non_null_index(frame: pd.DataFrame, *, context: str) -> None:
    if frame.index.isna().any():
        msg = f"{context} contains null index entries"
        raise TableSchemaError(msg)
