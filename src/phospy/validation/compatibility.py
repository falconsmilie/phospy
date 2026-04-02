from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from numbers import Real

import pandas as pd

from .errors import InputCompatibilityError, PhospyValidationError, TableSchemaError
from .normalization import normalize_identifier_series


@dataclass(frozen=True, slots=True)
class ProteinCorrectionMatchSummary:
    """Describe phosphosite-to-protein matching before correction."""

    input_rows: int
    matched_rows: int
    unmatched_rows: int
    unmatched_fraction: float
    unmatched_gene_preview: tuple[str, ...]


def _validate_fraction(
    value: float,
    *,
    name: str,
) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )

    resolved = float(value)
    if (
        not 0.0 <= resolved <= 1.0
        or resolved != resolved
        or resolved in {float("inf"), float("-inf")}
    ):
        raise PhospyValidationError(
            f"{name} must be a finite numeric value between 0 and 1"
        )
    return resolved


def _require_columns(
    frame: pd.DataFrame,
    *,
    required_columns: Sequence[str],
    context: str,
) -> None:
    missing_columns = [
        column for column in required_columns if column not in frame.columns
    ]
    if missing_columns:
        joined_columns = ", ".join(missing_columns)
        raise TableSchemaError(
            f"{context} is missing required columns: {joined_columns}"
        )


def validate_core_column_alignment(
    total_cols: Sequence[str],
    phospho_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    *,
    context: str = "Core preprocessing inputs",
) -> None:
    """Validate that paired preprocessing column groups align by width."""

    if len(total_cols) != len(phospho_cols):
        msg = f"{context} require the same number of total and phospho value columns"
        raise InputCompatibilityError(msg)
    if corrected_cols is not None and len(corrected_cols) != len(total_cols):
        msg = (
            f"{context} require corrected value columns to align with total and "
            "phospho value columns"
        )
        raise InputCompatibilityError(msg)


def validate_protein_correction_inputs(
    phospho_df: pd.DataFrame,
    total_df: pd.DataFrame,
    *,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    max_unmatched_fraction: float = 0.0,
    context: str = "Protein correction inputs",
) -> ProteinCorrectionMatchSummary:
    """Validate phosphosite/protein correction inputs and matching coverage."""

    resolved_max_unmatched_fraction = _validate_fraction(
        max_unmatched_fraction,
        name="max_unmatched_fraction",
    )
    if len(phospho_cols) != len(protein_cols):
        msg = f"{context} require the same number of phospho and protein columns"
        raise InputCompatibilityError(msg)
    if phospho_df.empty:
        msg = f"{context} contain no phosphosite rows after filtering"
        raise InputCompatibilityError(msg)
    if total_df.empty:
        msg = f"{context} contain no protein rows after filtering"
        raise InputCompatibilityError(msg)

    _require_columns(
        phospho_df,
        required_columns=[phospho_gene_col, *phospho_cols],
        context=f"{context} phospho input",
    )
    _require_columns(
        total_df,
        required_columns=[total_gene_col, *protein_cols],
        context=f"{context} total input",
    )

    total_gene_series = normalize_identifier_series(total_df[total_gene_col])
    if total_gene_series.duplicated().any():
        msg = (
            f"{context}: {total_gene_col} must be unique before protein "
            "correction to avoid duplicating phosphosite rows"
        )
        raise InputCompatibilityError(msg)

    phospho_genes = normalize_identifier_series(phospho_df[phospho_gene_col])
    total_gene_values = set(total_gene_series)
    matched_mask = phospho_genes.isin(total_gene_values)
    matched_rows = int(matched_mask.sum())

    if matched_rows == 0:
        msg = (
            f"{context} have no overlapping gene identifiers between "
            f"{phospho_gene_col} and {total_gene_col}"
        )
        raise InputCompatibilityError(msg)

    input_rows = int(len(phospho_df))
    unmatched_rows = input_rows - matched_rows
    unmatched_fraction = unmatched_rows / input_rows
    unmatched_genes = pd.unique(phospho_genes.loc[~matched_mask].dropna())
    unmatched_gene_preview = tuple(str(gene) for gene in unmatched_genes[:5])

    if unmatched_rows > 0 and unmatched_fraction > resolved_max_unmatched_fraction:
        unmatched_preview = ", ".join(unmatched_gene_preview)
        percent = unmatched_fraction * 100.0
        msg = (
            f"{context} would drop {unmatched_rows} of {input_rows} phosphosite "
            f"rows ({percent:.1f}%) due to missing protein matches in "
            f"{total_gene_col}: {unmatched_preview}"
        )
        raise InputCompatibilityError(msg)

    return ProteinCorrectionMatchSummary(
        input_rows=input_rows,
        matched_rows=matched_rows,
        unmatched_rows=unmatched_rows,
        unmatched_fraction=unmatched_fraction,
        unmatched_gene_preview=unmatched_gene_preview,
    )


def validate_pred_mat_overlap(
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
    *,
    pred_context: str = "pred_mat",
    matrix_context: str = "phospho_matrix",
    min_overlap: int = 1,
    min_fraction: float = 0.1,
) -> None:
    """Validate overlap between a prediction matrix and a phosphosite matrix."""

    overlap = pred_mat.index.intersection(phospho_matrix.index)
    overlap_count = len(overlap)
    if overlap_count == 0:
        msg = f"{pred_context} and {matrix_context} have no overlapping phosphosite IDs"
        raise InputCompatibilityError(msg)

    matrix_rows = len(phospho_matrix.index)
    overlap_fraction = overlap_count / max(matrix_rows, 1)
    if overlap_count < min_overlap or overlap_fraction < min_fraction:
        percent = overlap_fraction * 100.0
        msg = (
            f"{pred_context} and {matrix_context} have insufficient overlapping "
            f"phosphosite IDs: {overlap_count} row(s) ({percent:.1f}%)"
        )
        raise InputCompatibilityError(msg)


__all__ = [
    "ProteinCorrectionMatchSummary",
    "validate_core_column_alignment",
    "validate_pred_mat_overlap",
    "validate_protein_correction_inputs",
]
