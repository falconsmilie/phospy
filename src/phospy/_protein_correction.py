from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Literal, overload

import pandas as pd

from ._preprocessing_primitives import _require_columns, _require_numeric_series
from .validation.compatibility import (
    validate_core_column_alignment,
    validate_protein_correction_inputs,
)
from .validation.errors import InputCompatibilityError, PhospyValidationError
from .validation.normalization import normalize_identifier_series


@dataclass(frozen=True, slots=True)
class ProteinCorrectionSummary:
    """Describe a phospho-to-protein correction pass."""

    input_rows: int
    matched_rows: int
    unmatched_rows: int
    unmatched_fraction: float
    phospho_gene_col: str
    total_gene_col: str
    unmatched_gene_preview: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ProteinCorrectionResult:
    """Corrected phosphosite rows together with protein-match metadata."""

    corrected: pd.DataFrame
    summary: ProteinCorrectionSummary


def _resolve_required_columns(
    columns: Iterable[str],
    *,
    argument_name: str,
    context: str,
) -> list[str]:
    resolved_columns = list(columns)
    if not resolved_columns:
        raise PhospyValidationError(
            f"{context} requires at least one column name in '{argument_name}'"
        )
    return resolved_columns


def _require_numeric_columns(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    context: str,
) -> None:
    for column in columns:
        df[column] = _require_numeric_series(
            df[column],
            column=column,
            context=context,
        )


@overload
def run_protein_correction(
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    output_prefix: str = "phospho_corrected_",
    *,
    max_unmatched_fraction: float = 1.0,
    return_summary: Literal[False] = False,
) -> pd.DataFrame: ...


@overload
def run_protein_correction(
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    output_prefix: str = "phospho_corrected_",
    *,
    max_unmatched_fraction: float = 1.0,
    return_summary: Literal[True],
) -> ProteinCorrectionResult: ...


def run_protein_correction(
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str] | None = None,
    output_prefix: str = "phospho_corrected_",
    *,
    max_unmatched_fraction: float = 1.0,
    return_summary: bool = False,
) -> pd.DataFrame | ProteinCorrectionResult:
    resolved_phospho_cols = _resolve_required_columns(
        phospho_cols,
        argument_name="phospho_cols",
        context="correct_phospho_to_protein()",
    )
    resolved_protein_cols = _resolve_required_columns(
        protein_cols,
        argument_name="protein_cols",
        context="correct_phospho_to_protein()",
    )
    resolved_corrected_cols = (
        list(corrected_cols)
        if corrected_cols is not None
        else [
            f"{output_prefix}{idx}" for idx in range(1, len(resolved_phospho_cols) + 1)
        ]
    )

    validate_core_column_alignment(
        resolved_protein_cols,
        resolved_phospho_cols,
        resolved_corrected_cols,
        context="correct_phospho_to_protein() inputs",
    )
    _require_columns(
        df_phospho,
        required_columns=[phospho_gene_col, *resolved_phospho_cols],
        context="correct_phospho_to_protein() phospho input",
    )
    _require_columns(
        df_total,
        required_columns=[total_gene_col, *resolved_protein_cols],
        context="correct_phospho_to_protein() total input",
    )

    phospho_join_col = "__phospy_normalized_phospho_gene_key"
    total_join_col = "__phospy_normalized_total_gene_key"

    phospho_work = df_phospho.copy()
    total_work = df_total.copy()
    _require_numeric_columns(
        phospho_work,
        columns=resolved_phospho_cols,
        context="correct_phospho_to_protein() phospho input",
    )
    _require_numeric_columns(
        total_work,
        columns=resolved_protein_cols,
        context="correct_phospho_to_protein() total input",
    )
    match_summary = validate_protein_correction_inputs(
        phospho_work,
        total_work,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        phospho_cols=resolved_phospho_cols,
        protein_cols=resolved_protein_cols,
        max_unmatched_fraction=max_unmatched_fraction,
        context="correct_phospho_to_protein() inputs",
    )

    phospho_work[phospho_join_col] = normalize_identifier_series(
        phospho_work[phospho_gene_col]
    )
    total_work[total_join_col] = normalize_identifier_series(total_work[total_gene_col])

    if total_work[total_join_col].duplicated().any():
        msg = (
            f"{total_gene_col} must be unique before protein correction to avoid "
            "duplicating phosphosite rows during the merge"
        )
        raise InputCompatibilityError(msg)

    merged = phospho_work.merge(
        total_work[[total_join_col, total_gene_col, *resolved_protein_cols]],
        left_on=phospho_join_col,
        right_on=total_join_col,
        how="inner",
    )

    drop_columns: list[str] = [phospho_join_col, total_join_col]
    if total_gene_col != phospho_gene_col and total_gene_col in merged.columns:
        drop_columns.append(total_gene_col)
    merged = merged.drop(columns=drop_columns, errors="ignore")

    for corrected_col, p_col, t_col in zip(
        resolved_corrected_cols,
        resolved_phospho_cols,
        resolved_protein_cols,
        strict=True,
    ):
        merged[corrected_col] = merged[p_col] - merged[t_col]

    if not return_summary:
        return merged

    summary = ProteinCorrectionSummary(
        input_rows=match_summary.input_rows,
        matched_rows=match_summary.matched_rows,
        unmatched_rows=match_summary.unmatched_rows,
        unmatched_fraction=match_summary.unmatched_fraction,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        unmatched_gene_preview=match_summary.unmatched_gene_preview,
    )
    return ProteinCorrectionResult(corrected=merged, summary=summary)


__all__ = [
    "ProteinCorrectionResult",
    "ProteinCorrectionSummary",
    "run_protein_correction",
]
