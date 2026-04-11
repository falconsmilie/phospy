from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal, overload

import pandas as pd

from ..validation.compatibility import (
    ProteinCorrectionMatchSummary,
    validate_core_column_alignment,
    validate_protein_correction_inputs,
)
from ..validation.schema.frames import require_columns, require_numeric_columns
from ..validation.values.collections import resolve_required_columns
from ..validation.values.identifiers import normalize_identifier_series


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


@dataclass(slots=True)
class ProteinCorrectionResult:
    """Protein-correction result bundle with corrected phosphosite rows.

    The corrected table is mutable pandas state; this wrapper does not imply
    immutable snapshot semantics.
    """

    corrected: pd.DataFrame
    summary: ProteinCorrectionSummary


def _resolve_correction_columns(
    *,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str] | None,
    output_prefix: str,
    context: str,
) -> tuple[list[str], list[str], list[str]]:
    resolved_phospho_cols = resolve_required_columns(
        phospho_cols,
        argument_name="phospho_cols",
        context=context,
    )
    resolved_protein_cols = resolve_required_columns(
        protein_cols,
        argument_name="protein_cols",
        context=context,
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
        context=f"{context} inputs",
    )
    return resolved_phospho_cols, resolved_protein_cols, resolved_corrected_cols


def _build_correction_summary(
    *,
    match_summary: ProteinCorrectionMatchSummary,
    phospho_gene_col: str,
    total_gene_col: str,
) -> ProteinCorrectionSummary:
    return ProteinCorrectionSummary(
        input_rows=match_summary.input_rows,
        matched_rows=match_summary.matched_rows,
        unmatched_rows=match_summary.unmatched_rows,
        unmatched_fraction=match_summary.unmatched_fraction,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        unmatched_gene_preview=match_summary.unmatched_gene_preview,
    )


def _run_protein_correction_from_numeric_inputs(
    *,
    df_phospho: pd.DataFrame,
    df_total: pd.DataFrame,
    phospho_gene_col: str,
    total_gene_col: str,
    phospho_cols: Sequence[str],
    protein_cols: Sequence[str],
    corrected_cols: Sequence[str],
    max_unmatched_fraction: float,
    context: str,
) -> tuple[pd.DataFrame, ProteinCorrectionMatchSummary]:
    match_summary = validate_protein_correction_inputs(
        df_phospho,
        df_total,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        phospho_cols=phospho_cols,
        protein_cols=protein_cols,
        max_unmatched_fraction=max_unmatched_fraction,
        context=f"{context} inputs",
    )

    phospho_join_key = normalize_identifier_series(df_phospho[phospho_gene_col])
    total_join_key = normalize_identifier_series(df_total[total_gene_col])

    total_lookup = df_total.loc[:, list(protein_cols)].set_index(total_join_key)
    matched_mask = phospho_join_key.isin(total_lookup.index)
    matched_phospho = df_phospho.loc[matched_mask].copy()
    matched_total = total_lookup.loc[phospho_join_key.loc[matched_mask].tolist()]

    for protein_col in protein_cols:
        matched_phospho[protein_col] = matched_total[protein_col].to_numpy(copy=False)

    for corrected_col, p_col, t_col in zip(
        corrected_cols,
        phospho_cols,
        protein_cols,
        strict=True,
    ):
        matched_phospho[corrected_col] = matched_phospho[p_col].to_numpy(
            copy=False
        ) - matched_total[t_col].to_numpy(copy=False)

    return matched_phospho, match_summary


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
    context = "correct_phospho_to_protein()"
    resolved_phospho_cols, resolved_protein_cols, resolved_corrected_cols = (
        _resolve_correction_columns(
            phospho_cols=phospho_cols,
            protein_cols=protein_cols,
            corrected_cols=corrected_cols,
            output_prefix=output_prefix,
            context=context,
        )
    )

    require_columns(
        df_phospho,
        required_columns=[phospho_gene_col, *resolved_phospho_cols],
        context=f"{context} phospho input",
    )
    require_columns(
        df_total,
        required_columns=[total_gene_col, *resolved_protein_cols],
        context=f"{context} total input",
    )

    phospho_work = df_phospho.copy()
    total_work = df_total.copy()
    require_numeric_columns(
        phospho_work,
        columns=resolved_phospho_cols,
        context=f"{context} phospho input",
    )
    require_numeric_columns(
        total_work,
        columns=resolved_protein_cols,
        context=f"{context} total input",
    )

    corrected, match_summary = _run_protein_correction_from_numeric_inputs(
        df_phospho=phospho_work,
        df_total=total_work,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        phospho_cols=resolved_phospho_cols,
        protein_cols=resolved_protein_cols,
        corrected_cols=resolved_corrected_cols,
        max_unmatched_fraction=max_unmatched_fraction,
        context=context,
    )

    if not return_summary:
        return corrected

    return ProteinCorrectionResult(
        corrected=corrected,
        summary=_build_correction_summary(
            match_summary=match_summary,
            phospho_gene_col=phospho_gene_col,
            total_gene_col=total_gene_col,
        ),
    )


def run_protein_correction_owned(
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
) -> pd.DataFrame:
    """Correct trusted numeric phospho and total tables without copying inputs.

    This internal helper assumes the caller already owns the input frames and
    has prepared their numeric columns. It still validates column alignment and
    phospho/protein matching coverage, but it does not take additional full-frame
    defensive copies.
    """

    context = "correct_phospho_to_protein()"
    resolved_phospho_cols, resolved_protein_cols, resolved_corrected_cols = (
        _resolve_correction_columns(
            phospho_cols=phospho_cols,
            protein_cols=protein_cols,
            corrected_cols=corrected_cols,
            output_prefix=output_prefix,
            context=context,
        )
    )
    require_columns(
        df_phospho,
        required_columns=[phospho_gene_col, *resolved_phospho_cols],
        context=f"{context} phospho input",
    )
    require_columns(
        df_total,
        required_columns=[total_gene_col, *resolved_protein_cols],
        context=f"{context} total input",
    )
    corrected, _ = _run_protein_correction_from_numeric_inputs(
        df_phospho=df_phospho,
        df_total=df_total,
        phospho_gene_col=phospho_gene_col,
        total_gene_col=total_gene_col,
        phospho_cols=resolved_phospho_cols,
        protein_cols=resolved_protein_cols,
        corrected_cols=resolved_corrected_cols,
        max_unmatched_fraction=max_unmatched_fraction,
        context=context,
    )
    return corrected


__all__ = [
    "ProteinCorrectionResult",
    "ProteinCorrectionSummary",
    "run_protein_correction",
    "run_protein_correction_owned",
]
