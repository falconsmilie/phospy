from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from ...errors import InputCompatibilityError
from ..schema.frames import require_columns
from ..values.identifiers import normalize_identifier_series
from ..values.numeric import validate_fraction


@dataclass(frozen=True, slots=True)
class ProteinCorrectionMatchSummary:
    """Describe phosphosite-to-protein matching before correction."""

    input_rows: int
    matched_rows: int
    unmatched_rows: int
    unmatched_fraction: float
    unmatched_gene_preview: tuple[str, ...]


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

    resolved_max_unmatched_fraction = validate_fraction(
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

    require_columns(
        phospho_df,
        required_columns=[phospho_gene_col, *phospho_cols],
        context=f"{context} phospho input",
    )
    require_columns(
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


__all__ = ["ProteinCorrectionMatchSummary", "validate_protein_correction_inputs"]
