from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from .matrices import build_site_matrix
from .validation.tables import SiteMatrixSchema


@dataclass(slots=True)
class SiteMatrixResult:
    phosr_input: pd.DataFrame
    matrix: pd.DataFrame
    sequences: pd.Series
    row_drop_stats: dict[str, int]


class SiteMatrixBuilder:
    """Build PhosR-style site matrices from corrected phosphosite rows."""

    def __init__(self, *, value_cols: Sequence[str]) -> None:
        self.value_cols = list(value_cols)

    def build(
        self,
        corrected_df: pd.DataFrame,
        *,
        gene_p_site_col: str = "gene_p_site",
        sequence_col: str = "centralized_sequence",
    ) -> SiteMatrixResult:
        phosr_input, matrix, sequences = build_site_matrix(
            df=corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
            value_cols=self.value_cols,
        )
        matrix = SiteMatrixSchema.validate(matrix, context="site matrix")
        row_drop_stats = dict(phosr_input.attrs.get("row_drop_stats", {}))
        return SiteMatrixResult(
            phosr_input=phosr_input,
            matrix=matrix,
            sequences=sequences,
            row_drop_stats=row_drop_stats,
        )
