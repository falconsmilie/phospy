from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from ..errors import TableSchemaError
from ..internal.constants import CENTRALIZED_SEQUENCE_COLUMN, GENE_P_SITE_COLUMN
from ..matrices import SiteMatrixPolicy, build_site_matrix, format_row_drop_diagnostics
from ..validation.schema.tables import SiteMatrixSchema


@dataclass(slots=True)
class SiteMatrixResult:
    """Derived site-matrix tables built from corrected phosphosite rows."""

    phosr_input: pd.DataFrame
    matrix: pd.DataFrame
    sequences: pd.Series
    row_drop_stats: dict[str, int | str]


class SiteMatrixBuilder:
    """Build PhosR-style site matrices from corrected phosphosite rows."""

    def __init__(self, *, value_cols: Sequence[str]) -> None:
        self.value_cols = list(value_cols)

    def build(
        self,
        corrected_df: pd.DataFrame,
        *,
        gene_p_site_col: str = GENE_P_SITE_COLUMN,
        sequence_col: str = CENTRALIZED_SEQUENCE_COLUMN,
        policy: SiteMatrixPolicy | None = None,
    ) -> SiteMatrixResult:
        phosr_input, matrix, sequences = build_site_matrix(
            df=corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
            value_cols=self.value_cols,
            policy=policy,
        )
        row_drop_stats = dict(phosr_input.attrs.get("row_drop_stats", {}))
        if matrix.empty:
            diagnostics = format_row_drop_diagnostics(row_drop_stats)
            raise TableSchemaError(
                f"site matrix must contain at least one phosphosite row; {diagnostics}"
            )
        matrix = SiteMatrixSchema.validate(matrix, context="site matrix")
        return SiteMatrixResult(
            phosr_input=phosr_input,
            matrix=matrix,
            sequences=sequences,
            row_drop_stats=row_drop_stats,
        )
