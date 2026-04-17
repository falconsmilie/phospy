from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from ..errors import TableSchemaError
from ..internal.constants import (
    CENTRALIZED_SEQUENCE_COLUMN,
    GENE_P_SITE_COLUMN,
    ROW_DROP_STATS_ATTR,
)
from ..validation.schema.tables import SiteMatrixSchema
from .matrices import SiteMatrixPolicy, build_site_matrix, format_row_drop_diagnostics


@dataclass(slots=True)
class SiteMatrixResult:
    """Derived site-matrix tables built from corrected phosphosite rows."""

    phosr_input: pd.DataFrame
    matrix: pd.DataFrame
    sequences: pd.Series
    row_drop_stats: dict[str, int | str]


class SiteMatrixBuilder:
    """Build PhosR-style site matrices from corrected phosphosite rows.

    ``build()`` is the public defensive boundary. ``build_owned()`` is the
    internal fast path for callers that already own the corrected table.
    """

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
        return self._build(
            corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
            policy=policy,
            copy_frame=True,
        )

    def build_owned(
        self,
        corrected_df: pd.DataFrame,
        *,
        gene_p_site_col: str = GENE_P_SITE_COLUMN,
        sequence_col: str = CENTRALIZED_SEQUENCE_COLUMN,
        policy: SiteMatrixPolicy | None = None,
    ) -> SiteMatrixResult:
        """Build a site matrix from an already-owned corrected table without copying."""
        return self._build(
            corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
            policy=policy,
            copy_frame=False,
        )

    def _build(
        self,
        corrected_df: pd.DataFrame,
        *,
        gene_p_site_col: str,
        sequence_col: str,
        policy: SiteMatrixPolicy | None,
        copy_frame: bool,
    ) -> SiteMatrixResult:
        phosr_input, matrix, sequences = build_site_matrix(
            df=corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
            value_cols=self.value_cols,
            policy=policy,
            copy_frame=copy_frame,
        )
        row_drop_stats = dict(phosr_input.attrs.get(ROW_DROP_STATS_ATTR, {}))
        if matrix.empty:
            diagnostics = format_row_drop_diagnostics(row_drop_stats)
            raise TableSchemaError(
                f"site matrix must contain at least one phosphosite row; {diagnostics}"
            )
        matrix = SiteMatrixSchema.validate(
            matrix,
            context="site matrix",
            copy_frame=False,
        )
        return SiteMatrixResult(
            phosr_input=phosr_input,
            matrix=matrix,
            sequences=sequences,
            row_drop_stats=row_drop_stats,
        )
