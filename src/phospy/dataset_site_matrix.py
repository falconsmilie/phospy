from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from .constants import CENTRALIZED_SEQUENCE_COLUMN, GENE_P_SITE_COLUMN
from .dataset_schema import DatasetSchema
from .site_matrix_builder import SiteMatrixBuilder, SiteMatrixResult


@dataclass(frozen=True, slots=True)
class DatasetSiteMatrix:
    """Bound site-matrix facade for a validated phosphoproteomics dataset."""

    schema: DatasetSchema

    def _builder(self) -> SiteMatrixBuilder:
        return SiteMatrixBuilder(value_cols=self.schema.corrected_cols)

    def build(
        self,
        corrected_df: pd.DataFrame,
        *,
        gene_p_site_col: str = GENE_P_SITE_COLUMN,
        sequence_col: str = CENTRALIZED_SEQUENCE_COLUMN,
    ) -> SiteMatrixResult:
        return self._builder().build(
            corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
        )
