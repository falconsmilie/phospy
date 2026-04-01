from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

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
        gene_p_site_col: str = "gene_p_site",
        sequence_col: str = "centralized_sequence",
    ) -> SiteMatrixResult:
        return self._builder().build(
            corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
        )
