from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from .constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from .core_processing import (
    CorePreprocessingConfig,
    CoreProcessingResult,
    CoreProcessor,
)
from .dataset_schema import DatasetSchema
from .site_matrix_builder import SiteMatrixBuilder, SiteMatrixResult


@dataclass(frozen=True, slots=True)
class DatasetPreprocessing:
    """Bound preprocessing facade for a validated phosphoproteomics dataset."""

    total_df: pd.DataFrame
    phospho_df: pd.DataFrame
    schema: DatasetSchema
    comparisons: tuple[ComparisonSpec, ...] | None = None

    def __init__(
        self,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> None:
        object.__setattr__(self, "total_df", total_df)
        object.__setattr__(self, "phospho_df", phospho_df)
        object.__setattr__(self, "schema", schema)
        object.__setattr__(
            self,
            "comparisons",
            tuple(comparisons) if comparisons is not None else None,
        )

    def _core_processor(self) -> CoreProcessor:
        return CoreProcessor(
            schema=self.schema,
            comparisons=self.comparisons,
        )

    def _site_matrix_builder(self) -> SiteMatrixBuilder:
        return SiteMatrixBuilder(value_cols=self.schema.corrected_cols)

    def prepare_total(
        self,
        gene_col: str = "genes",
        sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        min_observed: int = 4,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        return self._core_processor().prepare_total(
            self.total_df,
            gene_col=gene_col,
            sentinel=sentinel,
            min_observed=min_observed,
        )

    def prepare_phospho(
        self,
        gene_col: str = "gene_names",
        site_col: str = "gene_p_site",
        localization_col: str = "localization_prob",
        localization_threshold: float = 0.75,
        sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        min_observed: int = 4,
    ) -> pd.DataFrame:
        return self._core_processor().prepare_phospho(
            self.phospho_df,
            gene_col=gene_col,
            site_col=site_col,
            localization_col=localization_col,
            localization_threshold=localization_threshold,
            sentinel=sentinel,
            min_observed=min_observed,
        )

    def correct_to_protein(
        self,
        phospho_df: pd.DataFrame,
        total_df: pd.DataFrame,
        phospho_gene_col: str = "gene_names",
        total_gene_col: str = "genes",
        max_unmatched_fraction: float = 0.0,
    ) -> pd.DataFrame:
        return self._core_processor().correct_to_protein(
            phospho_df,
            total_df,
            phospho_gene_col=phospho_gene_col,
            total_gene_col=total_gene_col,
            max_unmatched_fraction=max_unmatched_fraction,
        )

    def add_pairwise_comparisons(
        self,
        corrected_df: pd.DataFrame,
        output_prefix: str = "p_",
    ) -> pd.DataFrame:
        return self._core_processor().add_pairwise_comparisons(
            corrected_df,
            output_prefix=output_prefix,
        )

    def build_site_matrix(
        self,
        corrected_df: pd.DataFrame,
        gene_p_site_col: str = "gene_p_site",
        sequence_col: str = "centralized_sequence",
    ) -> SiteMatrixResult:
        return self._site_matrix_builder().build(
            corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
        )

    def run(
        self,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        config: CorePreprocessingConfig | None = None,
    ) -> CoreProcessingResult:
        resolved = config or CorePreprocessingConfig(
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            total_sentinel=float(total_sentinel),
            phospho_sentinel=float(phospho_sentinel),
            max_unmatched_fraction=max_unmatched_fraction,
        )
        return self._core_processor().process(
            self.total_df,
            self.phospho_df,
            config=resolved,
        )
