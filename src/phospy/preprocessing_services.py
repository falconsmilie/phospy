from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from .constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from .dataset_schema import DatasetSchema
from .preprocessing import (
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_localized_sites,
    filter_min_observed,
    replace_sentinel_with_nan,
)
from .validation.compatibility import validate_protein_correction_inputs


@dataclass(frozen=True, slots=True)
class TotalPreprocessor:
    """Prepare total proteome rows for downstream phosphosite correction."""

    schema: DatasetSchema

    def prepare(
        self,
        total_df: pd.DataFrame,
        *,
        gene_col: str = "genes",
        sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        min_observed: int = 4,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        total = total_df.copy()
        total[gene_col] = total[gene_col].astype("string")
        total = replace_sentinel_with_nan(
            total,
            self.schema.total_cols,
            sentinel=sentinel,
        )
        total_unique = collapse_duplicate_genes(
            total,
            gene_col=gene_col,
            value_cols=self.schema.total_cols,
        )
        total_filtered = filter_min_observed(
            total_unique,
            self.schema.total_cols,
            min_observed=min_observed,
        )
        return total_unique, total_filtered


@dataclass(frozen=True, slots=True)
class PhosphoPreprocessor:
    """Prepare phosphosite rows before protein-level correction."""

    schema: DatasetSchema

    def prepare(
        self,
        phospho_df: pd.DataFrame,
        *,
        gene_col: str = "gene_names",
        site_col: str = "gene_p_site",
        localization_col: str = "localization_prob",
        localization_threshold: float = 0.75,
        sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        min_observed: int = 4,
    ) -> pd.DataFrame:
        phospho = phospho_df.copy()
        phospho[gene_col] = phospho[gene_col].astype("string").str.upper()
        phospho[site_col] = phospho[site_col].astype("string")

        phospho = filter_localized_sites(
            phospho,
            localization_col=localization_col,
            threshold=localization_threshold,
        )
        phospho = replace_sentinel_with_nan(
            phospho,
            self.schema.phospho_cols,
            sentinel=sentinel,
        )
        return filter_min_observed(
            phospho,
            self.schema.phospho_cols,
            min_observed=min_observed,
        )


@dataclass(frozen=True, slots=True)
class ProteinCorrectionService:
    """Correct phosphosite intensities against total protein and add contrasts."""

    schema: DatasetSchema
    comparisons: tuple[ComparisonSpec, ...] | None = None

    def __init__(
        self,
        *,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> None:
        object.__setattr__(self, "schema", schema)
        object.__setattr__(
            self,
            "comparisons",
            tuple(comparisons) if comparisons is not None else None,
        )

    def correct(
        self,
        phospho_df: pd.DataFrame,
        total_df: pd.DataFrame,
        *,
        phospho_gene_col: str = "gene_names",
        total_gene_col: str = "genes",
        max_unmatched_fraction: float = 0.0,
    ) -> pd.DataFrame:
        validate_protein_correction_inputs(
            phospho_df,
            total_df,
            phospho_gene_col=phospho_gene_col,
            total_gene_col=total_gene_col,
            phospho_cols=self.schema.phospho_cols,
            protein_cols=self.schema.total_cols,
            max_unmatched_fraction=max_unmatched_fraction,
        )
        return correct_phospho_to_protein(
            df_phospho=phospho_df,
            df_total=total_df,
            phospho_gene_col=phospho_gene_col,
            total_gene_col=total_gene_col,
            phospho_cols=self.schema.phospho_cols,
            protein_cols=self.schema.total_cols,
            corrected_cols=self.schema.corrected_cols,
        )

    def add_pairwise_comparisons(
        self,
        corrected_df: pd.DataFrame,
        *,
        output_prefix: str = "p_",
    ) -> pd.DataFrame:
        if not self.comparisons:
            return corrected_df.copy()

        return add_pairwise_comparisons(
            corrected_df,
            comparisons=self.comparisons,
            group_to_corrected_col=self.schema.group_to_corrected_col,
            output_prefix=output_prefix,
        )
