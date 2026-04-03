from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from ._preprocessing_primitives import (
    _add_pairwise_comparisons_in_place,
    _collapse_duplicate_genes_owned,
    _filter_localized_sites_without_copy,
    _filter_min_observed_without_copy,
    _replace_sentinel_with_nan_in_place,
)
from .constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    GENE_P_SITE_COLUMN,
    LOCALIZATION_PROB_COLUMN,
    PHOSPHO_GENE_COLUMN,
    TOTAL_GENE_COLUMN,
    ComparisonSpec,
)
from .dataset_schema import DatasetSchema
from .preprocessing import correct_phospho_to_protein

"""Internal preprocessing service layer.

These classes back the preferred dataset-bound preprocessing path and the core
processor. They remain available for advanced use, but they are not the
preferred public entrypoint for routine preprocessing.
"""


@dataclass(frozen=True, slots=True)
class TotalPreprocessor:
    """Prepare total proteome rows for downstream phosphosite correction."""

    schema: DatasetSchema

    def prepare(
        self,
        total_df: pd.DataFrame,
        *,
        gene_col: str = TOTAL_GENE_COLUMN,
        sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        min_observed: int = 4,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        total = total_df.copy()
        total[gene_col] = total[gene_col].astype("string")
        _replace_sentinel_with_nan_in_place(
            total,
            self.schema.total_cols,
            sentinel=sentinel,
        )
        total_unique = _collapse_duplicate_genes_owned(
            total,
            gene_col=gene_col,
            value_cols=self.schema.total_cols,
        )
        total_filtered = _filter_min_observed_without_copy(
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
        gene_col: str = PHOSPHO_GENE_COLUMN,
        site_col: str = GENE_P_SITE_COLUMN,
        localization_col: str = LOCALIZATION_PROB_COLUMN,
        localization_threshold: float = 0.75,
        sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        min_observed: int = 4,
    ) -> pd.DataFrame:
        phospho = phospho_df.copy()
        phospho[gene_col] = phospho[gene_col].astype("string").str.upper()
        phospho[site_col] = phospho[site_col].astype("string")

        phospho = _filter_localized_sites_without_copy(
            phospho,
            localization_col=localization_col,
            threshold=localization_threshold,
        )
        _replace_sentinel_with_nan_in_place(
            phospho,
            self.schema.phospho_cols,
            sentinel=sentinel,
        )
        return _filter_min_observed_without_copy(
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
            schema.validate_comparisons(
                comparisons,
                context="Protein correction service",
            ),
        )

    def correct(
        self,
        phospho_df: pd.DataFrame,
        total_df: pd.DataFrame,
        *,
        phospho_gene_col: str = PHOSPHO_GENE_COLUMN,
        total_gene_col: str = TOTAL_GENE_COLUMN,
        max_unmatched_fraction: float = 0.0,
    ) -> pd.DataFrame:
        return correct_phospho_to_protein(
            df_phospho=phospho_df,
            df_total=total_df,
            phospho_gene_col=phospho_gene_col,
            total_gene_col=total_gene_col,
            phospho_cols=self.schema.phospho_cols,
            protein_cols=self.schema.total_cols,
            corrected_cols=self.schema.corrected_cols,
            max_unmatched_fraction=max_unmatched_fraction,
        )

    def add_pairwise_comparisons(
        self,
        corrected_df: pd.DataFrame,
        *,
        output_prefix: str = "p_",
    ) -> pd.DataFrame:
        if not self.comparisons:
            return corrected_df

        return _add_pairwise_comparisons_in_place(
            corrected_df,
            comparisons=self.comparisons,
            group_to_corrected_col=self.schema.group_to_corrected_col,
            output_prefix=output_prefix,
            schema=self.schema,
        )
