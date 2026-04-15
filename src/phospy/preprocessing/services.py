from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from ..datasets.schema import DatasetSchema
from ..internal.constants import (
    GENE_P_SITE_COLUMN,
    LOCALIZATION_PROB_COLUMN,
    PHOSPHO_GENE_COLUMN,
    TOTAL_GENE_COLUMN,
    ComparisonSpec,
)
from ..internal.defaults import (
    DEFAULT_LOCALIZATION_THRESHOLD,
    DEFAULT_MAX_UNMATCHED_FRACTION,
    DEFAULT_MIN_OBSERVED_VALUES,
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
)
from .primitives import (
    _add_pairwise_comparisons_in_place,
    _collapse_duplicate_genes_owned,
    _filter_localized_sites_without_copy,
    _filter_min_observed_without_copy,
    _replace_sentinel_with_nan_in_place,
)
from .protein_correction import run_protein_correction, run_protein_correction_owned

"""Concrete preprocessing transform services.

These classes do the real table-level work for the reduced preprocessing stack:

- `DatasetPreprocessing` binds a dataset workspace to the core path
- `CoreProcessor` orchestrates full and phospho-only runs
- these services perform the concrete table transforms used by that path

They remain available for advanced use, but they are not the preferred public
entrypoint for routine preprocessing.
"""

_OWNED_FRAME_ATTR = "_phospy_owned_frame"


def _mark_owned_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Mark a frame as owned by the preprocessing pipeline."""
    frame.attrs[_OWNED_FRAME_ATTR] = True
    return frame


def _is_owned_numeric_frame(
    frame: pd.DataFrame,
    *,
    required_columns: tuple[str, ...],
    numeric_columns: tuple[str, ...],
) -> bool:
    """Return whether a frame can safely take the owned no-copy correction path."""
    if not bool(frame.attrs.get(_OWNED_FRAME_ATTR, False)):
        return False
    if any(column not in frame.columns for column in required_columns):
        return False
    return all(
        pd.api.types.is_numeric_dtype(frame[column]) for column in numeric_columns
    )


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
        min_observed: int = DEFAULT_MIN_OBSERVED_VALUES,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        total = total_df.copy()
        return self.prepare_owned(
            total,
            gene_col=gene_col,
            sentinel=sentinel,
            min_observed=min_observed,
        )

    def prepare_owned(
        self,
        total_df: pd.DataFrame,
        *,
        gene_col: str = TOTAL_GENE_COLUMN,
        sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        min_observed: int = DEFAULT_MIN_OBSERVED_VALUES,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        total_df[gene_col] = total_df[gene_col].astype("string")
        _replace_sentinel_with_nan_in_place(
            total_df,
            self.schema.total_cols,
            sentinel=sentinel,
        )
        total_unique = _collapse_duplicate_genes_owned(
            total_df,
            gene_col=gene_col,
            value_cols=self.schema.total_cols,
        )
        total_filtered = _filter_min_observed_without_copy(
            total_unique,
            self.schema.total_cols,
            min_observed=min_observed,
        )
        return _mark_owned_frame(total_unique), _mark_owned_frame(total_filtered)


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
        localization_threshold: float = DEFAULT_LOCALIZATION_THRESHOLD,
        sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        min_observed: int = DEFAULT_MIN_OBSERVED_VALUES,
    ) -> pd.DataFrame:
        phospho = phospho_df.copy()
        return self.prepare_owned(
            phospho,
            gene_col=gene_col,
            site_col=site_col,
            localization_col=localization_col,
            localization_threshold=localization_threshold,
            sentinel=sentinel,
            min_observed=min_observed,
        )

    def prepare_owned(
        self,
        phospho_df: pd.DataFrame,
        *,
        gene_col: str = PHOSPHO_GENE_COLUMN,
        site_col: str = GENE_P_SITE_COLUMN,
        localization_col: str = LOCALIZATION_PROB_COLUMN,
        localization_threshold: float = DEFAULT_LOCALIZATION_THRESHOLD,
        sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        min_observed: int = DEFAULT_MIN_OBSERVED_VALUES,
    ) -> pd.DataFrame:
        phospho_df[gene_col] = phospho_df[gene_col].astype("string").str.upper()
        phospho_df[site_col] = phospho_df[site_col].astype("string")

        _replace_sentinel_with_nan_in_place(
            phospho_df,
            self.schema.phospho_cols,
            sentinel=sentinel,
        )
        phospho_filtered = _filter_localized_sites_without_copy(
            phospho_df,
            localization_col=localization_col,
            threshold=localization_threshold,
        )
        return _mark_owned_frame(
            _filter_min_observed_without_copy(
                phospho_filtered,
                self.schema.phospho_cols,
                min_observed=min_observed,
            )
        )


@dataclass(slots=True)
class ProteinCorrectionService:
    """Correct phosphosite intensities against total protein and add contrasts."""

    schema: DatasetSchema
    comparisons: tuple[ComparisonSpec, ...] | None = None

    def __post_init__(self) -> None:
        self.comparisons = self.schema.validate_comparisons(
            self.comparisons,
            context="Protein correction service",
        )

    def correct(
        self,
        phospho_df: pd.DataFrame,
        total_df: pd.DataFrame,
        *,
        phospho_gene_col: str = PHOSPHO_GENE_COLUMN,
        total_gene_col: str = TOTAL_GENE_COLUMN,
        max_unmatched_fraction: float = DEFAULT_MAX_UNMATCHED_FRACTION,
    ) -> pd.DataFrame:
        required_phospho_columns = (phospho_gene_col, *self.schema.phospho_cols)
        required_total_columns = (total_gene_col, *self.schema.total_cols)
        if _is_owned_numeric_frame(
            phospho_df,
            required_columns=required_phospho_columns,
            numeric_columns=tuple(self.schema.phospho_cols),
        ) and _is_owned_numeric_frame(
            total_df,
            required_columns=required_total_columns,
            numeric_columns=tuple(self.schema.total_cols),
        ):
            return self.correct_owned(
                phospho_df,
                total_df,
                phospho_gene_col=phospho_gene_col,
                total_gene_col=total_gene_col,
                max_unmatched_fraction=max_unmatched_fraction,
            )

        return run_protein_correction(
            df_phospho=phospho_df,
            df_total=total_df,
            phospho_gene_col=phospho_gene_col,
            total_gene_col=total_gene_col,
            phospho_cols=self.schema.phospho_cols,
            protein_cols=self.schema.total_cols,
            corrected_cols=self.schema.corrected_cols,
            max_unmatched_fraction=max_unmatched_fraction,
        )

    def correct_owned(
        self,
        phospho_df: pd.DataFrame,
        total_df: pd.DataFrame,
        *,
        phospho_gene_col: str = PHOSPHO_GENE_COLUMN,
        total_gene_col: str = TOTAL_GENE_COLUMN,
        max_unmatched_fraction: float = DEFAULT_MAX_UNMATCHED_FRACTION,
    ) -> pd.DataFrame:
        return _mark_owned_frame(
            run_protein_correction_owned(
                df_phospho=phospho_df,
                df_total=total_df,
                phospho_gene_col=phospho_gene_col,
                total_gene_col=total_gene_col,
                phospho_cols=self.schema.phospho_cols,
                protein_cols=self.schema.total_cols,
                corrected_cols=self.schema.corrected_cols,
                max_unmatched_fraction=max_unmatched_fraction,
            )
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
