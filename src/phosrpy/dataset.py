from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pandas as pd

from .constants import DEFAULT_COMPARISONS, DEFAULT_CORRECTED_COLS, DEFAULT_PHOSPHO_COLS, DEFAULT_TOTAL_COLS
from .io import read_table
from .matrices import build_site_matrix
from .preprocessing import (
    add_pairwise_comparisons,
    collapse_duplicate_genes,
    correct_phospho_to_protein,
    filter_min_observed,
    replace_sentinel_with_nan,
)


@dataclass(slots=True)
class SiteMatrixResult:
    phosr_input: pd.DataFrame
    matrix: pd.DataFrame
    sequences: pd.Series


@dataclass(slots=True)
class CoreProcessingResult:
    total_unique: pd.DataFrame
    total_filtered: pd.DataFrame
    phospho_filtered: pd.DataFrame
    phospho_corrected: pd.DataFrame
    site_matrix: SiteMatrixResult


class PhosphoDataset:
    """Owns phosphoproteomics inputs and preprocessing steps.

    The class provides a domain-oriented API around the core data wrangling
    that feeds downstream PhosR-style analysis.
    """

    def __init__(
        self,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        total_cols: Sequence[str] | None = None,
        phospho_cols: Sequence[str] | None = None,
        corrected_cols: Sequence[str] | None = None,
        comparisons: Sequence[tuple[str, str, str, str]] | None = None,
    ) -> None:
        self.total_df = total_df.copy()
        self.phospho_df = phospho_df.copy()
        self.total_cols = list(total_cols or DEFAULT_TOTAL_COLS)
        self.phospho_cols = list(phospho_cols or DEFAULT_PHOSPHO_COLS)
        self.corrected_cols = list(corrected_cols or DEFAULT_CORRECTED_COLS)
        self.comparisons = list(comparisons or DEFAULT_COMPARISONS)

    @classmethod
    def from_files(
        cls,
        total_path: str | Path,
        phospho_path: str | Path,
        phospho_encoding: str | None = None,
    ) -> "PhosphoDataset":
        total_df = read_table(total_path)
        phospho_df = read_table(phospho_path, encoding=phospho_encoding)
        return cls(total_df=total_df, phospho_df=phospho_df)

    def prepare_total(
        self,
        gene_col: str = "genes",
        sentinel: float | int = 10,
        min_observed: int = 4,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        total = self.total_df.copy()
        total[gene_col] = total[gene_col].astype("string")
        for col in self.total_cols:
            total[col] = pd.to_numeric(total[col], errors="coerce")
        total = replace_sentinel_with_nan(total, self.total_cols, sentinel=sentinel)
        total_unique = collapse_duplicate_genes(total, gene_col=gene_col, value_cols=self.total_cols)
        total_filtered = filter_min_observed(total_unique, self.total_cols, min_observed=min_observed)
        return total_unique, total_filtered

    def prepare_phospho(
        self,
        gene_col: str = "gene_names",
        site_col: str = "gene_p_site",
        uid_col: str = "uid",
        localization_col: str = "localization_prob",
        localization_threshold: float = 0.75,
        sentinel: float | int = 12,
        min_observed: int = 4,
    ) -> pd.DataFrame:
        phospho = self.phospho_df.copy()
        phospho[gene_col] = phospho[gene_col].astype("string").str.upper()
        phospho[site_col] = phospho[site_col].astype("string")
        phospho[localization_col] = pd.to_numeric(phospho[localization_col], errors="coerce")
        for col in self.phospho_cols:
            phospho[col] = pd.to_numeric(phospho[col], errors="coerce")

        phospho = phospho.loc[phospho[uid_col].notna() & phospho[gene_col].notna()].copy()
        phospho = phospho.loc[phospho[localization_col] >= localization_threshold].copy()
        phospho = replace_sentinel_with_nan(phospho, self.phospho_cols, sentinel=sentinel)
        return filter_min_observed(phospho, self.phospho_cols, min_observed=min_observed)

    def correct_to_protein(
        self,
        phospho_df: pd.DataFrame,
        total_df: pd.DataFrame,
        phospho_gene_col: str = "gene_names",
        total_gene_col: str = "genes",
        output_prefix: str = "phospho_corrected_",
    ) -> pd.DataFrame:
        return correct_phospho_to_protein(
            df_phospho=phospho_df,
            df_total=total_df,
            phospho_gene_col=phospho_gene_col,
            total_gene_col=total_gene_col,
            phospho_cols=self.phospho_cols,
            protein_cols=self.total_cols,
            output_prefix=output_prefix,
        )

    def add_pairwise_comparisons(
        self,
        corrected_df: pd.DataFrame,
        output_prefix: str = "p_",
    ) -> pd.DataFrame:
        group_map = {f"group{i}": self.corrected_cols[i - 1] for i in range(1, len(self.corrected_cols) + 1)}
        return add_pairwise_comparisons(
            corrected_df,
            comparisons=self.comparisons,
            group_to_corrected_col=group_map,
            output_prefix=output_prefix,
        )

    def build_site_matrix(
        self,
        corrected_df: pd.DataFrame,
        gene_p_site_col: str = "gene_p_site",
        sequence_col: str = "centralized_sequence",
    ) -> SiteMatrixResult:
        phosr_input, matrix, sequences = build_site_matrix(
            df=corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
            value_cols=self.corrected_cols,
        )
        return SiteMatrixResult(phosr_input=phosr_input, matrix=matrix, sequences=sequences)

    def process_core(self) -> CoreProcessingResult:
        total_unique, total_filtered = self.prepare_total()
        phospho_filtered = self.prepare_phospho()
        phospho_corrected = self.correct_to_protein(phospho_filtered, total_filtered)
        phospho_corrected = self.add_pairwise_comparisons(phospho_corrected)
        site_matrix = self.build_site_matrix(phospho_corrected)
        return CoreProcessingResult(
            total_unique=total_unique,
            total_filtered=total_filtered,
            phospho_filtered=phospho_filtered,
            phospho_corrected=phospho_corrected,
            site_matrix=site_matrix,
        )

    @staticmethod
    def write_core_outputs(result: CoreProcessingResult, outdir: str | Path) -> None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        result.total_unique.to_csv(outdir / "df_total_unique.csv", index=False)
        result.total_filtered.to_csv(outdir / "df_total_filtered.csv", index=False)
        result.phospho_filtered.to_csv(outdir / "df_phospho_filtered.csv", index=False)
        result.phospho_corrected.to_csv(outdir / "df_phospho_corrected.csv", index=False)
        result.site_matrix.phosr_input.to_csv(outdir / "phosr_input.csv", index=False)
        result.site_matrix.matrix.to_csv(outdir / "mat_phospho_corrected.csv")
        result.site_matrix.sequences.rename("centralized_sequence").to_csv(outdir / "site_sequences.csv")
