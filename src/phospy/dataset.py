from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from .constants import (
    DEFAULT_CORRECTED_COLS,
    DEFAULT_PHOSPHO_COLS,
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_COLS,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from .core_processing import (
    CorePreprocessingConfig,
    CoreProcessingResult,
    CoreProcessor,
)
from .dataset_loader import DatasetLoader
from .site_matrix_builder import SiteMatrixBuilder, SiteMatrixResult
from .validation.compatibility import validate_core_column_alignment


class PhosphoDataset:
    """Facade around validated phosphoproteomics inputs and core processing."""

    def __init__(
        self,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        total_cols: Sequence[str] | None = None,
        phospho_cols: Sequence[str] | None = None,
        corrected_cols: Sequence[str] | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> None:
        resolved_total_cols = list(total_cols or DEFAULT_TOTAL_COLS)
        resolved_phospho_cols = list(phospho_cols or DEFAULT_PHOSPHO_COLS)
        loader = DatasetLoader(
            total_cols=resolved_total_cols,
            phospho_cols=resolved_phospho_cols,
        )
        validated_total, validated_phospho = loader.validate(
            total_df=total_df,
            phospho_df=phospho_df,
        )
        self._initialize(
            total_df=validated_total,
            phospho_df=validated_phospho,
            total_cols=resolved_total_cols,
            phospho_cols=resolved_phospho_cols,
            corrected_cols=corrected_cols,
            comparisons=comparisons,
            loader=loader,
        )

    def _initialize(
        self,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        total_cols: Sequence[str],
        phospho_cols: Sequence[str],
        corrected_cols: Sequence[str] | None,
        comparisons: Sequence[ComparisonSpec] | None,
        loader: DatasetLoader,
    ) -> None:
        self.total_cols = list(total_cols)
        self.phospho_cols = list(phospho_cols)
        self.corrected_cols = list(corrected_cols or DEFAULT_CORRECTED_COLS)
        validate_core_column_alignment(
            self.total_cols,
            self.phospho_cols,
            self.corrected_cols,
        )

        self.loader = loader
        self.total_df = total_df
        self.phospho_df = phospho_df
        self.comparisons = list(comparisons) if comparisons is not None else None
        self.site_matrix_builder = SiteMatrixBuilder(value_cols=self.corrected_cols)
        self.core_processor = CoreProcessor(
            total_cols=self.total_cols,
            phospho_cols=self.phospho_cols,
            corrected_cols=self.corrected_cols,
            comparisons=self.comparisons,
            site_matrix_builder=self.site_matrix_builder,
        )

    @classmethod
    def _from_validated_frames(
        cls,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        total_cols: Sequence[str] | None = None,
        phospho_cols: Sequence[str] | None = None,
        corrected_cols: Sequence[str] | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> PhosphoDataset:
        instance = cls.__new__(cls)
        resolved_total_cols = list(total_cols or DEFAULT_TOTAL_COLS)
        resolved_phospho_cols = list(phospho_cols or DEFAULT_PHOSPHO_COLS)
        instance._initialize(
            total_df=total_df,
            phospho_df=phospho_df,
            total_cols=resolved_total_cols,
            phospho_cols=resolved_phospho_cols,
            corrected_cols=corrected_cols,
            comparisons=comparisons,
            loader=DatasetLoader(
                total_cols=resolved_total_cols,
                phospho_cols=resolved_phospho_cols,
            ),
        )
        return instance

    @classmethod
    def from_files(
        cls,
        total_path: str | Path,
        phospho_path: str | Path,
        phospho_encoding: str | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> PhosphoDataset:
        loader = DatasetLoader()
        total_df, phospho_df = loader.load(
            total_path,
            phospho_path,
            phospho_encoding=phospho_encoding,
        )
        return cls._from_validated_frames(
            total_df=total_df,
            phospho_df=phospho_df,
            comparisons=comparisons,
        )

    def prepare_total(
        self,
        gene_col: str = "genes",
        sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        min_observed: int = 4,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        return self.core_processor.prepare_total(
            self.total_df,
            gene_col=gene_col,
            sentinel=sentinel,
            min_observed=min_observed,
        )

    def prepare_phospho(
        self,
        gene_col: str = "gene_names",
        site_col: str = "gene_p_site",
        uid_col: str = "uid",
        localization_col: str = "localization_prob",
        localization_threshold: float = 0.75,
        sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        min_observed: int = 4,
    ) -> pd.DataFrame:
        del uid_col
        return self.core_processor.prepare_phospho(
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
        return self.core_processor.correct_to_protein(
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
        return self.core_processor.add_pairwise_comparisons(
            corrected_df,
            output_prefix=output_prefix,
        )

    def build_site_matrix(
        self,
        corrected_df: pd.DataFrame,
        gene_p_site_col: str = "gene_p_site",
        sequence_col: str = "centralized_sequence",
    ) -> SiteMatrixResult:
        return self.site_matrix_builder.build(
            corrected_df,
            gene_p_site_col=gene_p_site_col,
            sequence_col=sequence_col,
        )

    def process_core(
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
        return self.core_processor.process(
            self.total_df,
            self.phospho_df,
            config=resolved,
        )

    @staticmethod
    def write_core_outputs(result: CoreProcessingResult, outdir: str | Path) -> None:
        from .writers import CoreOutputWriter

        CoreOutputWriter.write(result, outdir)
