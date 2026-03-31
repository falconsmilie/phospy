from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
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
from .writers import CoreOutputWriter, CoreProcessingResultWriter


@dataclass(frozen=True, slots=True)
class CoreInputs:
    """Validated in-memory tables used by the core preprocessing pipeline."""

    total_df: pd.DataFrame
    phospho_df: pd.DataFrame


@dataclass(frozen=True, slots=True, init=False)
class PhosphoDataset:
    """Thin immutable holder around validated phosphoproteomics inputs."""

    inputs: CoreInputs
    total_cols: tuple[str, ...]
    phospho_cols: tuple[str, ...]
    corrected_cols: tuple[str, ...]
    comparisons: tuple[ComparisonSpec, ...] | None

    def __init__(
        self,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        total_cols: Sequence[str] | None = None,
        phospho_cols: Sequence[str] | None = None,
        corrected_cols: Sequence[str] | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> None:
        resolved_total_cols = self._resolve_total_cols(total_cols)
        resolved_phospho_cols = self._resolve_phospho_cols(phospho_cols)
        loader = DatasetLoader(
            total_cols=resolved_total_cols,
            phospho_cols=resolved_phospho_cols,
        )
        validated_total, validated_phospho = loader.validate(
            total_df=total_df,
            phospho_df=phospho_df,
        )
        self._set_state(
            inputs=CoreInputs(total_df=validated_total, phospho_df=validated_phospho),
            total_cols=resolved_total_cols,
            phospho_cols=resolved_phospho_cols,
            corrected_cols=corrected_cols,
            comparisons=comparisons,
        )

    def _set_state(
        self,
        *,
        inputs: CoreInputs,
        total_cols: Sequence[str],
        phospho_cols: Sequence[str],
        corrected_cols: Sequence[str] | None,
        comparisons: Sequence[ComparisonSpec] | None,
    ) -> None:
        resolved_corrected_cols = self._resolve_corrected_cols(corrected_cols)
        validate_core_column_alignment(
            total_cols,
            phospho_cols,
            resolved_corrected_cols,
        )

        object.__setattr__(self, "inputs", inputs)
        object.__setattr__(self, "total_cols", tuple(total_cols))
        object.__setattr__(self, "phospho_cols", tuple(phospho_cols))
        object.__setattr__(self, "corrected_cols", tuple(resolved_corrected_cols))
        object.__setattr__(
            self,
            "comparisons",
            tuple(comparisons) if comparisons is not None else None,
        )

    @staticmethod
    def _resolve_total_cols(total_cols: Sequence[str] | None) -> list[str]:
        return list(total_cols or DEFAULT_TOTAL_COLS)

    @staticmethod
    def _resolve_phospho_cols(phospho_cols: Sequence[str] | None) -> list[str]:
        return list(phospho_cols or DEFAULT_PHOSPHO_COLS)

    @staticmethod
    def _resolve_corrected_cols(corrected_cols: Sequence[str] | None) -> list[str]:
        return list(corrected_cols or DEFAULT_CORRECTED_COLS)

    @property
    def total_df(self) -> pd.DataFrame:
        return self.inputs.total_df

    @property
    def phospho_df(self) -> pd.DataFrame:
        return self.inputs.phospho_df

    def _core_processor(self) -> CoreProcessor:
        return CoreProcessor(
            total_cols=self.total_cols,
            phospho_cols=self.phospho_cols,
            corrected_cols=self.corrected_cols,
            comparisons=self.comparisons,
        )

    def _site_matrix_builder(self) -> SiteMatrixBuilder:
        return SiteMatrixBuilder(value_cols=self.corrected_cols)

    @classmethod
    def from_validated_inputs(
        cls,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        total_cols: Sequence[str] | None = None,
        phospho_cols: Sequence[str] | None = None,
        corrected_cols: Sequence[str] | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> PhosphoDataset:
        """Build a dataset from already-validated in-memory inputs."""
        instance = cls.__new__(cls)
        instance._set_state(
            inputs=CoreInputs(total_df=total_df, phospho_df=phospho_df),
            total_cols=cls._resolve_total_cols(total_cols),
            phospho_cols=cls._resolve_phospho_cols(phospho_cols),
            corrected_cols=corrected_cols,
            comparisons=comparisons,
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
        return cls.from_validated_inputs(
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
        return self._core_processor().process(
            self.total_df,
            self.phospho_df,
            config=resolved,
        )

    @staticmethod
    def write_core_outputs(
        result: CoreProcessingResult,
        outdir: str | Path,
        *,
        writer: CoreProcessingResultWriter = CoreOutputWriter,
    ) -> None:
        writer.write(result, outdir)
