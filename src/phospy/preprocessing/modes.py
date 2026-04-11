from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from ..datasets.builders import DatasetSiteMatrix
from ..datasets.loaders import DatasetLoader, LoadedDatasetInputs
from ..datasets.models import AnalysisReadyPhosphoDataset, PhosphoDataset
from ..datasets.schema import DatasetSchema
from ..internal.constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    TOTAL_GENE_COLUMN,
    ComparisonSpec,
)
from .core import (
    CorePreprocessingConfig,
    CoreProcessingResult,
    resolve_core_preprocessing_config,
)
from .services import PhosphoPreprocessor, ProteinCorrectionService


class FullAnalysisReadyPreprocessor:
    """Build analysis-ready datasets from full phospho and total inputs."""

    def __init__(
        self,
        *,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None = None,
        dataset_loader: DatasetLoader | None = None,
    ) -> None:
        self.schema = schema
        self.comparisons = comparisons
        self.dataset_loader = dataset_loader or DatasetLoader(schema=schema)

    def build(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        total: pd.DataFrame | str | Path,
        phospho_encoding: str | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        source: str = "analysis ready dataset builder",
    ) -> AnalysisReadyPhosphoDataset:
        loaded_inputs = _resolve_full_dataset_inputs(
            total=total,
            phospho=phospho,
            dataset_loader=self.dataset_loader,
            phospho_encoding=phospho_encoding,
        )
        dataset = PhosphoDataset.from_loaded_inputs(
            loaded_inputs,
            comparisons=self.comparisons,
        )
        return dataset.run_analysis_ready(
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            config=preprocessing_config,
            source=source,
        )


class PhosphoOnlyAnalysisReadyPreprocessor:
    """Build analysis-ready datasets from phosphosite-only inputs."""

    def __init__(
        self,
        *,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None = None,
        dataset_loader: DatasetLoader | None = None,
    ) -> None:
        self.schema = schema
        self.comparisons = comparisons
        self.dataset_loader = dataset_loader or DatasetLoader(schema=schema)

    def build(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        phospho_encoding: str | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        source: str = "analysis ready dataset builder (phospho only)",
    ) -> AnalysisReadyPhosphoDataset:
        phospho_df = _resolve_phospho_input(
            phospho,
            dataset_loader=self.dataset_loader,
            phospho_encoding=phospho_encoding,
        )
        resolved_config = resolve_core_preprocessing_config(
            config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            phospho_sentinel=phospho_sentinel,
            context="build_analysis_ready_dataset()",
            config_param_name="preprocessing_config",
        )
        phospho_filtered = PhosphoPreprocessor(schema=self.schema).prepare(
            phospho_df,
            localization_threshold=resolved_config.localization_threshold,
            sentinel=resolved_config.phospho_sentinel,
            min_observed=resolved_config.min_observed,
        )
        phospho_corrected = phospho_filtered.rename(
            columns=dict(
                zip(self.schema.phospho_cols, self.schema.corrected_cols, strict=True)
            )
        )
        phospho_corrected = ProteinCorrectionService(
            schema=self.schema,
            comparisons=self.comparisons,
        ).add_pairwise_comparisons(phospho_corrected)
        site_matrix = DatasetSiteMatrix(schema=self.schema).build(phospho_corrected)
        core_result = CoreProcessingResult(
            total_unique=pd.DataFrame(
                columns=[TOTAL_GENE_COLUMN, *self.schema.total_cols]
            ),
            total_filtered=pd.DataFrame(
                columns=[TOTAL_GENE_COLUMN, *self.schema.total_cols]
            ),
            phospho_filtered=phospho_filtered,
            phospho_corrected=phospho_corrected,
            site_matrix=site_matrix,
        )
        return AnalysisReadyPhosphoDataset.from_core_processing_result(
            core_result,
            schema=self.schema,
            comparisons=self.comparisons,
            source=source,
        )


class AnalysisReadyDatasetBuilder:
    """Stable analysis-ready preprocessing boundary used by public workflows."""

    def build(
        self,
        *,
        phospho: pd.DataFrame | str | Path,
        total: pd.DataFrame | str | Path | None = None,
        phospho_encoding: str | None = None,
        schema: DatasetSchema | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
        preprocessing_config: CorePreprocessingConfig | None = None,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        source: str = "analysis ready dataset builder",
        phospho_only_source: str = "analysis ready dataset builder (phospho only)",
    ) -> AnalysisReadyPhosphoDataset:
        resolved_schema = DatasetSchema() if schema is None else schema
        dataset_loader = DatasetLoader(schema=resolved_schema)
        if total is None:
            return PhosphoOnlyAnalysisReadyPreprocessor(
                schema=resolved_schema,
                comparisons=comparisons,
                dataset_loader=dataset_loader,
            ).build(
                phospho=phospho,
                phospho_encoding=phospho_encoding,
                preprocessing_config=preprocessing_config,
                localization_threshold=localization_threshold,
                min_observed=min_observed,
                phospho_sentinel=phospho_sentinel,
                source=phospho_only_source,
            )
        return FullAnalysisReadyPreprocessor(
            schema=resolved_schema,
            comparisons=comparisons,
            dataset_loader=dataset_loader,
        ).build(
            phospho=phospho,
            total=total,
            phospho_encoding=phospho_encoding,
            preprocessing_config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            source=source,
        )


def _resolve_total_input(
    total: pd.DataFrame | str | Path,
    *,
    dataset_loader: DatasetLoader,
) -> pd.DataFrame:
    if isinstance(total, pd.DataFrame):
        return dataset_loader.validate_total(total)
    return dataset_loader.load_total(total)


def _resolve_phospho_input(
    phospho: pd.DataFrame | str | Path,
    *,
    dataset_loader: DatasetLoader,
    phospho_encoding: str | None,
) -> pd.DataFrame:
    if isinstance(phospho, pd.DataFrame):
        return dataset_loader.validate_phospho(phospho)
    return dataset_loader.load_phospho(phospho, encoding=phospho_encoding)


def _resolve_full_dataset_inputs(
    *,
    total: pd.DataFrame | str | Path,
    phospho: pd.DataFrame | str | Path,
    dataset_loader: DatasetLoader,
    phospho_encoding: str | None,
) -> LoadedDatasetInputs:
    if isinstance(total, pd.DataFrame) and isinstance(phospho, pd.DataFrame):
        return dataset_loader.validate_inputs(total_df=total, phospho_df=phospho)
    if not isinstance(total, pd.DataFrame) and not isinstance(phospho, pd.DataFrame):
        return dataset_loader.load(
            total,
            phospho,
            phospho_encoding=phospho_encoding,
        )

    total_df = _resolve_total_input(total, dataset_loader=dataset_loader)
    phospho_df = _resolve_phospho_input(
        phospho,
        dataset_loader=dataset_loader,
        phospho_encoding=phospho_encoding,
    )
    return LoadedDatasetInputs(
        total_df=total_df,
        phospho_df=phospho_df,
        schema=dataset_loader.schema,
    )
