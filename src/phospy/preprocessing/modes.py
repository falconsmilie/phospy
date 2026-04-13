from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..datasets.loaders import DatasetLoader
from ..datasets.models import AnalysisReadyPhosphoDataset
from ..datasets.schema import DatasetSchema
from ..errors import InputCompatibilityError
from ..internal.constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from .core import (
    CorePreprocessingConfig,
    CoreProcessor,
    resolve_core_preprocessing_config,
)
from .dataset import DatasetPreprocessing

"""Analysis-ready preprocessing builder.

This module intentionally keeps one high-level builder instead of separate
full-input and phospho-only orchestration wrappers. The reduced flow is:

1. resolve validated dataset or phospho inputs through `DatasetLoader`
2. reuse `DatasetPreprocessing` or `CoreProcessor.process_phospho_only()`
3. adapt the core result into `AnalysisReadyPhosphoDataset`
"""


@dataclass(slots=True)
class AnalysisReadyDatasetBuilder:
    """Stable analysis-ready preprocessing boundary used by public workflows."""

    dataset_loader: DatasetLoader | None = None

    def _resolve_schema(self, schema: DatasetSchema | None) -> DatasetSchema:
        return DatasetSchema() if schema is None else schema

    def _resolve_loader(self, schema: DatasetSchema) -> DatasetLoader:
        if self.dataset_loader is None:
            return DatasetLoader(schema=schema)
        if self.dataset_loader.schema != schema:
            msg = (
                "AnalysisReadyDatasetBuilder dataset loader schema must match the "
                "requested preprocessing schema"
            )
            raise InputCompatibilityError(msg)
        return self.dataset_loader

    def _build_from_full_inputs(
        self,
        *,
        loader: DatasetLoader,
        schema: DatasetSchema,
        total: pd.DataFrame | str | Path,
        phospho: pd.DataFrame | str | Path,
        phospho_encoding: str | None,
        comparisons: Sequence[ComparisonSpec] | None,
        preprocessing_config: CorePreprocessingConfig | None,
        localization_threshold: float,
        min_observed: int,
        max_unmatched_fraction: float,
        total_sentinel: float | int,
        phospho_sentinel: float | int,
        source: str,
    ) -> AnalysisReadyPhosphoDataset:
        loaded_inputs = loader.resolve_inputs(
            total=total,
            phospho=phospho,
            phospho_encoding=phospho_encoding,
        )
        preprocessing = DatasetPreprocessing(
            total_df=loaded_inputs.total_df,
            phospho_df=loaded_inputs.phospho_df,
            schema=schema,
            comparisons=comparisons,
        )
        return preprocessing.run_analysis_ready(
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            config=preprocessing_config,
            source=source,
        )

    def _build_from_phospho_only(
        self,
        *,
        loader: DatasetLoader,
        schema: DatasetSchema,
        phospho: pd.DataFrame | str | Path,
        phospho_encoding: str | None,
        comparisons: Sequence[ComparisonSpec] | None,
        preprocessing_config: CorePreprocessingConfig | None,
        localization_threshold: float,
        min_observed: int,
        phospho_sentinel: float | int,
        source: str,
    ) -> AnalysisReadyPhosphoDataset:
        phospho_df = loader.resolve_phospho(
            phospho,
            encoding=phospho_encoding,
        )
        resolved_config = resolve_core_preprocessing_config(
            config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            phospho_sentinel=phospho_sentinel,
            context="build_analysis_ready_dataset()",
            config_param_name="preprocessing_config",
        )
        core_result = CoreProcessor(
            schema=schema,
            comparisons=comparisons,
        ).process_phospho_only(
            phospho_df,
            config=resolved_config,
        )
        return AnalysisReadyPhosphoDataset.from_core_processing_result(
            core_result,
            schema=schema,
            comparisons=comparisons,
            source=source,
        )

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
        resolved_schema = self._resolve_schema(schema)
        loader = self._resolve_loader(resolved_schema)
        if total is None:
            return self._build_from_phospho_only(
                loader=loader,
                schema=resolved_schema,
                phospho=phospho,
                phospho_encoding=phospho_encoding,
                comparisons=comparisons,
                preprocessing_config=preprocessing_config,
                localization_threshold=localization_threshold,
                min_observed=min_observed,
                phospho_sentinel=phospho_sentinel,
                source=phospho_only_source,
            )
        return self._build_from_full_inputs(
            loader=loader,
            schema=resolved_schema,
            total=total,
            phospho=phospho,
            phospho_encoding=phospho_encoding,
            comparisons=comparisons,
            preprocessing_config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            source=source,
        )
