from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from ..constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    TOTAL_GENE_COLUMN,
    ComparisonSpec,
)
from ..core_processing import (
    CorePreprocessingConfig,
    CoreProcessingResult,
    resolve_core_preprocessing_config,
)
from ..datasets import DatasetSiteMatrix
from ..datasets.loaders import DatasetLoader
from ..datasets.models import AnalysisReadyPhosphoDataset, PhosphoDataset
from ..datasets.schema import DatasetSchema
from ..preprocessing_services import PhosphoPreprocessor, ProteinCorrectionService


def build_analysis_ready_dataset(
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
    """Build the analysis-ready phosphosite boundary from user-shaped inputs."""

    resolved_schema = DatasetSchema() if schema is None else schema
    dataset_loader = DatasetLoader(schema=resolved_schema)

    if total is None:
        return _build_phospho_only_analysis_ready_dataset(
            phospho=phospho,
            phospho_encoding=phospho_encoding,
            schema=resolved_schema,
            comparisons=comparisons,
            preprocessing_config=preprocessing_config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            phospho_sentinel=phospho_sentinel,
            source=phospho_only_source,
            dataset_loader=dataset_loader,
        )

    phospho_df = _resolve_phospho_input(
        phospho,
        dataset_loader=dataset_loader,
        phospho_encoding=phospho_encoding,
    )
    total_df = _resolve_total_input(total, dataset_loader=dataset_loader)
    dataset = PhosphoDataset.from_loaded_inputs(
        dataset_loader.validate_inputs(total_df=total_df, phospho_df=phospho_df),
        comparisons=comparisons,
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


def _build_phospho_only_analysis_ready_dataset(
    *,
    phospho: pd.DataFrame | str | Path,
    phospho_encoding: str | None,
    schema: DatasetSchema,
    comparisons: Sequence[ComparisonSpec] | None,
    preprocessing_config: CorePreprocessingConfig | None,
    localization_threshold: float,
    min_observed: int,
    phospho_sentinel: float | int,
    source: str,
    dataset_loader: DatasetLoader,
) -> AnalysisReadyPhosphoDataset:
    phospho_df = _resolve_phospho_input(
        phospho,
        dataset_loader=dataset_loader,
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
    phospho_filtered = PhosphoPreprocessor(schema=schema).prepare(
        phospho_df,
        localization_threshold=resolved_config.localization_threshold,
        sentinel=resolved_config.phospho_sentinel,
        min_observed=resolved_config.min_observed,
    )
    phospho_corrected = phospho_filtered.rename(
        columns=dict(zip(schema.phospho_cols, schema.corrected_cols, strict=True))
    )
    phospho_corrected = ProteinCorrectionService(
        schema=schema,
        comparisons=comparisons,
    ).add_pairwise_comparisons(phospho_corrected)
    site_matrix = DatasetSiteMatrix(schema=schema).build(phospho_corrected)
    core_result = CoreProcessingResult(
        total_unique=pd.DataFrame(columns=[TOTAL_GENE_COLUMN, *schema.total_cols]),
        total_filtered=pd.DataFrame(columns=[TOTAL_GENE_COLUMN, *schema.total_cols]),
        phospho_filtered=phospho_filtered,
        phospho_corrected=phospho_corrected,
        site_matrix=site_matrix,
    )
    return AnalysisReadyPhosphoDataset.from_core_processing_result(
        core_result,
        schema=schema,
        comparisons=comparisons,
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
