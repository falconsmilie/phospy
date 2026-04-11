from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from ..datasets.schema import DatasetSchema
from ..internal.constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from .core import CorePreprocessingConfig
from .modes import AnalysisReadyDatasetBuilder


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
):
    """Build the analysis-ready phosphosite boundary from user-shaped inputs."""

    return AnalysisReadyDatasetBuilder().build(
        phospho=phospho,
        total=total,
        phospho_encoding=phospho_encoding,
        schema=schema,
        comparisons=comparisons,
        preprocessing_config=preprocessing_config,
        localization_threshold=localization_threshold,
        min_observed=min_observed,
        max_unmatched_fraction=max_unmatched_fraction,
        total_sentinel=total_sentinel,
        phospho_sentinel=phospho_sentinel,
        source=source,
        phospho_only_source=phospho_only_source,
    )
