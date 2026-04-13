from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from ..datasets.schema import DatasetSchema
from ..internal.constants import ComparisonSpec
from .core import CorePreprocessingConfig
from .modes import AnalysisReadyDatasetBuilder


def build_analysis_ready_dataset(
    *,
    phospho: pd.DataFrame | str | Path,
    preprocessing_config: CorePreprocessingConfig,
    total: pd.DataFrame | str | Path | None = None,
    phospho_encoding: str | None = None,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[ComparisonSpec] | None = None,
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
        source=source,
        phospho_only_source=phospho_only_source,
    )
