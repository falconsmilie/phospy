from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import pandas as pd

from .constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from .core_processing import (
    CorePreprocessingConfig,
    CoreProcessingResult,
    CoreProcessor,
    resolve_core_preprocessing_config,
)
from .datasets.schema import DatasetSchema

if TYPE_CHECKING:
    from .datasets.models import AnalysisReadyPhosphoDataset

"""Bound dataset preprocessing facade.

`DatasetPreprocessing` is the preferred public entrypoint for running the core
preprocessing path. Lower-level step services remain available in
`phospy.core_processing` and `phospy.preprocessing_services` for advanced use,
but are intentionally not mirrored as separate bound public methods here.
"""


@dataclass(frozen=True, slots=True)
class DatasetPreprocessing:
    """Bound preprocessing facade for a validated phosphoproteomics dataset.

    `run()` and `run_analysis_ready()` are the preferred public entrypoints
    for dataset-bound preprocessing. Advanced stepwise orchestration lives in
    the lower-level processing modules.

    When created from :class:`phospy.PhosphoDataset`, this facade is bound to
    the dataset's explicit `total_df_live` and `phospho_df_live` accessors so
    in-memory processing works against the owned workspace state intentionally.
    """

    total_df: pd.DataFrame
    phospho_df: pd.DataFrame
    schema: DatasetSchema
    comparisons: tuple[ComparisonSpec, ...] | None = None

    def __init__(
        self,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> None:
        object.__setattr__(self, "total_df", total_df)
        object.__setattr__(self, "phospho_df", phospho_df)
        object.__setattr__(self, "schema", schema)
        object.__setattr__(
            self,
            "comparisons",
            tuple(comparisons) if comparisons is not None else None,
        )

    def _core_processor(self) -> CoreProcessor:
        return CoreProcessor(
            schema=self.schema,
            comparisons=self.comparisons,
        )

    def run(
        self,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        config: CorePreprocessingConfig | None = None,
    ) -> CoreProcessingResult:
        resolved = resolve_core_preprocessing_config(
            config=config,
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            max_unmatched_fraction=max_unmatched_fraction,
            context="DatasetPreprocessing.run()",
            config_param_name="config",
        )
        return self._core_processor().process(
            self.total_df,
            self.phospho_df,
            config=resolved,
        )

    def to_analysis_ready(
        self,
        result: CoreProcessingResult,
        *,
        source: str = "dataset preprocessing",
    ) -> AnalysisReadyPhosphoDataset:
        """Adapt one preprocessing result into an analysis-ready dataset.

        This is the supported bound adapter from the dataset preprocessing lane
        into :class:`phospy.AnalysisReadyPhosphoDataset`. It reuses the current
        core preprocessing and site-matrix behaviour and binds provenance to the
        schema and comparisons already owned by this preprocessing facade.
        """
        from .datasets.models import AnalysisReadyPhosphoDataset

        return AnalysisReadyPhosphoDataset.from_core_processing_result(
            result,
            schema=self.schema,
            comparisons=self.comparisons,
            source=source,
        )

    def run_analysis_ready(
        self,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        config: CorePreprocessingConfig | None = None,
        source: str = "dataset preprocessing",
    ) -> AnalysisReadyPhosphoDataset:
        """Run preprocessing and return the supported analysis-ready boundary."""
        core_result = self.run(
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            config=config,
        )
        return self.to_analysis_ready(core_result, source=source)
