from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import pandas as pd

from ..datasets.schema import DatasetSchema
from ..internal.constants import ComparisonSpec
from .core import (
    CorePreprocessingConfig,
    CoreProcessingResult,
    CoreProcessor,
)

if TYPE_CHECKING:
    from ..datasets.models import AnalysisReadyPhosphoDataset

"""Bound dataset preprocessing facade.

`DatasetPreprocessing` is the preferred public entrypoint for running the core
preprocessing path over an owned dataset workspace. The reduced flow is:

1. resolve one preprocessing config
2. run `CoreProcessor` over the bound dataset tables
3. optionally adapt the result into `AnalysisReadyPhosphoDataset`
"""


@dataclass(slots=True)
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
    _inputs_owned: bool = field(default=False, repr=False, compare=False)

    @classmethod
    def from_owned(
        cls,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> DatasetPreprocessing:
        """Build a facade over already-owned mutable dataset tables."""
        return cls(
            total_df=total_df,
            phospho_df=phospho_df,
            schema=schema,
            comparisons=comparisons,
            _inputs_owned=True,
        )

    def __post_init__(self) -> None:
        self.comparisons = (
            tuple(self.comparisons) if self.comparisons is not None else None
        )

    def run(
        self,
        *,
        config: CorePreprocessingConfig,
    ) -> CoreProcessingResult:
        """Run preprocessing while respecting the facade's ownership boundary.

        Default construction routes through the defensive `CoreProcessor.process()`
        boundary. `from_owned()` opts into `process_owned()` for trusted
        already-owned dataset tables.
        """
        if not isinstance(config, CorePreprocessingConfig):
            msg = "DatasetPreprocessing.run(): config must be a CorePreprocessingConfig instance"
            raise TypeError(msg)
        processor = CoreProcessor(
            schema=self.schema,
            comparisons=self.comparisons,
        )
        if self._inputs_owned:
            return processor.process_owned(
                self.total_df,
                self.phospho_df,
                config=config,
            )
        return processor.process(
            self.total_df,
            self.phospho_df,
            config=config,
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
        from ..datasets.models import AnalysisReadyPhosphoDataset

        return AnalysisReadyPhosphoDataset.from_core_processing_result(
            result,
            schema=self.schema,
            comparisons=self.comparisons,
            source=source,
        )

    def run_analysis_ready(
        self,
        *,
        config: CorePreprocessingConfig,
        source: str = "dataset preprocessing",
    ) -> AnalysisReadyPhosphoDataset:
        """Run preprocessing and return the supported analysis-ready boundary."""
        core_result = self.run(config=config)
        return self.to_analysis_ready(core_result, source=source)
