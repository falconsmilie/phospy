from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

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
)
from .dataset_schema import DatasetSchema

"""Bound dataset preprocessing facade.

`DatasetPreprocessing` is the preferred public entrypoint for running the core
preprocessing path. Lower-level step services remain available in
`phospy.core_processing` and `phospy.preprocessing_services` for advanced use,
but are intentionally not mirrored as separate bound public methods here.
"""


@dataclass(frozen=True, slots=True)
class DatasetPreprocessing:
    """Bound preprocessing facade for a validated phosphoproteomics dataset.

    `run()` is the single preferred public entrypoint for dataset-bound core
    preprocessing. Advanced stepwise orchestration lives in the lower-level
    processing modules rather than being re-exposed here as overlapping bound
    methods.
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
