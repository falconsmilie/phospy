"""Builder-facing adapter for the internal preprocessing subsystem."""

from __future__ import annotations

import pandas as pd

from phospy.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.datasets.preprocessing.models import PreprocessingPlan, PreprocessingState
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline


class DatasetPreprocessor:
    """Translate builder input into internal preprocessing state and run it."""

    def __init__(self, *, pipeline: PreprocessingPipeline | None = None) -> None:
        self._pipeline = pipeline or PreprocessingPipeline()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        plan: PreprocessingPlan,
    ) -> PreprocessedDatasetBuildTables:
        preprocessed_state = self._pipeline.run(
            PreprocessingState(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                plan=plan,
            )
        )
        return PreprocessedDatasetBuildTables(
            phospho=preprocessed_state.phospho,
            site_metadata=preprocessed_state.site_metadata,
            sample_metadata=preprocessed_state.sample_metadata,
            total=preprocessed_state.total,
            comparisons=preprocessed_state.comparisons,
        )
