"""Placeholder comparison-building stage for future science restoration."""

from __future__ import annotations

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    PreprocessingState,
)


class ComparisonsStage:
    """Reserved stage boundary for comparison-building science."""

    stage_key = DATASET_PREPROCESSING_STAGE_COMPARISONS

    def run(self, state: PreprocessingState) -> PreprocessingState:
        return state


__all__ = ["ComparisonsStage"]
