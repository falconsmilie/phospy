"""Placeholder total-protein correction stage for future science restoration."""

from __future__ import annotations

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION,
    PreprocessingState,
)


class TotalProteinCorrectionStage:
    """Reserved stage boundary for total-protein correction science."""

    stage_key = DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION

    def run(self, state: PreprocessingState) -> PreprocessingState:
        return state


__all__ = ["TotalProteinCorrectionStage"]
