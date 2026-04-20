"""Placeholder site-matrix stage for future science restoration."""

from __future__ import annotations

from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_SITE_MATRIX,
    PreprocessingState,
)


class SiteMatrixStage:
    """Reserved stage boundary for site-matrix construction science."""

    stage_key = DATASET_PREPROCESSING_STAGE_SITE_MATRIX

    def run(self, state: PreprocessingState) -> PreprocessingState:
        return state


__all__ = ["SiteMatrixStage"]
