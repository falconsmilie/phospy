"""Differential-domain defensive internal views for workflow collaborators."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from phospy.science.differential.models.fit import DifferentialComputationResult
from phospy.science.differential.models.protein_aware import (
    ProteinAwareDifferentialComputationResult,
)


class DifferentialComputationResultInternalView:
    """Narrow access to mutation-isolated computation result tables."""

    __slots__ = ("_result",)

    def __init__(self, result: DifferentialComputationResult) -> None:
        self._result = result

    @property
    def contrast_tables(self) -> Mapping[str, pd.DataFrame]:
        return self._result._borrow_contrast_tables()


class ProteinAwareDifferentialComputationResultInternalView:
    """Narrow access to mutation-isolated protein-aware computation tables."""

    __slots__ = ("_result",)

    def __init__(self, result: ProteinAwareDifferentialComputationResult) -> None:
        self._result = result

    @property
    def contrast_tables(self) -> Mapping[str, pd.DataFrame]:
        return self._result.contrast_tables


__all__ = [
    "DifferentialComputationResultInternalView",
    "ProteinAwareDifferentialComputationResultInternalView",
]
