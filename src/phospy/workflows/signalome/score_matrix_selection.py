"""Downstream score matrix selection for signalome interpretation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.contracts.results import KinaseScoringResult
from phospy.science.prediction.scoring import (
    SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY,
    DownstreamScoreSelectionPolicy,
    select_downstream_score_matrix,
)
from phospy.science.scoring.policy_models import DownstreamScoreSource


@dataclass(frozen=True, slots=True)
class SignalomeScoreMatrixSelection:
    downstream_score_matrix: pd.DataFrame
    downstream_score_source: DownstreamScoreSource
    downstream_score_selection_policy: DownstreamScoreSelectionPolicy


class _SelectDownstreamScoreMatrixContract(Protocol):
    def __call__(
        self,
        *,
        profile_scores: pd.DataFrame,
        rank_weighted_fusion_scores: pd.DataFrame | None,
    ) -> tuple[pd.DataFrame, DownstreamScoreSource]: ...


class SignalomeScoreMatrixSelector:
    """Select the downstream score matrix lane for signalome execution."""

    def __init__(
        self,
        *,
        select_matrix: _SelectDownstreamScoreMatrixContract = (
            select_downstream_score_matrix
        ),
        selection_policy: DownstreamScoreSelectionPolicy = (
            SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY
        ),
    ) -> None:
        self._select_matrix = select_matrix
        self._selection_policy = selection_policy

    def run(self, scoring_result: KinaseScoringResult) -> SignalomeScoreMatrixSelection:
        downstream_score_matrix, downstream_score_source = self._select_matrix(
            profile_scores=scoring_result._borrow_profile_scores_frame(),
            rank_weighted_fusion_scores=(
                scoring_result._borrow_rank_weighted_fusion_scores_frame()
            ),
        )
        return SignalomeScoreMatrixSelection(
            downstream_score_matrix=downstream_score_matrix,
            downstream_score_source=downstream_score_source,
            downstream_score_selection_policy=self._selection_policy,
        )


__all__ = [
    "SignalomeScoreMatrixSelection",
    "SignalomeScoreMatrixSelector",
]
