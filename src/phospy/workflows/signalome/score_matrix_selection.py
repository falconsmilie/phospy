"""Downstream score matrix selection for signalome interpretation."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import pandas as pd

from phospy.contracts.results import KinaseScoringResult
from phospy.science.prediction.scoring import (
    SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY,
    DownstreamScoreSelectionPolicy,
)
from phospy.science.scoring.policy_models import DownstreamScoreSource


@dataclass(frozen=True, slots=True)
class SignalomeScoreMatrixSelection:
    downstream_score_matrix: pd.DataFrame
    downstream_score_source: DownstreamScoreSource
    downstream_score_selection_policy: DownstreamScoreSelectionPolicy


class SignalomeScoreMatrixSelector:
    """Select the downstream score matrix lane for signalome execution."""

    def __init__(
        self,
        *,
        select_matrix: Callable[..., tuple[pd.DataFrame, object]] | None = None,
        selection_policy: DownstreamScoreSelectionPolicy = (
            SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY
        ),
    ) -> None:
        self._select_matrix = select_matrix
        self._selection_policy = selection_policy

    def run(self, scoring_result: KinaseScoringResult) -> SignalomeScoreMatrixSelection:
        if self._select_matrix is not None:
            downstream_score_matrix, raw_source = self._select_matrix(
                profile_scores=scoring_result._borrow_profile_scores_frame(),
                rank_weighted_fusion_scores=(
                    scoring_result._borrow_rank_weighted_fusion_scores_frame()
                ),
            )
            source_value = (
                raw_source.value
                if isinstance(raw_source, DownstreamScoreSource)
                else str(raw_source)
            )
            downstream_score_source = DownstreamScoreSource.parse(
                source_value,
                field_name="signalome downstream score source",
            )
            return SignalomeScoreMatrixSelection(
                downstream_score_matrix=downstream_score_matrix,
                downstream_score_source=downstream_score_source,
                downstream_score_selection_policy=self._selection_policy,
            )
        downstream_score_source = DownstreamScoreSource.parse(
            scoring_result.score_source,
            field_name="signalome downstream score source",
        )
        downstream_score_matrix = scoring_result._borrow_authoritative_scores_frame()
        return SignalomeScoreMatrixSelection(
            downstream_score_matrix=downstream_score_matrix,
            downstream_score_source=downstream_score_source,
            downstream_score_selection_policy=self._selection_policy,
        )


__all__ = [
    "SignalomeScoreMatrixSelection",
    "SignalomeScoreMatrixSelector",
]
