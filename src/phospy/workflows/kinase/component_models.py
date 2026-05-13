"""Shared component dataclasses for kinase workflow execution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.prediction.models import KinaseScoringResult
from phospy.prediction.scoring import DownstreamScoreSelectionPolicy
from phospy.scoring.policy_models import DownstreamScoreSource

CANDIDATE_SCORE_THRESHOLD = 0.0
CANDIDATE_MIN_INCLUSION = 1


@dataclass(frozen=True, slots=True)
class KinaseScoringRunResult:
    scoring_result: KinaseScoringResult
    downstream_score_matrix: pd.DataFrame
    downstream_score_source: DownstreamScoreSource
    quantified_substrates: dict[str, list[str]]
    downstream_score_selection_policy: DownstreamScoreSelectionPolicy | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "downstream_score_source",
            DownstreamScoreSource.parse(
                self.downstream_score_source,
                field_name="kinase scoring run downstream_score_source",
            ),
        )


__all__ = [
    "CANDIDATE_MIN_INCLUSION",
    "CANDIDATE_SCORE_THRESHOLD",
    "KinaseScoringRunResult",
]
