"""Shared component dataclasses for kinase workflow execution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.prediction.models import KinaseScoringResult

CANDIDATE_SCORE_THRESHOLD = 0.0
CANDIDATE_MIN_INCLUSION = 1


@dataclass(frozen=True, slots=True)
class KinaseScoringRunResult:
    scoring_result: KinaseScoringResult
    downstream_score_matrix: pd.DataFrame
    downstream_score_source: str
    quantified_substrates: dict[str, list[str]]


__all__ = [
    "CANDIDATE_MIN_INCLUSION",
    "CANDIDATE_SCORE_THRESHOLD",
    "KinaseScoringRunResult",
]
