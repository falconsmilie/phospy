"""Activity stage result models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True, slots=True)
class KinaseActivityResult:
    """Activity-stage outputs."""

    activity_scores: pd.DataFrame

    def __post_init__(self) -> None:
        if not isinstance(self.activity_scores, pd.DataFrame):
            raise TypeError(
                "activity_result.activity_scores must be a pandas DataFrame"
            )
        object.__setattr__(
            self, "activity_scores", self.activity_scores.copy(deep=True)
        )
