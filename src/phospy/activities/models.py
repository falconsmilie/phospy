"""Activity stage result models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True, slots=True)
class KinaseActivityResult:
    """Activity-stage outputs.

    Expected columns include:
    - ``activity_score``: mean positive-support kinase prediction score
    - ``weighted_signal``: phosphosite signal weighted by positive-support scores
    - ``n_predicted_sites``: number of positive-support substrate sites
    - ``is_active``: threshold flag from ``activity_config.threshold``
    """

    activity_scores: pd.DataFrame

    def __post_init__(self) -> None:
        if not isinstance(self.activity_scores, pd.DataFrame):
            raise TypeError(
                "activity_result.activity_scores must be a pandas DataFrame"
            )
        object.__setattr__(
            self, "activity_scores", self.activity_scores.copy(deep=True)
        )
