"""Activity stage result models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True, slots=True)
class KinaseActivityResult:
    """Activity-stage outputs."""

    activity_scores: pd.DataFrame
