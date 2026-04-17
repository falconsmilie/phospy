from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from ..validation.compatibility import PredMatOverlapSummary

__all__ = ["KinaseActivityResult"]


@dataclass(slots=True)
class KinaseActivityResult:
    """Kinase activity tables produced by one analyzer run."""

    weighted_activity: pd.DataFrame
    ksea_scores: pd.DataFrame
    ksea_counts: pd.Series
    target_counts: pd.Series
    target_table: pd.DataFrame
    overlap_summary: PredMatOverlapSummary
