"""Quantification-depth validation for differential empirical-Bayes trends."""

from __future__ import annotations

from typing import Final

import numpy as np
import pandas as pd

from phospy.evidence.quantification_depth import validate_quantification_depth_series

QUANTIFICATION_DEPTH_LOG2_TREND_COVARIATE_NAME: Final[str] = "log2_quantification_depth"
QUANTIFICATION_DEPTH_TREND_TRANSFORMATION: Final[str] = "log2"


def log2_quantification_depth_series(
    depth: pd.Series,
    *,
    field_name: str,
    expected_index: pd.Index | None = None,
) -> pd.Series:
    """Return ``log2(quantification_depth)`` after enforcing count validity."""

    validated = validate_quantification_depth_series(
        depth,
        field_name=field_name,
        expected_index=expected_index,
    )
    return pd.Series(
        np.log2(validated.to_numpy(dtype=float, copy=False)),
        index=validated.index.copy(),
        name=QUANTIFICATION_DEPTH_LOG2_TREND_COVARIATE_NAME,
    )


__all__ = [
    "QUANTIFICATION_DEPTH_LOG2_TREND_COVARIATE_NAME",
    "QUANTIFICATION_DEPTH_TREND_TRANSFORMATION",
    "log2_quantification_depth_series",
]
