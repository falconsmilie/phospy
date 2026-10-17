from __future__ import annotations

import numpy as np
import pandas as pd


def constant_column_mask(mat: pd.DataFrame) -> pd.Series:
    """Return True for columns with no measurable motif-score contrast."""

    numeric = mat.astype(float)
    return (numeric.max(axis=0) - numeric.min(axis=0)) == 0.0


def minmax_scale_columns(mat: pd.DataFrame) -> pd.DataFrame:
    """Apply column-wise min-max scaling used by baseline motif scoring."""

    scaled = mat.astype(float).copy()
    constant_mask = constant_column_mask(scaled)
    for column in scaled.columns:
        values = scaled.loc[:, column]
        min_value = float(values.min())
        max_value = float(values.max())
        denominator = max_value - min_value
        if bool(constant_mask.loc[column]) or denominator == 0.0:
            scaled.loc[:, column] = np.nan
        else:
            scaled.loc[:, column] = (values - min_value) / denominator
    return scaled
