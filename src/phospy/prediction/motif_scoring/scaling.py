from __future__ import annotations

import numpy as np
import pandas as pd


def minmax_scale_columns(mat: pd.DataFrame) -> pd.DataFrame:
    """Apply column-wise min-max scaling used by baseline motif scoring."""

    scaled = mat.astype(float).copy()
    for column in scaled.columns:
        values = scaled.loc[:, column]
        min_value = float(values.min())
        max_value = float(values.max())
        denominator = max_value - min_value
        if denominator == 0.0:
            scaled.loc[:, column] = np.nan
        else:
            scaled.loc[:, column] = (values - min_value) / denominator
    return scaled
