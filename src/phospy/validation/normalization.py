from __future__ import annotations

import pandas as pd


def normalize_identifier_series(series: pd.Series) -> pd.Series:
    """Normalize identifier values for case/whitespace-insensitive joins.

    Values are converted to pandas ``string`` dtype, trimmed, and uppercased so
    validation and execution share identical join semantics.
    """
    return series.astype("string").str.strip().str.upper()
