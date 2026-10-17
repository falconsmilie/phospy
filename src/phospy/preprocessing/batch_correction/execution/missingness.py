"""Array-oriented internals shared by batch-correction preparation paths."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def materialize_originally_missing(
    *,
    index: pd.Index,
    columns: pd.Index,
    feature_ids: Sequence[str],
    sample_ids: Sequence[str],
    missing_cells: Sequence[tuple[str, str]],
) -> pd.DataFrame:
    """Materialize governed missing coordinates after callers validate axes."""
    cell_count = len(missing_cells)
    if cell_count == 0:
        return pd.DataFrame(False, index=index.copy(), columns=columns.copy())
    feature_positions = {
        feature_id: position for position, feature_id in enumerate(feature_ids)
    }
    sample_positions = {
        sample_id: position for position, sample_id in enumerate(sample_ids)
    }
    column_count = len(columns)
    flat_positions = np.fromiter(
        (
            feature_positions[feature_id] * column_count + sample_positions[sample_id]
            for feature_id, sample_id in missing_cells
        ),
        dtype=np.intp,
        count=cell_count,
    )
    materialized = np.zeros((len(index), len(columns)), dtype=bool)
    materialized.ravel()[flat_positions] = True
    return pd.DataFrame(materialized, index=index.copy(), columns=columns.copy())
