from __future__ import annotations

import numpy as np
import pandas as pd

from tests.support.performance_contracts import with_missing_fraction


def test_with_missing_fraction_sets_requested_cells_for_fortran_backed_frame() -> None:
    values = np.asfortranarray(np.arange(24, dtype=float).reshape(4, 6))
    matrix = pd.DataFrame(values)

    result = with_missing_fraction(matrix, missing_fraction=0.25, seed=123)

    assert int(result.isna().sum().sum()) == 6
    assert int(matrix.isna().sum().sum()) == 0
    assert result.shape == matrix.shape
    assert result.index.equals(matrix.index)
    assert result.columns.equals(matrix.columns)
