"""Strict missing-data rejection policy."""

from __future__ import annotations

from typing import cast

import pandas as pd

from phospy.datasets.processing_state import MissingDataDiagnosticsV1
from phospy.errors.input import PhosPyInputError

from .diagnostics import label_preview


def fail_if_forbid_policy_has_missing_values(
    phospho: pd.DataFrame,
    *,
    diagnostics: MissingDataDiagnosticsV1 | None = None,
) -> None:
    """Raise if any missing values are present for forbid policy."""

    missing_mask = phospho.isna()
    missing_cell_count = int(missing_mask.to_numpy().sum())
    if missing_cell_count == 0:
        return

    affected_row_count = int(missing_mask.any(axis=1).sum())
    affected_column_count = int(missing_mask.any(axis=0).sum())
    affected_row_preview = label_preview(
        [
            cast(object, str(label))
            for label in phospho.index[missing_mask.any(axis=1)].tolist()
        ]
    )
    affected_column_preview = label_preview(
        [
            cast(object, str(label))
            for label in phospho.columns[missing_mask.any(axis=0)].tolist()
        ]
    )
    raise PhosPyInputError(
        "dataset preprocessing stage 'missing_data' rejected missing values because "
        "missing_data.policy='forbid'; "
        f"found {missing_cell_count} missing values across "
        f"{affected_row_count} rows and {affected_column_count} columns. "
        f"affected row labels (preview): {affected_row_preview}. "
        f"affected column labels (preview): {affected_column_preview}. "
        "choose missing_data.policy='impute_row_median' or clean the input data.",
        diagnostics=diagnostics,
    )
