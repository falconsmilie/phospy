"""Shared numeric DataFrame validation helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.validation.common.dataframes import require_dataframe
from phospy.validation.common.missing_values import (
    MissingValuePolicy,
    require_missing_value_policy,
)

ValidationErrorType = type[PhosPyValidationError]


def require_numeric_matrix(
    value: object,
    *,
    field_name: str,
    allow_empty: bool,
    missing_value_policy: MissingValuePolicy,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require a DataFrame with numeric values and a missing-value policy."""

    frame = require_dataframe(
        value,
        field_name=field_name,
        allow_empty=allow_empty,
        error_type=error_type,
    )
    try:
        numeric_frame = frame.astype(float)
    except (TypeError, ValueError) as exc:
        raise error_type(f"{field_name} must contain numeric values") from exc
    require_missing_value_policy(
        numeric_frame,
        field_name=field_name,
        policy=missing_value_policy,
        error_type=error_type,
    )
    return numeric_frame


def require_numeric_unit_interval(
    value: pd.DataFrame,
    *,
    field_name: str,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Require finite numeric values in ``[0.0, 1.0]`` for non-missing entries."""

    values = value.to_numpy(dtype=float, copy=False)
    finite_mask = np.isfinite(values)
    invalid_mask = finite_mask & ((values < 0.0) | (values > 1.0))
    if not invalid_mask.any():
        return value
    bad_positions = np.argwhere(invalid_mask)
    preview_items: list[str] = []
    for row_idx, col_idx in bad_positions[:3]:
        preview_items.append(
            f"({value.index[row_idx]!r}, {value.columns[col_idx]!r})={values[row_idx, col_idx]}"
        )
    suffix = (
        "" if bad_positions.shape[0] <= 3 else f", +{bad_positions.shape[0] - 3} more"
    )
    raise error_type(
        f"{field_name} must contain scores between 0.0 and 1.0; "
        f"found at {', '.join(preview_items)}{suffix}"
    )
