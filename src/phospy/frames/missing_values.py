"""Shared missing-value policy validation."""

from __future__ import annotations

from enum import Enum

import numpy as np
import pandas as pd

from phospy.errors.validation import PhosPyValidationError

ValidationErrorType = type[PhosPyValidationError]


class MissingValuePolicy(str, Enum):
    """Policy for missing/non-finite numeric values."""

    ALLOW = "allow"
    FORBID = "forbid"


def require_missing_value_policy(
    value: pd.DataFrame,
    *,
    field_name: str,
    policy: MissingValuePolicy,
    error_type: ValidationErrorType,
) -> pd.DataFrame:
    """Validate missing-value policy against a numeric DataFrame.

    ``ALLOW`` permits missing values but still rejects +/-inf.
    ``FORBID`` rejects both missing values and +/-inf.
    """

    missing_mask = value.isna()
    if policy is MissingValuePolicy.FORBID and missing_mask.to_numpy().any():
        raise error_type(
            f"{field_name} must not contain missing values; "
            f"{_invalid_location_preview(missing_mask)}"
        )

    infinite_mask = pd.DataFrame(
        np.isinf(value.to_numpy(dtype=float, copy=False)),
        index=value.index,
        columns=value.columns,
    )
    if infinite_mask.to_numpy().any():
        raise error_type(
            f"{field_name} must contain finite numeric values; "
            f"{_invalid_location_preview(infinite_mask)}"
        )
    return value


def _invalid_location_preview(mask: pd.DataFrame, *, max_items: int = 3) -> str:
    locations = np.argwhere(mask.to_numpy())
    count = int(locations.shape[0])
    preview = [
        f"({mask.index[row_idx]!r}, {mask.columns[col_idx]!r})"
        for row_idx, col_idx in locations[:max_items]
    ]
    suffix = "" if count <= max_items else f", +{count - max_items} more"
    return f"found at {', '.join(preview)}{suffix}"
