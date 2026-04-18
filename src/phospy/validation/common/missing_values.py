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
    """Validate missing-value policy against a numeric DataFrame."""

    if policy is MissingValuePolicy.ALLOW:
        return value
    if value.isna().to_numpy().any():
        raise error_type(f"{field_name} must not contain missing values")
    if not np.isfinite(value.to_numpy(dtype=float, copy=False)).all():
        raise error_type(f"{field_name} must contain finite numeric values")
    return value
