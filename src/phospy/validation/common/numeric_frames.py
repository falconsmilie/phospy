"""Shared numeric DataFrame validation helpers."""

from __future__ import annotations

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
