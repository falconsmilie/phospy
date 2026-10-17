"""Private validation shared by concrete Kinase matrix owners."""

from __future__ import annotations

import pandas as pd

from phospy.errors.validation import PhosPyValidationError
from phospy.frames.numeric import require_numeric_unit_interval
from phospy.frames.table_schema import require_canonical_label_index
from phospy.frames.validation import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
    require_unique_index,
)
from phospy.phosphosites.validation import require_site_key_index


def _validate_kinase_score_like_matrix(
    *,
    frame: pd.DataFrame,
    field_name: str,
    error_type: type[PhosPyValidationError],
    allow_missing: bool,
    enforce_unit_interval: bool,
    require_site_index: bool,
) -> None:
    require_dataframe(
        frame,
        field_name=field_name,
        allow_empty=False,
        error_type=error_type,
    )
    require_unique_columns(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    require_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    require_finite_numeric_dataframe(
        frame,
        field_name=field_name,
        error_type=error_type,
        allow_missing=allow_missing,
    )
    require_unique_index(
        frame,
        field_name=field_name,
        error_type=error_type,
    )
    if require_site_index:
        require_site_key_index(
            frame.index,
            field_name=f"{field_name}.index",
            error_type=error_type,
        )
    else:
        require_canonical_label_index(
            frame.index,
            field_name=f"{field_name}.index",
            error_type=error_type,
        )
    require_canonical_label_index(
        frame.columns,
        field_name=f"{field_name}.columns",
        error_type=error_type,
    )
    if enforce_unit_interval:
        require_numeric_unit_interval(
            frame,
            field_name=field_name,
            error_type=error_type,
        )
