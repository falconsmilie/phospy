"""Reusable dataframe column validation mechanics."""

from __future__ import annotations

import math
from typing import TypeVar, cast

import numpy as np
import pandas as pd

ErrorType = TypeVar("ErrorType", bound=Exception)


def enforce_required_non_empty_string_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    workflow_name: str,
    column_name: str,
    error_type: type[ErrorType],
) -> None:
    """Require one dataframe column to be non-missing, non-empty strings."""

    if column_name not in site_metadata.columns:
        site_examples = _site_id_examples(site_metadata.index)
        raise error_type(
            f"{field_name} is missing required columns: {column_name}; "
            f"{workflow_name} requires {field_name}.{column_name}; "
            f"missing required column={field_name}.{column_name}; "
            f"affected_rows={int(site_metadata.shape[0])}; "
            f"example_site_ids={site_examples}"
        )
    column = site_metadata[column_name]
    invalid_mask = pd.Series(False, index=pd.Index(column.index), dtype="boolean")
    for site_id, raw_value in column.items():
        if _is_missing(raw_value):
            invalid_mask.at[site_id] = True
            continue
        if not isinstance(raw_value, str):
            invalid_mask.at[site_id] = True
            continue
        if raw_value.strip() == "":
            invalid_mask.at[site_id] = True
    invalid_count = int(invalid_mask.sum())
    if invalid_count == 0:
        return
    invalid_sites = _site_id_examples(
        _index_by_boolean_mask(site_metadata.index, invalid_mask)
    )
    raise error_type(
        f"{workflow_name} requires {field_name}.{column_name} to contain non-empty "
        f"string values; affected_rows={invalid_count}; "
        f"example_site_ids={invalid_sites}"
    )


def _site_id_examples(index: pd.Index, *, limit: int = 5) -> str:
    labels = [str(value) for value in index.tolist()]
    if not labels:
        return "(none)"
    preview = ", ".join(repr(label) for label in labels[:limit])
    suffix = "" if len(labels) <= limit else f", +{len(labels) - limit} more"
    return f"[{preview}{suffix}]"


def _index_by_boolean_mask(index: pd.Index, mask: pd.Series) -> pd.Index:
    labels = index.tolist()
    mask_values = cast(list[object], mask.tolist())
    return pd.Index(
        label
        for label, include in zip(labels, mask_values, strict=True)
        if bool(include)
    )


def _is_missing(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value = cast(object, value)
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value = cast(object, value)
        return str(temporal_value) == "NaT"
    return False


__all__ = ["enforce_required_non_empty_string_column"]
