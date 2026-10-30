"""Typed pandas content-comparison helpers.

These helpers intentionally compare only pandas leaves. Domain containers own
the decision about which fields are scientifically meaningful.
"""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd


def dataframe_equals(left: pd.DataFrame, right: pd.DataFrame) -> bool:
    """Return ``True`` when two DataFrames have identical pandas content."""

    return left.equals(right)


def optional_dataframe_equals(
    left: pd.DataFrame | None,
    right: pd.DataFrame | None,
) -> bool:
    """Return ``True`` when optional DataFrames are both absent or equal."""

    if left is None or right is None:
        return left is right
    return left.equals(right)


def series_equals(left: pd.Series, right: pd.Series) -> bool:
    """Return ``True`` when two Series have identical pandas content."""

    return left.equals(right)


def optional_series_equals(left: pd.Series | None, right: pd.Series | None) -> bool:
    """Return ``True`` when optional Series are both absent or equal."""

    if left is None or right is None:
        return left is right
    return left.equals(right)


def index_equals(left: pd.Index, right: pd.Index) -> bool:
    """Return ``True`` when two Index objects have identical pandas content."""

    return left.equals(right)


def optional_index_equals(left: pd.Index | None, right: pd.Index | None) -> bool:
    """Return ``True`` when optional Index objects are both absent or equal."""

    if left is None or right is None:
        return left is right
    return left.equals(right)


def dataframe_mapping_equals(
    left: Mapping[str, pd.DataFrame],
    right: Mapping[str, pd.DataFrame],
) -> bool:
    """Return ``True`` when string-keyed DataFrame mappings have equal content."""

    if set(left.keys()) != set(right.keys()):
        return False
    return all(left[key].equals(right[key]) for key in left)


__all__ = [
    "dataframe_equals",
    "dataframe_mapping_equals",
    "index_equals",
    "optional_dataframe_equals",
    "optional_index_equals",
    "optional_series_equals",
    "series_equals",
]
