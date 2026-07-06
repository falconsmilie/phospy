"""Typed pandas boundary helpers for strict workflow orchestration modules."""

from __future__ import annotations

from collections.abc import Hashable
from typing import cast

import pandas as pd


def dataframe_column(frame: pd.DataFrame, column: Hashable) -> pd.Series:
    result = frame.loc[:, column]  # pyright: ignore[reportUnknownMemberType, reportCallIssue, reportArgumentType] - pandas-stubs cannot express scalar-column DataFrame.loc selection for generic frames.
    return cast(pd.Series, result)


def dataframe_copy(frame: pd.DataFrame, *, deep: bool = True) -> pd.DataFrame:
    result = frame.copy(deep=deep)  # pyright: ignore[reportUnknownMemberType] - pandas-stubs loses the concrete DataFrame type for copy on unparameterized frames.
    return cast(pd.DataFrame, result)


def dataframe_loc(
    frame: pd.DataFrame,
    rows: object = slice(None),
    columns: object = slice(None),
) -> pd.DataFrame:
    result = frame.loc[rows, columns]  # pyright: ignore[reportUnknownMemberType, reportCallIssue, reportArgumentType] - pandas-stubs cannot express runtime-valid generic DataFrame.loc indexers.
    return cast(pd.DataFrame, result)


def dataframe_reindex(frame: pd.DataFrame, index: pd.Index) -> pd.DataFrame:
    result = frame.reindex(index)  # pyright: ignore[reportUnknownMemberType] - pandas-stubs returns Unknown for generic DataFrame.reindex.
    return cast(pd.DataFrame, result)


def dataframe_reset_index(frame: pd.DataFrame, *, drop: bool = False) -> pd.DataFrame:
    result = frame.reset_index(drop=drop)  # pyright: ignore[reportUnknownMemberType] - pandas-stubs loses the concrete DataFrame type for reset_index on generic frames.
    return cast(pd.DataFrame, result)


def index_as_strings(index: pd.Index) -> list[str]:
    return [str(value) for value in index]


def index_snapshot(index: pd.Index, *, name: Hashable | None = None) -> pd.Index:
    return pd.Index(list(index), name=index.name if name is None else name)


def series_as_strings(
    series: pd.Series,
    *,
    fill_missing: str | None = None,
    strip: bool = False,
) -> list[str]:
    values: list[str] = []
    for value in series.to_numpy(dtype=object):
        text = (
            fill_missing
            if fill_missing is not None and bool(pd.isna(value))
            else str(value)
        )
        values.append(text.strip() if strip else text)
    return values


def series_copy(series: pd.Series, *, deep: bool = True) -> pd.Series:
    result = series.copy(deep=deep)  # pyright: ignore[reportUnknownMemberType] - pandas-stubs loses the concrete Series type for copy on unparameterized series.
    return cast(pd.Series, result)


__all__ = [
    "dataframe_column",
    "dataframe_copy",
    "dataframe_loc",
    "dataframe_reindex",
    "dataframe_reset_index",
    "index_as_strings",
    "index_snapshot",
    "series_as_strings",
    "series_copy",
]
