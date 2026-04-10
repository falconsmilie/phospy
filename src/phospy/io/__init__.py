"""Structured data access domain.

This package owns shared table reading, mapping-file loading, and path-oriented
I/O helpers that are not specific scientific behaviours. Domain logic should
not accumulate here."""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

from ..validation.schema.tables import (
    PhosphoInputSchema,
    PredMatSchema,
    TotalInputSchema,
)

DEFAULT_TEXT_ENCODING = "utf-8"


def clean_columns(columns: Iterable[str]) -> list[str]:
    cleaned: list[str] = []
    for col in columns:
        value = col.strip().lower()
        value = re.sub(r"[^0-9a-zA-Z]+", "_", value)
        value = re.sub(r"_+", "_", value)
        value = value.strip("_")
        cleaned.append(value)
    return cleaned


def default_text_encoding(path: str | Path | None = None) -> str:
    """Return the package default text encoding.

    The loader does not infer encodings from file contents. Callers should pass
    an explicit encoding when they need something other than the package default.
    The optional ``path`` argument is accepted for API convenience and backward
    compatibility with earlier helper usage.
    """

    _ = path
    return DEFAULT_TEXT_ENCODING


def infer_text_encoding(path: str | Path) -> str:
    """Backward-compatible wrapper around :func:`default_text_encoding`."""

    return default_text_encoding(path)


def read_table_raw(
    path: str | Path,
    *,
    sep: str = "\t",
    encoding: str | None = None,
    index_col: int | str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    resolved_encoding = encoding or DEFAULT_TEXT_ENCODING
    return pd.read_csv(
        path,
        sep=sep,
        encoding=resolved_encoding,
        low_memory=False,
        index_col=index_col,
        usecols=usecols,
    )


def read_table(
    path: str | Path,
    encoding: str | None = None,
    *,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = read_table_raw(path, encoding=encoding, usecols=usecols)
    return _clean_table_columns(frame)


def load_total_table(
    path: str | Path,
    *,
    encoding: str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = _clean_table_columns(
        read_table_raw(path, encoding=encoding, usecols=usecols)
    )
    return TotalInputSchema.validate(frame, context=f"total input table ({path})")


def load_phospho_table(
    path: str | Path,
    *,
    encoding: str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = _clean_table_columns(
        read_table_raw(path, encoding=encoding, usecols=usecols)
    )
    return PhosphoInputSchema.validate(frame, context=f"phospho input table ({path})")


def load_pred_mat(
    path: str | Path,
    *,
    encoding: str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = read_table_raw(
        path,
        sep=",",
        encoding=encoding,
        index_col=0,
        usecols=usecols,
    )
    frame.index = frame.index.map(str)
    frame.columns = [str(column).strip() for column in frame.columns]
    return PredMatSchema.validate(frame, context=f"pred_mat ({path})")


def load_grouped_mapping(
    path: str | Path,
    *,
    group_column: str,
    value_column: str,
    sep: str = ",",
    encoding: str | None = None,
) -> dict[str, tuple[str, ...]]:
    frame = _clean_table_columns(read_table_raw(path, sep=sep, encoding=encoding))
    group_key = clean_columns([group_column])[0]
    value_key = clean_columns([value_column])[0]
    missing = [
        column for column in (group_key, value_key) if column not in frame.columns
    ]
    if missing:
        msg = f"Grouped mapping file is missing required columns: {', '.join(missing)}"
        raise ValueError(msg)
    grouped: dict[str, list[str]] = {}
    for group, value in frame.loc[:, [group_key, value_key]].itertuples(index=False):
        grouped.setdefault(str(group).strip(), []).append(str(value).strip())
    return {key: tuple(values) for key, values in grouped.items()}


def load_string_mapping(
    path: str | Path,
    *,
    key_column: str,
    value_column: str,
    sep: str = ",",
    encoding: str | None = None,
) -> dict[str, str]:
    frame = _clean_table_columns(read_table_raw(path, sep=sep, encoding=encoding))
    key_name = clean_columns([key_column])[0]
    value_name = clean_columns([value_column])[0]
    missing = [
        column for column in (key_name, value_name) if column not in frame.columns
    ]
    if missing:
        msg = f"String mapping file is missing required columns: {', '.join(missing)}"
        raise ValueError(msg)
    return {
        str(key).strip(): str(value).strip()
        for key, value in frame.loc[:, [key_name, value_name]].itertuples(index=False)
    }


def _clean_table_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame.columns = clean_columns(str(column) for column in frame.columns)
    return frame
