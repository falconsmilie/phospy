from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

from .validation.tables import PhosphoInputSchema, PredMatSchema, TotalInputSchema


def clean_columns(columns: Iterable[str]) -> list[str]:
    cleaned: list[str] = []
    for col in columns:
        value = col.strip().lower()
        value = re.sub(r"[^0-9a-zA-Z]+", "_", value)
        value = re.sub(r"_+", "_", value)
        value = value.strip("_")
        cleaned.append(value)
    return cleaned


def infer_text_encoding(path: str | Path) -> str:
    with Path(path).open("rb") as handle:
        raw = handle.read(4096)
    if raw.startswith(b"\xff\xfe"):
        return "utf-16le"
    if raw.startswith(b"\xfe\xff"):
        return "utf-16be"
    if b"\x00" in raw:
        return "utf-16le"
    return "utf-8"


def read_table_raw(
    path: str | Path,
    *,
    sep: str = "\t",
    encoding: str | None = None,
    index_col: int | str | None = None,
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    resolved_encoding = encoding or infer_text_encoding(path)
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
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = _clean_table_columns(read_table_raw(path, usecols=usecols))
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
    usecols: Sequence[str | int] | None = None,
) -> pd.DataFrame:
    frame = read_table_raw(path, sep=",", index_col=0, usecols=usecols)
    frame.index = frame.index.map(str)
    frame.columns = [str(column).strip() for column in frame.columns]
    return PredMatSchema.validate(frame, context=f"pred_mat ({path})")


def _clean_table_columns(frame: pd.DataFrame) -> pd.DataFrame:
    frame.columns = clean_columns(str(column) for column in frame.columns)
    return frame
