from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

import pandas as pd


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
    raw = Path(path).read_bytes()[:4096]
    if raw.startswith(b"\xff\xfe"):
        return "utf-16le"
    if raw.startswith(b"\xfe\xff"):
        return "utf-16be"
    if b"\x00" in raw:
        return "utf-16le"
    return "utf-8"


def read_table(path: str | Path, encoding: str | None = None) -> pd.DataFrame:
    resolved_encoding = encoding or infer_text_encoding(path)
    frame = pd.read_csv(path, sep="\t", encoding=resolved_encoding, low_memory=False)
    frame.columns = clean_columns(frame.columns)
    return frame
