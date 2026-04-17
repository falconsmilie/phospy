from __future__ import annotations

from pathlib import Path

from .readers import clean_columns, clean_table_columns, read_table_raw


def load_grouped_mapping(
    path: str | Path,
    *,
    group_column: str,
    value_column: str,
    sep: str = ",",
    encoding: str | None = None,
) -> dict[str, tuple[str, ...]]:
    frame = clean_table_columns(read_table_raw(path, sep=sep, encoding=encoding))
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
    frame = clean_table_columns(read_table_raw(path, sep=sep, encoding=encoding))
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


__all__ = ["load_grouped_mapping", "load_string_mapping"]
