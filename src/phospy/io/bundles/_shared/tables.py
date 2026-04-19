"""Table-level bundle read/write utilities."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.paths import resolve_bundle_relative_path
from phospy.io.bundles._shared.primitives import require_str
from phospy.io.readers.tables import read_table, write_table


def write_bundle_table(
    *,
    table,
    bundle_root: Path,
    relative_path: Path,
    written: dict[str, Path],
    written_key: str,
) -> str:
    """Write a table and return its manifest-safe relative path."""

    output_path = bundle_root / relative_path
    write_table(table, output_path)
    written[written_key] = output_path
    return relative_path.as_posix()


def write_optional_bundle_table(
    *,
    table,
    bundle_root: Path,
    relative_path: Path,
    written: dict[str, Path],
    written_key: str,
) -> str | None:
    """Write optional table and return manifest-safe relative path when present."""

    if table is None:
        return None
    return write_bundle_table(
        table=table,
        bundle_root=bundle_root,
        relative_path=relative_path,
        written=written,
        written_key=written_key,
    )


def read_required_table(
    *,
    bundle_root: Path,
    tables: Mapping[str, object],
    table_key: str,
    field_name: str,
):
    """Read required table declared in manifest tables section."""

    table_path = resolve_bundle_relative_path(
        bundle_root,
        require_str(tables.get(table_key), field_name=field_name),
        field_name=field_name,
    )
    return read_table(table_path)


def read_optional_table(
    *,
    bundle_root: Path,
    tables: Mapping[str, object],
    table_key: str,
    field_name: str,
):
    """Read optional table declared in manifest tables section."""

    if table_key not in tables:
        raise PhosPyInputError(f"{field_name} must be declared in the bundle manifest")
    raw_value = tables[table_key]
    if raw_value is None:
        return None
    table_path = resolve_bundle_relative_path(
        bundle_root,
        require_str(raw_value, field_name=field_name),
        field_name=field_name,
    )
    return read_table(table_path)


def read_optional_series(
    *,
    bundle_root: Path,
    tables: Mapping[str, object],
    table_key: str,
    field_name: str,
    series_name: str,
) -> pd.Series | None:
    """Read an optional single-column table and coerce it to a named Series."""

    frame = read_optional_table(
        bundle_root=bundle_root,
        tables=tables,
        table_key=table_key,
        field_name=field_name,
    )
    if frame is None:
        return None
    if frame.shape[1] != 1:
        raise PhosPyInputError(
            f"{field_name} must resolve to a single-column table for series '{series_name}'"
        )
    series = frame.iloc[:, 0].copy(deep=True)
    series.name = series_name
    return series
