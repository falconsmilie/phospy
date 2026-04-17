"""Table-level file loading and writing for supported CLI I/O."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError

_CSV = "csv"
_TSV = "tsv"
_PARQUET = "parquet"
_SUPPORTED_FORMATS = (_CSV, _TSV, _PARQUET)


def read_table(path: Path) -> pd.DataFrame:
    """Load a DataFrame from a supported table file path."""

    normalized_path = Path(path)
    table_format = table_format_from_path(normalized_path)
    try:
        if table_format == _CSV:
            return pd.read_csv(normalized_path, index_col=0)
        if table_format == _TSV:
            return pd.read_csv(normalized_path, sep="\t", index_col=0)
        return pd.read_parquet(normalized_path)
    except FileNotFoundError as exc:
        raise PhosPyInputError(f"input file does not exist: {normalized_path}") from exc
    except PermissionError as exc:
        raise PhosPyInputError(
            f"permission denied while reading input file: {normalized_path}"
        ) from exc
    except ImportError as exc:
        raise UnsupportedInputFormatError(
            "parquet input requires optional parquet dependencies (for example pyarrow)"
        ) from exc
    except (
        OSError,
        UnicodeDecodeError,
        ValueError,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
    ) as exc:
        raise PhosPyInputError(
            f"failed to parse supported table input '{normalized_path}': {exc}"
        ) from exc


def write_table(table: pd.DataFrame, path: Path) -> None:
    """Write a DataFrame to a supported table path."""

    normalized_path = Path(path)
    table_format = table_format_from_path(normalized_path)
    normalized_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        if table_format == _CSV:
            table.to_csv(normalized_path)
            return
        if table_format == _TSV:
            table.to_csv(normalized_path, sep="\t")
            return
        table.to_parquet(normalized_path)
    except ImportError as exc:
        raise UnsupportedInputFormatError(
            "parquet output requires optional parquet dependencies (for example pyarrow)"
        ) from exc
    except (OSError, ValueError) as exc:
        raise PhosPyInputError(
            f"failed to write output table '{normalized_path}': {exc}"
        ) from exc


def table_format_from_path(path: Path) -> str:
    """Infer supported table format from file suffix."""

    suffix = path.suffix.lower()
    if suffix == ".csv":
        return _CSV
    if suffix in {".tsv", ".txt"}:
        return _TSV
    if suffix == ".parquet":
        return _PARQUET
    supported = ", ".join(_SUPPORTED_FORMATS)
    raise UnsupportedInputFormatError(
        f"unsupported input format for '{path}'. supported formats: {supported}"
    )


def table_suffix_for_format(output_format: str) -> str:
    """Resolve an output format token to a file suffix."""

    normalized = output_format.strip().lower()
    if normalized == _CSV:
        return ".csv"
    if normalized == _TSV:
        return ".tsv"
    if normalized == _PARQUET:
        return ".parquet"
    supported = ", ".join(_SUPPORTED_FORMATS)
    raise UnsupportedInputFormatError(
        f"unsupported output format '{output_format}'. supported formats: {supported}"
    )
