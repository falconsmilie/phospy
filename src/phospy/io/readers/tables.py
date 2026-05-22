"""Table-level file loading and writing for supported CLI I/O."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, NoReturn

import pandas as pd

from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError

_CSV = "csv"
_TSV = "tsv"
_PARQUET = "parquet"
_SUPPORTED_OUTPUT_FORMATS = (_CSV, _TSV, _PARQUET)
_SUPPORTED_INPUT_FORMATS_LABEL = (
    "csv (.csv), tsv (.tsv), txt as tab-separated tsv (.txt), parquet (.parquet)"
)
_MISSING_NUMERIC_TOKENS = frozenset({"", "na", "n/a", "nan", "null"})


def supported_table_input_formats() -> str:
    """Return the user-facing list of supported table input formats."""

    return _SUPPORTED_INPUT_FORMATS_LABEL


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
            f"failed to parse table input '{normalized_path}': {exc}"
        ) from exc


def read_phospho_matrix(path: Path) -> pd.DataFrame:
    """Load phospho matrix with explicit numeric parsing and missing support."""

    return _read_numeric_table(
        path,
        table_role="phospho_matrix",
        allow_missing=True,
    )


def read_total_matrix(path: Path) -> pd.DataFrame:
    """Load total-proteome matrix with explicit numeric parsing and missing support."""

    return _read_numeric_table(
        path,
        table_role="total_matrix",
        allow_missing=True,
    )


def read_design_matrix(path: Path) -> pd.DataFrame:
    """Load differential design matrix as strict finite numeric values."""

    return _read_numeric_table(
        path,
        table_role="design_matrix",
        allow_missing=False,
    )


def read_contrast_matrix(path: Path) -> pd.DataFrame:
    """Load differential contrast matrix as strict finite numeric values."""

    return _read_numeric_table(
        path,
        table_role="contrast_matrix",
        allow_missing=False,
    )


def read_site_metadata(path: Path) -> pd.DataFrame:
    """Load site metadata while preserving identifiers as explicit strings."""

    return _read_metadata_table(path, table_role="site_metadata")


def read_sample_metadata(path: Path) -> pd.DataFrame:
    """Load sample metadata while preserving identifiers as explicit strings."""

    return _read_metadata_table(path, table_role="sample_metadata")


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
    raise UnsupportedInputFormatError(
        "unsupported table file format for "
        f"'{path}'. supported formats: {supported_table_input_formats()}"
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
    supported = ", ".join(_SUPPORTED_OUTPUT_FORMATS)
    raise UnsupportedInputFormatError(
        f"unsupported output format '{output_format}'. supported formats: {supported}"
    )


def _read_metadata_table(path: Path, *, table_role: str) -> pd.DataFrame:
    frame = _read_table_with_policy(
        path,
        table_role=table_role,
        dtype=str,
        keep_default_na=False,
        na_values=[],
    )
    return _stringify_dataframe(frame)


def _read_numeric_table(
    path: Path,
    *,
    table_role: str,
    allow_missing: bool,
) -> pd.DataFrame:
    raw = _read_table_with_policy(
        path,
        table_role=table_role,
        dtype=str,
        keep_default_na=False,
        na_values=[],
    )
    raw = _stringify_index_and_columns(raw)
    return _parse_numeric_cells(
        raw,
        source_path=Path(path),
        table_role=table_role,
        allow_missing=allow_missing,
    )


def _read_table_with_policy(
    path: Path,
    *,
    table_role: str,
    dtype: type[str] | None,
    keep_default_na: bool,
    na_values: list[str],
) -> pd.DataFrame:
    normalized_path = Path(path)
    table_format = table_format_from_path(normalized_path)
    try:
        if table_format == _CSV:
            return _read_delimited_table_with_policy(
                normalized_path,
                sep=",",
                dtype=dtype,
                keep_default_na=keep_default_na,
                na_values=na_values,
            )
        if table_format == _TSV:
            return _read_delimited_table_with_policy(
                normalized_path,
                sep="\t",
                dtype=dtype,
                keep_default_na=keep_default_na,
                na_values=na_values,
            )
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
            f"failed to parse {table_role} table input '{normalized_path}': {exc}"
        ) from exc


def _read_delimited_table_with_policy(
    path: Path,
    *,
    sep: str,
    dtype: type[str] | None,
    keep_default_na: bool,
    na_values: list[str],
) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        sep=sep,
        dtype=dtype,
        keep_default_na=keep_default_na,
        na_values=na_values,
    )
    if frame.shape[1] == 0:
        return frame
    index_column = frame.columns[0]
    normalized_index_name = (
        None
        if isinstance(index_column, str) and index_column.startswith("Unnamed:")
        else index_column
    )
    index = pd.Index(frame.iloc[:, 0].tolist(), name=normalized_index_name)
    return pd.DataFrame(
        frame.iloc[:, 1:].to_numpy(copy=True),
        index=index,
        columns=frame.columns[1:].copy(),
    )


def _stringify_dataframe(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = _stringify_index_and_columns(frame)
    return normalized.map(_stringify_value)


def _stringify_index_and_columns(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.copy(deep=True)
    normalized.index = pd.Index(
        [_stringify_value(label) for label in normalized.index],
        name=normalized.index.name,
    )
    normalized.columns = pd.Index(
        [_stringify_value(label) for label in normalized.columns],
        name=normalized.columns.name,
    )
    return normalized


def _stringify_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        missing = False
    if isinstance(missing, bool) and missing:
        return ""
    return str(value)


def _parse_numeric_cells(
    frame: pd.DataFrame,
    *,
    source_path: Path,
    table_role: str,
    allow_missing: bool,
) -> pd.DataFrame:
    missing_tokens = _MISSING_NUMERIC_TOKENS
    parsed_rows: list[list[float]] = []
    for row_label, row in frame.iterrows():
        parsed_row: list[float] = []
        for column_label, raw_value in row.items():
            value = _stringify_value(raw_value)
            normalized_value = value.strip().lower()
            if allow_missing and normalized_value in missing_tokens:
                parsed_row.append(float("nan"))
                continue
            try:
                numeric_value = float(value)
            except ValueError as exc:
                _raise_numeric_cell_error(
                    source_path=source_path,
                    table_role=table_role,
                    row_label=row_label,
                    column_label=column_label,
                    offending_value=value,
                    allow_missing=allow_missing,
                    original_error=exc,
                )
            if not math.isfinite(numeric_value):
                _raise_numeric_cell_error(
                    source_path=source_path,
                    table_role=table_role,
                    row_label=row_label,
                    column_label=column_label,
                    offending_value=value,
                    allow_missing=allow_missing,
                    original_error=ValueError("non-finite numeric value"),
                )
            parsed_row.append(numeric_value)
        parsed_rows.append(parsed_row)
    return pd.DataFrame(parsed_rows, index=frame.index.copy(), columns=frame.columns)


def _raise_numeric_cell_error(
    *,
    source_path: Path,
    table_role: str,
    row_label: object,
    column_label: object,
    offending_value: str,
    allow_missing: bool,
    original_error: Exception,
) -> NoReturn:
    expected_type = "finite numeric value"
    if allow_missing:
        expected_type = (
            "finite numeric value or allowed missing marker "
            f"{tuple(sorted(_MISSING_NUMERIC_TOKENS))}"
        )
    raise PhosPyInputError(
        "failed to parse numeric cell: "
        f"path='{source_path}', table_role='{table_role}', "
        f"row_label='{row_label}', column_label='{column_label}', "
        f"offending_value={offending_value!r}, expected_type='{expected_type}'"
    ) from original_error
