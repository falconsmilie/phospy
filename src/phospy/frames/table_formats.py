"""Pure table-format recognition shared by file adapters."""

from __future__ import annotations

from pathlib import Path

from phospy.errors.input import UnsupportedInputFormatError

TABLE_FORMAT_CSV = "csv"
TABLE_FORMAT_TSV = "tsv"
TABLE_FORMAT_PARQUET = "parquet"
SUPPORTED_TABLE_OUTPUT_FORMATS = (
    TABLE_FORMAT_CSV,
    TABLE_FORMAT_TSV,
    TABLE_FORMAT_PARQUET,
)
_SUPPORTED_INPUT_FORMATS_LABEL = (
    "csv (.csv), tsv (.tsv), txt as tab-separated tsv (.txt), parquet (.parquet)"
)


def supported_table_input_formats() -> str:
    """Return the user-facing list of supported table input formats."""

    return _SUPPORTED_INPUT_FORMATS_LABEL


def table_format_from_path(path: Path) -> str:
    """Infer a supported table format from a file suffix."""

    suffix = path.suffix.lower()
    if suffix == ".csv":
        return TABLE_FORMAT_CSV
    if suffix in {".tsv", ".txt"}:
        return TABLE_FORMAT_TSV
    if suffix == ".parquet":
        return TABLE_FORMAT_PARQUET
    raise UnsupportedInputFormatError(
        "unsupported table file format for "
        f"'{path}'. supported formats: {supported_table_input_formats()}"
    )


__all__ = [
    "SUPPORTED_TABLE_OUTPUT_FORMATS",
    "TABLE_FORMAT_CSV",
    "TABLE_FORMAT_PARQUET",
    "TABLE_FORMAT_TSV",
    "supported_table_input_formats",
    "table_format_from_path",
]
