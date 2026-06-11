"""Local source-table loading for reference bundle construction."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from phospy.errors.input import PhosPyInputError, UnsupportedInputFormatError
from phospy.io.readers.tables import (
    supported_table_input_formats,
    table_format_from_path,
)


class ReferenceSourceTableReader:
    """Read local tabular reference source files without treating column 0 as index."""

    def run(self, path: Path, *, field_name: str) -> pd.DataFrame:
        normalized_path = Path(path)
        table_format = table_format_from_path(normalized_path)
        try:
            if table_format == "csv":
                return self._stringify_dataframe(
                    pd.read_csv(
                        normalized_path,
                        dtype=str,
                        keep_default_na=False,
                        na_values=[],
                    )
                )
            if table_format == "tsv":
                return self._stringify_dataframe(
                    pd.read_csv(
                        normalized_path,
                        sep="\t",
                        dtype=str,
                        keep_default_na=False,
                        na_values=[],
                    )
                )
            frame = pd.read_parquet(normalized_path)
            if not isinstance(frame.index, pd.RangeIndex):
                frame = frame.reset_index()
            return self._stringify_dataframe(frame)
        except FileNotFoundError as exc:
            raise PhosPyInputError(
                f"{field_name} source file does not exist: {normalized_path}"
            ) from exc
        except PermissionError as exc:
            raise PhosPyInputError(
                f"permission denied while reading {field_name} source file: "
                f"{normalized_path}"
            ) from exc
        except UnsupportedInputFormatError:
            raise
        except ImportError as exc:
            raise UnsupportedInputFormatError(
                "parquet reference source input requires optional parquet "
                "dependencies (for example pyarrow)"
            ) from exc
        except (
            OSError,
            UnicodeDecodeError,
            ValueError,
            pd.errors.EmptyDataError,
            pd.errors.ParserError,
        ) as exc:
            raise PhosPyInputError(
                f"failed to parse {field_name} source file '{normalized_path}': {exc}"
            ) from exc

    @staticmethod
    def supported_input_formats() -> str:
        return supported_table_input_formats()

    @classmethod
    def _stringify_dataframe(cls, frame: pd.DataFrame) -> pd.DataFrame:
        normalized = frame.copy(deep=True)
        normalized.columns = pd.Index(
            [cls._stringify_value(column) for column in normalized.columns]
        )
        return normalized.map(cls._stringify_value)

    @staticmethod
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
