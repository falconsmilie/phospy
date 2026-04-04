from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ._dataset_validation import _validate_dataset_file_paths, _validate_dataset_frames
from .dataset_schema import DatasetSchema
from .io import read_table
from .validation.errors import RequestValidationError, TableSchemaError


@dataclass(frozen=True, slots=True)
class _LoadedDatasetInputs:
    """Validated dataset tables produced by the internal file loader.

    The loader validates and materializes trusted in-memory frames, but it does
    not transfer ownership to downstream workspace objects automatically.
    Call :meth:`copy_frames` when you need detached caller-owned mutation.
    """

    schema: DatasetSchema
    total_df: pd.DataFrame
    phospho_df: pd.DataFrame

    def copy_frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return deep copies suitable for caller-owned mutation."""
        return self.total_df.copy(deep=True), self.phospho_df.copy(deep=True)


class _DatasetLoader:
    """Internal loader that reads and validates dataset frames from memory or disk."""

    def __init__(
        self,
        *,
        schema: DatasetSchema | None = None,
    ) -> None:
        self.schema = schema or DatasetSchema()

    def validate_inputs(
        self,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
    ) -> _LoadedDatasetInputs:
        validated_total, validated_phospho = _validate_dataset_frames(
            total_df=total_df,
            phospho_df=phospho_df,
            schema=self.schema,
        )
        return _LoadedDatasetInputs(
            total_df=validated_total,
            phospho_df=validated_phospho,
            schema=self.schema,
        )

    def load(
        self,
        total_path: str | Path,
        phospho_path: str | Path,
        *,
        phospho_encoding: str | None = None,
    ) -> _LoadedDatasetInputs:
        validated_paths = _validate_dataset_file_paths(
            total_path,
            phospho_path,
        )
        validated_total_path = validated_paths.total_path
        validated_phospho_path = validated_paths.phospho_path
        total_df = self._read_input_table(
            validated_total_path,
            context="total input table",
        )
        phospho_df = self._read_input_table(
            validated_phospho_path,
            context="phospho input table",
            encoding=phospho_encoding,
        )
        return self.validate_inputs(total_df=total_df, phospho_df=phospho_df)

    @staticmethod
    def _read_input_table(
        path: Path,
        *,
        context: str,
        encoding: str | None = None,
    ) -> pd.DataFrame:
        try:
            return read_table(path, encoding=encoding)
        except TableSchemaError:
            raise
        except (
            OSError,
            UnicodeError,
            pd.errors.ParserError,
            pd.errors.EmptyDataError,
        ) as error:
            msg = f"Invalid {context} ({path}): unable to read file: {error}"
            raise RequestValidationError(msg) from error
