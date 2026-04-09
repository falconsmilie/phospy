from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .dataset_schema import DatasetSchema
from .io import read_table
from .validation.errors import RequestValidationError, TableSchemaError
from .validation.requests import validate_dataset_file_paths, validate_dataset_frames

__all__ = ["DatasetLoader", "LoadedDatasetInputs"]


@dataclass(slots=True)
class LoadedDatasetInputs:
    """Package-internal bundle produced by :class:`DatasetLoader`.

    The loader materializes owned validated frames exactly once. Downstream
    package-internal constructors may take ownership of these frames directly
    instead of copying them again.
    """

    schema: DatasetSchema
    total_df: pd.DataFrame
    phospho_df: pd.DataFrame


class DatasetLoader:
    """Package-internal loader contract for reading validated dataset frames."""

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
    ) -> LoadedDatasetInputs:
        validated_total, validated_phospho = validate_dataset_frames(
            total_df=total_df,
            phospho_df=phospho_df,
            schema=self.schema,
        )
        return LoadedDatasetInputs(
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
    ) -> LoadedDatasetInputs:
        validated_paths = validate_dataset_file_paths(
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
