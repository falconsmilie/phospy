from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..errors import RequestValidationError, TableSchemaError
from ..io import read_table
from ..validation.requests import validate_dataset_file_paths, validate_dataset_frames
from ..validation.schema.files import validate_existing_file_path
from ..validation.schema.tables import PhosphoInputSchema, TotalInputSchema
from .schema import DatasetSchema

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

    def validate_total(self, total_df: pd.DataFrame) -> pd.DataFrame:
        """Validate one in-memory total-proteome input table."""

        return TotalInputSchema.validate(total_df, total_cols=self.schema.total_cols)

    def validate_phospho(self, phospho_df: pd.DataFrame) -> pd.DataFrame:
        """Validate one in-memory phosphoproteome input table."""

        return PhosphoInputSchema.validate(
            phospho_df,
            phospho_cols=self.schema.phospho_cols,
        )

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

    def _load_total_from_validated_path(
        self,
        path: Path,
        *,
        encoding: str | None = None,
    ) -> pd.DataFrame:
        frame = self._read_input_table(
            path,
            context="total input table",
            encoding=encoding,
        )
        return self.validate_total(frame)

    def _load_phospho_from_validated_path(
        self,
        path: Path,
        *,
        encoding: str | None = None,
    ) -> pd.DataFrame:
        frame = self._read_input_table(
            path,
            context="phospho input table",
            encoding=encoding,
        )
        return self.validate_phospho(frame)

    def load_total(
        self,
        total_path: str | Path,
        *,
        encoding: str | None = None,
    ) -> pd.DataFrame:
        """Load and validate one total-proteome input table from disk."""

        validated_path = validate_existing_file_path(
            total_path,
            context="total input table path",
        )
        return self._load_total_from_validated_path(
            validated_path,
            encoding=encoding,
        )

    def load_phospho(
        self,
        phospho_path: str | Path,
        *,
        encoding: str | None = None,
    ) -> pd.DataFrame:
        """Load and validate one phosphoproteome input table from disk."""

        validated_path = validate_existing_file_path(
            phospho_path,
            context="phospho input table path",
        )
        return self._load_phospho_from_validated_path(
            validated_path,
            encoding=encoding,
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
        total_df = self._load_total_from_validated_path(validated_paths.total_path)
        phospho_df = self._load_phospho_from_validated_path(
            validated_paths.phospho_path,
            encoding=phospho_encoding,
        )
        return LoadedDatasetInputs(
            total_df=total_df,
            phospho_df=phospho_df,
            schema=self.schema,
        )

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
