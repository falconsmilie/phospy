from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .dataset_schema import DatasetSchema
from .io import read_table
from .validation.errors import RequestValidationError, TableSchemaError
from .validation.paths import validate_existing_file_path
from .validation.tables import PhosphoInputSchema, TotalInputSchema


@dataclass(frozen=True, slots=True)
class ValidatedCoreInputs:
    """Validated dataset tables produced by :class:`DatasetLoader`.

    The loader is the boundary owner for validated dataset frames. The stored
    frames are the single trusted in-memory snapshot that downstream code may
    read directly. Call :meth:`copy_frames` before mutating them outside the
    dataset-bound processing path.
    """

    schema: DatasetSchema
    total_df: pd.DataFrame
    phospho_df: pd.DataFrame

    def copy_frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return deep copies suitable for caller-owned mutation."""
        return self.total_df.copy(deep=True), self.phospho_df.copy(deep=True)


class DatasetLoader:
    """Load and validate dataset frames from memory or disk."""

    def __init__(
        self,
        *,
        schema: DatasetSchema | None = None,
    ) -> None:
        self.schema = schema or DatasetSchema()

    def validate(
        self,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
    ) -> ValidatedCoreInputs:
        validated_total = TotalInputSchema.validate(
            total_df,
            total_cols=self.schema.total_cols,
        )
        validated_phospho = PhosphoInputSchema.validate(
            phospho_df,
            phospho_cols=self.schema.phospho_cols,
        )
        return ValidatedCoreInputs(
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
    ) -> ValidatedCoreInputs:
        validated_total_path = validate_existing_file_path(
            total_path,
            context="total input table path",
        )
        validated_phospho_path = validate_existing_file_path(
            phospho_path,
            context="phospho input table path",
        )
        total_df = self._read_input_table(
            validated_total_path,
            context="total input table",
        )
        phospho_df = self._read_input_table(
            validated_phospho_path,
            context="phospho input table",
            encoding=phospho_encoding,
        )
        return self.validate(total_df=total_df, phospho_df=phospho_df)

    def _read_input_table(
        self,
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
