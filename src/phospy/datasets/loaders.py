from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..errors import RequestValidationError, TableSchemaError
from ..io import read_table
from ..validation.requests.dataset import (
    validate_dataset_file_paths,
    validate_dataset_frames,
)
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

        return self._validate_table(
            total_df,
            validator=lambda frame: TotalInputSchema.validate(
                frame,
                total_cols=self.schema.total_cols,
            ),
        )

    def validate_phospho(self, phospho_df: pd.DataFrame) -> pd.DataFrame:
        """Validate one in-memory phosphoproteome input table."""

        return self._validate_table(
            phospho_df,
            validator=lambda frame: PhosphoInputSchema.validate(
                frame,
                phospho_cols=self.schema.phospho_cols,
            ),
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
        return self._build_loaded_inputs(
            total_df=validated_total,
            phospho_df=validated_phospho,
        )

    def resolve_total(self, total: pd.DataFrame | str | Path) -> pd.DataFrame:
        """Resolve one total input from memory or disk into a validated frame."""

        if isinstance(total, pd.DataFrame):
            return self.validate_total(total)
        return self.load_total(total)

    def resolve_phospho(
        self,
        phospho: pd.DataFrame | str | Path,
        *,
        encoding: str | None = None,
    ) -> pd.DataFrame:
        """Resolve one phospho input from memory or disk into a validated frame."""

        if isinstance(phospho, pd.DataFrame):
            return self.validate_phospho(phospho)
        return self.load_phospho(phospho, encoding=encoding)

    def resolve_inputs(
        self,
        *,
        total: pd.DataFrame | str | Path,
        phospho: pd.DataFrame | str | Path,
        phospho_encoding: str | None = None,
    ) -> LoadedDatasetInputs:
        """Resolve dataset inputs from file-backed, in-memory, or mixed sources."""

        if isinstance(total, pd.DataFrame) and isinstance(phospho, pd.DataFrame):
            return self.validate_inputs(total_df=total, phospho_df=phospho)
        if not isinstance(total, pd.DataFrame) and not isinstance(
            phospho, pd.DataFrame
        ):
            return self.load(
                total,
                phospho,
                phospho_encoding=phospho_encoding,
            )
        return self._build_loaded_inputs(
            total_df=self.resolve_total(total),
            phospho_df=self.resolve_phospho(
                phospho,
                encoding=phospho_encoding,
            ),
        )

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
        return self._load_from_validated_path(
            validated_path,
            context="total input table",
            validator=self.validate_total,
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
        return self._load_from_validated_path(
            validated_path,
            context="phospho input table",
            validator=self.validate_phospho,
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
        return self._build_loaded_inputs(
            total_df=self._load_from_validated_path(
                validated_paths.total_path,
                context="total input table",
                validator=self.validate_total,
            ),
            phospho_df=self._load_from_validated_path(
                validated_paths.phospho_path,
                context="phospho input table",
                validator=self.validate_phospho,
                encoding=phospho_encoding,
            ),
        )

    def _build_loaded_inputs(
        self,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
    ) -> LoadedDatasetInputs:
        return LoadedDatasetInputs(
            total_df=total_df,
            phospho_df=phospho_df,
            schema=self.schema,
        )

    @staticmethod
    def _validate_table(
        frame: pd.DataFrame,
        *,
        validator: Callable[[pd.DataFrame], pd.DataFrame],
    ) -> pd.DataFrame:
        return validator(frame)

    def _load_from_validated_path(
        self,
        path: Path,
        *,
        context: str,
        validator: Callable[[pd.DataFrame], pd.DataFrame],
        encoding: str | None = None,
    ) -> pd.DataFrame:
        frame = self._read_input_table(
            path,
            context=context,
            encoding=encoding,
        )
        return self._validate_table(frame, validator=validator)

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
