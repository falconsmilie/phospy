from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .dataset_schema import DatasetSchema
from .io import read_table
from .validation.tables import PhosphoInputSchema, TotalInputSchema


@dataclass(frozen=True, slots=True, init=False)
class ValidatedCoreInputs:
    """Validated dataset tables produced by :class:`DatasetLoader`.

    The loader owns construction of this type so downstream code can distinguish
    between arbitrary in-memory frames and frames that have passed the dataset
    boundary validators. Defensive copies are stored and returned to preserve the
    validated snapshot.
    """

    schema: DatasetSchema
    _total_df: pd.DataFrame
    _phospho_df: pd.DataFrame

    def __init__(
        self,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        schema: DatasetSchema,
    ) -> None:
        object.__setattr__(self, "schema", schema)
        object.__setattr__(self, "_total_df", total_df.copy(deep=True))
        object.__setattr__(self, "_phospho_df", phospho_df.copy(deep=True))

    @property
    def total_df(self) -> pd.DataFrame:
        return self._total_df.copy(deep=True)

    @property
    def phospho_df(self) -> pd.DataFrame:
        return self._phospho_df.copy(deep=True)


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
        total_df = read_table(total_path)
        phospho_df = read_table(phospho_path, encoding=phospho_encoding)
        return self.validate(total_df=total_df, phospho_df=phospho_df)
