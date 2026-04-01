from __future__ import annotations

from pathlib import Path

import pandas as pd

from .dataset_schema import DatasetSchema
from .io import read_table
from .validation.tables import PhosphoInputSchema, TotalInputSchema


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
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        validated_total = TotalInputSchema.validate(
            total_df,
            total_cols=self.schema.total_cols,
        )
        validated_phospho = PhosphoInputSchema.validate(
            phospho_df,
            phospho_cols=self.schema.phospho_cols,
        )
        return validated_total, validated_phospho

    def load(
        self,
        total_path: str | Path,
        phospho_path: str | Path,
        *,
        phospho_encoding: str | None = None,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        total_df = read_table(total_path)
        phospho_df = read_table(phospho_path, encoding=phospho_encoding)
        return self.validate(total_df=total_df, phospho_df=phospho_df)
