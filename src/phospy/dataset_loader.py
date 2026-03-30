from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from .constants import DEFAULT_PHOSPHO_COLS, DEFAULT_TOTAL_COLS
from .io import read_table
from .validation.tables import PhosphoInputSchema, TotalInputSchema


class DatasetLoader:
    """Load and validate dataset frames from memory or disk."""

    def __init__(
        self,
        *,
        total_cols: Sequence[str] | None = None,
        phospho_cols: Sequence[str] | None = None,
    ) -> None:
        self.total_cols = list(total_cols or DEFAULT_TOTAL_COLS)
        self.phospho_cols = list(phospho_cols or DEFAULT_PHOSPHO_COLS)

    def validate(
        self,
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        validated_total = TotalInputSchema.validate(
            total_df, total_cols=self.total_cols
        )
        validated_phospho = PhosphoInputSchema.validate(
            phospho_df,
            phospho_cols=self.phospho_cols,
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
