from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..dataset_schema import DatasetSchema
from .paths import validate_existing_file_path
from .tables import PhosphoInputSchema, TotalInputSchema


@dataclass(frozen=True, slots=True)
class ValidatedDatasetPaths:
    """Validated file-backed dataset input paths."""

    total_path: Path
    phospho_path: Path


def validate_dataset_file_paths(
    total_path: str | Path,
    phospho_path: str | Path,
) -> ValidatedDatasetPaths:
    """Validate dataset file paths before table loading."""

    return ValidatedDatasetPaths(
        total_path=validate_existing_file_path(
            total_path,
            context="total input table path",
        ),
        phospho_path=validate_existing_file_path(
            phospho_path,
            context="phospho input table path",
        ),
    )


def validate_dataset_frames(
    *,
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    schema: DatasetSchema,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate in-memory dataset tables against the configured schema."""

    validated_total = TotalInputSchema.validate(
        total_df,
        total_cols=schema.total_cols,
    )
    validated_phospho = PhosphoInputSchema.validate(
        phospho_df,
        phospho_cols=schema.phospho_cols,
    )
    return validated_total, validated_phospho


__all__ = [
    "ValidatedDatasetPaths",
    "validate_dataset_file_paths",
    "validate_dataset_frames",
]
