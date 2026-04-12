from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ...datasets.schema import DatasetSchema
from ...internal.constants import ComparisonSpec
from ..domain import validate_dataset_comparisons
from ..schema.files import validate_existing_file_path
from ..schema.tables import PhosphoInputSchema, TotalInputSchema


@dataclass(frozen=True, slots=True)
class DatasetFilePaths:
    """Validated file-backed dataset input paths."""

    total_path: Path
    phospho_path: Path


@dataclass(slots=True)
class DatasetInputs:
    """Trusted dataset inputs owned by the dataset boundary."""

    schema: DatasetSchema
    total_df: pd.DataFrame
    phospho_df: pd.DataFrame
    comparisons: tuple[ComparisonSpec, ...] | None = None


def validate_dataset_file_paths(
    total_path: str | Path,
    phospho_path: str | Path,
) -> DatasetFilePaths:
    """Validate dataset file paths before table loading."""

    return DatasetFilePaths(
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


def validate_dataset_request(
    *,
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[ComparisonSpec] | None = None,
    context: str = "PhosphoDataset",
) -> DatasetInputs:
    """Validate raw dataset inputs for the public dataset boundary."""

    resolved_schema = schema or DatasetSchema()
    validated_total, validated_phospho = validate_dataset_frames(
        total_df=total_df,
        phospho_df=phospho_df,
        schema=resolved_schema,
    )
    validated_comparisons = validate_dataset_comparisons(
        schema=resolved_schema,
        comparisons=comparisons,
        context=context,
    )
    return DatasetInputs(
        schema=resolved_schema,
        total_df=validated_total,
        phospho_df=validated_phospho,
        comparisons=validated_comparisons,
    )


def build_dataset_inputs(
    *,
    schema: DatasetSchema,
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    comparisons: Sequence[ComparisonSpec] | None = None,
    context: str = "PhosphoDataset",
) -> DatasetInputs:
    """Build dataset-owned inputs from already validated tables."""

    return DatasetInputs(
        schema=schema,
        total_df=total_df,
        phospho_df=phospho_df,
        comparisons=validate_dataset_comparisons(
            schema=schema,
            comparisons=comparisons,
            context=context,
        ),
    )
