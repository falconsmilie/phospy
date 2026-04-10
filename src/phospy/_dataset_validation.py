from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .constants import ComparisonSpec
from .dataset_schema import DatasetSchema
from .validation.errors import InputCompatibilityError
from .validation.schema.files import validate_existing_file_path
from .validation.schema.tables import PhosphoInputSchema, TotalInputSchema


@dataclass(frozen=True, slots=True)
class _ValidatedDatasetPaths:
    """Validated internal file-backed dataset input paths."""

    total_path: Path
    phospho_path: Path


@dataclass(slots=True)
class _ValidatedDatasetInputs:
    """Trusted internal dataset bundle for dataset construction boundaries.

    This internal helper carries validated pandas tables into dataset creation.
    It is trusted by convention, not structurally immutable.
    """

    schema: DatasetSchema
    total_df: pd.DataFrame
    phospho_df: pd.DataFrame
    comparisons: tuple[ComparisonSpec, ...] | None = None


def _validate_dataset_file_paths(
    total_path: str | Path,
    phospho_path: str | Path,
) -> _ValidatedDatasetPaths:
    """Validate dataset file paths before table loading."""

    return _ValidatedDatasetPaths(
        total_path=validate_existing_file_path(
            total_path,
            context="total input table path",
        ),
        phospho_path=validate_existing_file_path(
            phospho_path,
            context="phospho input table path",
        ),
    )


def _validate_dataset_frames(
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


def _validate_dataset_request(
    *,
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    schema: DatasetSchema | None = None,
    comparisons: Sequence[ComparisonSpec] | None = None,
    context: str = "PhosphoDataset",
) -> _ValidatedDatasetInputs:
    """Validate raw dataset inputs for internal dataset boundaries."""

    resolved_schema = schema or DatasetSchema()
    validated_total, validated_phospho = _validate_dataset_frames(
        total_df=total_df,
        phospho_df=phospho_df,
        schema=resolved_schema,
    )
    validated_comparisons = _validate_dataset_comparisons(
        schema=resolved_schema,
        comparisons=comparisons,
        context=context,
    )
    return _ValidatedDatasetInputs(
        schema=resolved_schema,
        total_df=validated_total,
        phospho_df=validated_phospho,
        comparisons=validated_comparisons,
    )


def _build_validated_dataset_inputs(
    *,
    schema: DatasetSchema,
    total_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    comparisons: Sequence[ComparisonSpec] | None = None,
    context: str = "PhosphoDataset",
) -> _ValidatedDatasetInputs:
    """Build an internal validated dataset request from trusted frames."""

    return _ValidatedDatasetInputs(
        schema=schema,
        total_df=total_df,
        phospho_df=phospho_df,
        comparisons=_validate_dataset_comparisons(
            schema=schema,
            comparisons=comparisons,
            context=context,
        ),
    )


def _validate_dataset_comparisons(
    *,
    schema: DatasetSchema,
    comparisons: Sequence[ComparisonSpec] | None,
    context: str,
) -> tuple[ComparisonSpec, ...] | None:
    try:
        return schema.validate_comparisons(comparisons, context=context)
    except (InputCompatibilityError, TypeError, ValueError):
        raise
