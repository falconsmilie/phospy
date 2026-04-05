from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .constants import ComparisonSpec
from .dataset_loader import DatasetLoader, LoadedDatasetInputs
from .dataset_preprocessing import DatasetPreprocessing
from .dataset_schema import DatasetSchema
from .dataset_site_matrix import DatasetSiteMatrix
from .validation.dataset import (
    ValidatedDatasetInputs,
    build_validated_dataset_inputs,
    validate_dataset_request,
)


@dataclass(slots=True)
class CoreInputs:
    """Owned mutable dataset tables used by the core preprocessing pipeline."""

    total_df: pd.DataFrame
    phospho_df: pd.DataFrame

    def copy_frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return deep copies suitable for caller-owned mutation."""
        return self.total_df.copy(deep=True), self.phospho_df.copy(deep=True)


class PhosphoDataset:
    """Mutable workspace that owns validated phosphoproteomics input tables.

    `PhosphoDataset` is an owned in-memory processing workspace, not an immutable
    snapshot. It validates raw constructor inputs at the boundary, stores owned
    mutable pandas tables, and exposes them through explicit `*_live` and
    `*_copy` accessors.

    Prefer `total_df_copy` / `phospho_df_copy` or `copy_inputs()` for caller-owned
    inspection, export, and other read-oriented work. Use `total_df_live` /
    `phospho_df_live` only when you intentionally want shared workspace state.
    """

    __slots__ = ("_inputs", "_schema", "_comparisons")

    def __init__(
        self,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        *,
        schema: DatasetSchema | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> None:
        """Validate raw inputs and take ownership of isolated mutable workspace tables.

        Raw caller-supplied frames are isolated at this boundary so later mutation of
        the caller's original inputs does not affect the dataset workspace.
        """
        validated_request = validate_dataset_request(
            total_df=total_df,
            phospho_df=phospho_df,
            schema=schema,
            comparisons=comparisons,
        )
        self._set_state(validated_request=validated_request)

    @property
    def inputs(self) -> CoreInputs:
        """Return the owned mutable input bundle for advanced internal workflows."""
        return self._inputs

    @property
    def schema(self) -> DatasetSchema:
        """Return the schema governing the owned workspace tables."""
        return self._schema

    @property
    def comparisons(self) -> tuple[ComparisonSpec, ...] | None:
        """Return the owned comparison specs bound to this dataset workspace."""
        return self._comparisons

    def _set_state(self, *, validated_request: ValidatedDatasetInputs) -> None:
        self._inputs = CoreInputs(
            total_df=validated_request.total_df,
            phospho_df=validated_request.phospho_df,
        )
        self._schema = validated_request.schema
        self._comparisons = validated_request.comparisons

    @classmethod
    def _from_owned_validated_request(
        cls,
        validated_request: ValidatedDatasetInputs,
    ) -> PhosphoDataset:
        """Build a dataset from already-owned validated frames without copying again."""
        instance = cls.__new__(cls)
        instance._set_state(validated_request=validated_request)
        return instance

    @property
    def total_df_live(self) -> pd.DataFrame:
        """Return the owned validated total-protein workspace table.

        Mutating the returned frame mutates this dataset's owned workspace state.
        Prefer `total_df_copy` for read-oriented caller work.
        """
        return self.inputs.total_df

    @property
    def phospho_df_live(self) -> pd.DataFrame:
        """Return the owned validated phosphoproteomics workspace table.

        Mutating the returned frame mutates this dataset's owned workspace state.
        Prefer `phospho_df_copy` for read-oriented caller work.
        """
        return self.inputs.phospho_df

    @property
    def total_df_copy(self) -> pd.DataFrame:
        """Return a detached deep copy of the validated total-protein table."""
        return self.inputs.total_df.copy(deep=True)

    @property
    def phospho_df_copy(self) -> pd.DataFrame:
        """Return a detached deep copy of the validated phosphoproteomics table."""
        return self.inputs.phospho_df.copy(deep=True)

    def copy_inputs(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return detached deep copies suitable for caller-owned mutation and read use."""
        return self.total_df_copy, self.phospho_df_copy

    @property
    def preprocessing(self) -> DatasetPreprocessing:
        """Return the bound preprocessing facade for this dataset workspace."""
        return DatasetPreprocessing(
            total_df=self.total_df_live,
            phospho_df=self.phospho_df_live,
            schema=self.schema,
            comparisons=self.comparisons,
        )

    @property
    def site_matrix(self) -> DatasetSiteMatrix:
        """Return the bound site-matrix facade for this dataset workspace."""
        return DatasetSiteMatrix(schema=self.schema)

    @classmethod
    def _from_loaded_inputs(
        cls,
        loaded_inputs: LoadedDatasetInputs,
        *,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> PhosphoDataset:
        """Build a dataset workspace from internal trusted loader output.

        This internal boundary transfers ownership of the loader-managed frames
        directly into the dataset workspace. Public callers should use
        ``PhosphoDataset(...)`` or ``PhosphoDataset.from_files(...)``.
        """
        if not isinstance(loaded_inputs, LoadedDatasetInputs):
            msg = (
                "_from_loaded_inputs requires an internal LoadedDatasetInputs instance."
            )
            raise TypeError(msg)
        validated_request = build_validated_dataset_inputs(
            schema=loaded_inputs.schema,
            total_df=loaded_inputs.total_df,
            phospho_df=loaded_inputs.phospho_df,
            comparisons=comparisons,
            context="PhosphoDataset",
        )
        return cls._from_owned_validated_request(validated_request)

    @classmethod
    def from_files(
        cls,
        total_path: str | Path,
        phospho_path: str | Path,
        phospho_encoding: str | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
        schema: DatasetSchema | None = None,
    ) -> PhosphoDataset:
        """Load files through the public dataset boundary and return an owned workspace."""
        loader = DatasetLoader(schema=schema)
        validated_inputs = loader.load(
            total_path,
            phospho_path,
            phospho_encoding=phospho_encoding,
        )
        return cls._from_loaded_inputs(
            validated_inputs,
            comparisons=comparisons,
        )
