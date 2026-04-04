from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .constants import ComparisonSpec
from .dataset_loader import DatasetLoader, ValidatedCoreInputs
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
    mutable pandas tables, and exposes those owned tables directly through
    `total_df` and `phospho_df`.

    Use `copy_inputs()` when you need detached caller-owned copies instead of the
    dataset's shared workspace state.
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
        """Validate raw inputs and take ownership of mutable workspace tables."""
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

    @property
    def total_df(self) -> pd.DataFrame:
        """Return the owned validated total-protein workspace table.

        Mutating the returned frame mutates this dataset's owned workspace state.
        """
        return self.inputs.total_df

    @property
    def phospho_df(self) -> pd.DataFrame:
        """Return the owned validated phosphoproteomics workspace table.

        Mutating the returned frame mutates this dataset's owned workspace state.
        """
        return self.inputs.phospho_df

    def copy_inputs(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return detached deep copies suitable for caller-owned mutation."""
        return self.inputs.copy_frames()

    @property
    def preprocessing(self) -> DatasetPreprocessing:
        """Return the bound preprocessing facade for this dataset workspace."""
        return DatasetPreprocessing(
            total_df=self.inputs.total_df,
            phospho_df=self.inputs.phospho_df,
            schema=self.schema,
            comparisons=self.comparisons,
        )

    @property
    def site_matrix(self) -> DatasetSiteMatrix:
        """Return the bound site-matrix facade for this dataset workspace."""
        return DatasetSiteMatrix(schema=self.schema)

    @classmethod
    def from_validated_inputs(
        cls,
        validated_inputs: ValidatedCoreInputs,
        *,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> PhosphoDataset:
        """Build a dataset workspace from trusted validated inputs."""
        if not isinstance(validated_inputs, ValidatedCoreInputs):
            msg = (
                "from_validated_inputs requires a ValidatedCoreInputs instance. "
                "Call DatasetLoader.validate(...) or DatasetLoader.load(...) first."
            )
            raise TypeError(msg)
        validated_request = build_validated_dataset_inputs(
            schema=validated_inputs.schema,
            total_df=validated_inputs.total_df,
            phospho_df=validated_inputs.phospho_df,
            comparisons=comparisons,
            context="PhosphoDataset",
        )
        instance = cls.__new__(cls)
        instance._set_state(validated_request=validated_request)
        return instance

    @classmethod
    def from_files(
        cls,
        total_path: str | Path,
        phospho_path: str | Path,
        phospho_encoding: str | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
        schema: DatasetSchema | None = None,
    ) -> PhosphoDataset:
        """Load validated files and return a dataset workspace owning their tables."""
        loader = DatasetLoader(schema=schema)
        validated_inputs = loader.load(
            total_path,
            phospho_path,
            phospho_encoding=phospho_encoding,
        )
        return cls.from_validated_inputs(
            validated_inputs,
            comparisons=comparisons,
        )
