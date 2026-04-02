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


@dataclass(frozen=True, slots=True)
class CoreInputs:
    """Owned dataset tables used by the core preprocessing pipeline.

    The frames stored here are the canonical validated in-memory dataset owned by
    :class:`PhosphoDataset`. They are exposed directly for read access. Call
    :meth:`copy_frames` before mutating them outside the processing boundary.
    """

    total_df: pd.DataFrame
    phospho_df: pd.DataFrame

    def copy_frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return deep copies suitable for caller-owned mutation."""
        return self.total_df.copy(deep=True), self.phospho_df.copy(deep=True)


@dataclass(frozen=True, slots=True, init=False)
class PhosphoDataset:
    """Explicit owner around validated phosphoproteomics inputs.

    ``PhosphoDataset`` keeps one trusted validated snapshot of the source frames.
    Public accessors return the owned in-memory frames directly rather than
    creating defensive copies on every read. Internal preprocessing remains safe
    because mutation happens on caller-owned copies inside the processing
    services. Call :meth:`copy_inputs` when you need a mutable copy of the
    dataset inputs.
    """

    inputs: CoreInputs
    schema: DatasetSchema
    comparisons: tuple[ComparisonSpec, ...] | None

    def __init__(
        self,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        *,
        schema: DatasetSchema | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> None:
        resolved_schema = schema or DatasetSchema()
        loader = DatasetLoader(schema=resolved_schema)
        validated_inputs = loader.validate(
            total_df=total_df,
            phospho_df=phospho_df,
        )
        validated_comparisons = validated_inputs.schema.validate_comparisons(
            comparisons,
            context="PhosphoDataset",
        )
        self._set_state(
            inputs=CoreInputs(
                total_df=validated_inputs.total_df,
                phospho_df=validated_inputs.phospho_df,
            ),
            schema=validated_inputs.schema,
            comparisons=validated_comparisons,
        )

    def _set_state(
        self,
        *,
        inputs: CoreInputs,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None,
    ) -> None:
        object.__setattr__(self, "inputs", inputs)
        object.__setattr__(self, "schema", schema)
        object.__setattr__(
            self,
            "comparisons",
            tuple(comparisons) if comparisons is not None else None,
        )

    @property
    def total_df(self) -> pd.DataFrame:
        """Return the owned validated total-protein table."""
        return self.inputs.total_df

    @property
    def phospho_df(self) -> pd.DataFrame:
        """Return the owned validated phosphoproteomics table."""
        return self.inputs.phospho_df

    def copy_inputs(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return deep copies suitable for caller-owned mutation."""
        return self.inputs.copy_frames()

    @property
    def preprocessing(self) -> DatasetPreprocessing:
        """Return the bound preprocessing facade for this dataset."""
        return DatasetPreprocessing(
            total_df=self.inputs.total_df,
            phospho_df=self.inputs.phospho_df,
            schema=self.schema,
            comparisons=self.comparisons,
        )

    @property
    def site_matrix(self) -> DatasetSiteMatrix:
        """Return the bound site-matrix facade for this dataset."""
        return DatasetSiteMatrix(schema=self.schema)

    @classmethod
    def from_validated_inputs(
        cls,
        validated_inputs: ValidatedCoreInputs,
        *,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> PhosphoDataset:
        """Build a dataset from validated inputs produced by ``DatasetLoader``."""
        if not isinstance(validated_inputs, ValidatedCoreInputs):
            msg = (
                "validated_inputs must be a ValidatedCoreInputs instance "
                "produced by DatasetLoader"
            )
            raise TypeError(msg)
        validated_comparisons = validated_inputs.schema.validate_comparisons(
            comparisons,
            context="PhosphoDataset",
        )
        instance = cls.__new__(cls)
        instance._set_state(
            inputs=CoreInputs(
                total_df=validated_inputs.total_df,
                phospho_df=validated_inputs.phospho_df,
            ),
            schema=validated_inputs.schema,
            comparisons=validated_comparisons,
        )
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
