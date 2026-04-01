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


@dataclass(frozen=True, slots=True, init=False)
class CoreInputs:
    """Validated in-memory tables used by the core preprocessing pipeline.

    Defensive copies are taken on construction and when exposing the stored
    tables so callers cannot mutate dataset state through the public API.
    """

    _total_df: pd.DataFrame
    _phospho_df: pd.DataFrame

    def __init__(self, total_df: pd.DataFrame, phospho_df: pd.DataFrame) -> None:
        object.__setattr__(self, "_total_df", total_df.copy(deep=True))
        object.__setattr__(self, "_phospho_df", phospho_df.copy(deep=True))

    @property
    def total_df(self) -> pd.DataFrame:
        return self._total_df.copy(deep=True)

    @property
    def phospho_df(self) -> pd.DataFrame:
        return self._phospho_df.copy(deep=True)


@dataclass(frozen=True, slots=True, init=False)
class PhosphoDataset:
    """Immutable holder around validated phosphoproteomics inputs.

    The dataset stores defensive copies of validated tables and returns deep
    copies from its public accessors so callers cannot mutate internal state.
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
        self._set_state(
            inputs=CoreInputs(
                total_df=validated_inputs.total_df,
                phospho_df=validated_inputs.phospho_df,
            ),
            schema=validated_inputs.schema,
            comparisons=comparisons,
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
        return self.inputs.total_df

    @property
    def phospho_df(self) -> pd.DataFrame:
        return self.inputs.phospho_df

    @property
    def preprocessing(self) -> DatasetPreprocessing:
        """Return the bound preprocessing facade for this dataset."""
        return DatasetPreprocessing(
            total_df=self.total_df,
            phospho_df=self.phospho_df,
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
        instance = cls.__new__(cls)
        instance._set_state(
            inputs=CoreInputs(
                total_df=validated_inputs.total_df,
                phospho_df=validated_inputs.phospho_df,
            ),
            schema=validated_inputs.schema,
            comparisons=comparisons,
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
