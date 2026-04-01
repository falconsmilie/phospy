from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .constants import ComparisonSpec
from .core_processing import CoreProcessingResult
from .dataset_loader import DatasetLoader
from .dataset_preprocessing import DatasetPreprocessing
from .dataset_schema import DatasetSchema
from .dataset_site_matrix import DatasetSiteMatrix
from .writers import CoreOutputWriter, CoreProcessingResultWriter


@dataclass(frozen=True, slots=True)
class CoreInputs:
    """Validated in-memory tables used by the core preprocessing pipeline."""

    total_df: pd.DataFrame
    phospho_df: pd.DataFrame


@dataclass(frozen=True, slots=True, init=False)
class PhosphoDataset:
    """Thin immutable holder around validated phosphoproteomics inputs."""

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
        validated_total, validated_phospho = loader.validate(
            total_df=total_df,
            phospho_df=phospho_df,
        )
        self._set_state(
            inputs=CoreInputs(total_df=validated_total, phospho_df=validated_phospho),
            schema=resolved_schema,
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
        *,
        total_df: pd.DataFrame,
        phospho_df: pd.DataFrame,
        schema: DatasetSchema | None = None,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> PhosphoDataset:
        """Build a dataset from validated in-memory inputs."""
        instance = cls.__new__(cls)
        instance._set_state(
            inputs=CoreInputs(total_df=total_df, phospho_df=phospho_df),
            schema=schema or DatasetSchema(),
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
        total_df, phospho_df = loader.load(
            total_path,
            phospho_path,
            phospho_encoding=phospho_encoding,
        )
        return cls.from_validated_inputs(
            total_df=total_df,
            phospho_df=phospho_df,
            comparisons=comparisons,
            schema=loader.schema,
        )

    @staticmethod
    def write_core_outputs(
        result: CoreProcessingResult,
        outdir: str | Path,
        *,
        writer: CoreProcessingResultWriter = CoreOutputWriter,
    ) -> None:
        writer.write(result, outdir)
