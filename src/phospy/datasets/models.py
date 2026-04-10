from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..constants import (
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
    ComparisonSpec,
)
from ..core_processing import CorePreprocessingConfig, CoreProcessingResult
from ..dataset_preprocessing import DatasetPreprocessing
from .builders import (
    DatasetSiteMatrix,
    _build_site_metadata,
    _validate_analysis_ready_alignment,
)
from .schema import DatasetSchema

if TYPE_CHECKING:
    from ..validation.requests import ValidatedDatasetInputs
    from .loaders import LoadedDatasetInputs

__all__ = [
    "AnalysisReadyPhosphoDataset",
    "AnalysisReadyPreprocessingProvenance",
    "AnalysisReadyRowCounts",
    "AnalysisReadySiteMatrixStats",
    "CoreInputs",
    "PhosphoDataset",
]


@dataclass(slots=True)
class CoreInputs:
    """Owned mutable dataset tables used by the core preprocessing pipeline."""

    total_df: pd.DataFrame
    phospho_df: pd.DataFrame

    def copy_frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return deep copies suitable for caller-owned mutation."""
        return self.total_df.copy(deep=True), self.phospho_df.copy(deep=True)


@dataclass(frozen=True, slots=True)
class AnalysisReadyRowCounts:
    """Row counts captured across the core preprocessing stages."""

    total_unique: int
    total_filtered: int
    phospho_filtered: int
    phospho_corrected: int
    phospho_matrix_sites: int


@dataclass(frozen=True, slots=True)
class AnalysisReadySiteMatrixStats:
    """Typed row-drop diagnostics from site-matrix construction."""

    input_rows: int
    dropped_missing_sequence: int
    dropped_incomplete_values: int
    deduplicated_site_rows: int
    retained_rows: int

    @classmethod
    def from_mapping(
        cls,
        row_drop_stats: dict[str, int],
    ) -> AnalysisReadySiteMatrixStats:
        return cls(
            input_rows=int(row_drop_stats.get("input_rows", 0)),
            dropped_missing_sequence=int(
                row_drop_stats.get("dropped_missing_sequence", 0)
            ),
            dropped_incomplete_values=int(
                row_drop_stats.get("dropped_incomplete_values", 0)
            ),
            deduplicated_site_rows=int(row_drop_stats.get("deduplicated_site_rows", 0)),
            retained_rows=int(row_drop_stats.get("retained_rows", 0)),
        )


@dataclass(frozen=True, slots=True)
class AnalysisReadyPreprocessingProvenance:
    """Preprocessing provenance for one analysis-ready phosphosite dataset."""

    source: str
    schema: DatasetSchema
    comparisons: tuple[ComparisonSpec, ...] | None
    row_counts: AnalysisReadyRowCounts
    site_matrix_stats: AnalysisReadySiteMatrixStats


@dataclass(frozen=True, slots=True, init=False)
class AnalysisReadyPhosphoDataset:
    """Owned analysis-ready phosphosite state between preprocessing and inference.

    This immutable boundary object carries the minimum phosphosite state needed
    after preprocessing and before workflow-specific kinase inputs are resolved.
    It intentionally separates:

    - the phosphosite analysis matrix used for inference
    - aligned site metadata keyed by site identifier
    - aligned site-centred sequences
    - the corrected phosphosite source table the matrix was derived from
    - preprocessing provenance describing how the boundary was produced
    """

    phospho_matrix: pd.DataFrame
    site_metadata: pd.DataFrame
    site_sequences: pd.Series
    phospho_corrected: pd.DataFrame
    provenance: AnalysisReadyPreprocessingProvenance

    def __init__(
        self,
        *,
        phospho_matrix: pd.DataFrame,
        site_metadata: pd.DataFrame,
        site_sequences: pd.Series,
        phospho_corrected: pd.DataFrame,
        provenance: AnalysisReadyPreprocessingProvenance,
    ) -> None:
        owned_phospho_matrix = phospho_matrix.copy(deep=True)
        owned_site_metadata = site_metadata.copy(deep=True)
        owned_site_sequences = site_sequences.copy(deep=True)
        owned_phospho_corrected = phospho_corrected.copy(deep=True)

        _validate_analysis_ready_alignment(
            phospho_matrix=owned_phospho_matrix,
            site_metadata=owned_site_metadata,
            site_sequences=owned_site_sequences,
        )

        object.__setattr__(self, "phospho_matrix", owned_phospho_matrix)
        object.__setattr__(self, "site_metadata", owned_site_metadata)
        object.__setattr__(self, "site_sequences", owned_site_sequences)
        object.__setattr__(self, "phospho_corrected", owned_phospho_corrected)
        object.__setattr__(self, "provenance", provenance)

    @classmethod
    def from_core_processing_result(
        cls,
        result: CoreProcessingResult,
        *,
        schema: DatasetSchema,
        comparisons: Sequence[ComparisonSpec] | None = None,
        source: str = "core preprocessing",
    ) -> AnalysisReadyPhosphoDataset:
        """Build an analysis-ready dataset from existing preprocessing output."""
        if not isinstance(result, CoreProcessingResult):
            msg = (
                "AnalysisReadyPhosphoDataset.from_core_processing_result() "
                "requires a CoreProcessingResult instance."
            )
            raise TypeError(msg)

        site_metadata = _build_site_metadata(
            phosr_input=result.site_matrix.phosr_input,
            corrected_cols=schema.corrected_cols,
        )
        row_counts = AnalysisReadyRowCounts(
            total_unique=len(result.total_unique),
            total_filtered=len(result.total_filtered),
            phospho_filtered=len(result.phospho_filtered),
            phospho_corrected=len(result.phospho_corrected),
            phospho_matrix_sites=len(result.site_matrix.matrix),
        )
        provenance = AnalysisReadyPreprocessingProvenance(
            source=str(source),
            schema=schema,
            comparisons=schema.validate_comparisons(
                comparisons,
                context="AnalysisReadyPhosphoDataset provenance comparisons",
            ),
            row_counts=row_counts,
            site_matrix_stats=AnalysisReadySiteMatrixStats.from_mapping(
                result.site_matrix.row_drop_stats
            ),
        )
        return cls(
            phospho_matrix=result.site_matrix.matrix,
            site_metadata=site_metadata,
            site_sequences=result.site_matrix.sequences,
            phospho_corrected=result.phospho_corrected,
            provenance=provenance,
        )


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
        from ..validation.requests import validate_dataset_request

        validated_request = validate_dataset_request(
            total_df=total_df,
            phospho_df=phospho_df,
            schema=schema,
            comparisons=comparisons,
            context="PhosphoDataset",
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

    def to_analysis_ready(
        self,
        result: CoreProcessingResult,
        *,
        source: str = "dataset preprocessing",
    ) -> AnalysisReadyPhosphoDataset:
        """Adapt a preprocessing result through the supported dataset-bound path."""
        return self.preprocessing.to_analysis_ready(result, source=source)

    def run_analysis_ready(
        self,
        localization_threshold: float = 0.75,
        min_observed: int = 4,
        max_unmatched_fraction: float = 0.0,
        total_sentinel: float | int = DEFAULT_TOTAL_SENTINEL,
        phospho_sentinel: float | int = DEFAULT_PHOSPHO_SENTINEL,
        config: CorePreprocessingConfig | None = None,
        source: str = "dataset preprocessing",
    ) -> AnalysisReadyPhosphoDataset:
        """Run preprocessing and return an analysis-ready phosphosite dataset."""
        return self.preprocessing.run_analysis_ready(
            localization_threshold=localization_threshold,
            min_observed=min_observed,
            max_unmatched_fraction=max_unmatched_fraction,
            total_sentinel=total_sentinel,
            phospho_sentinel=phospho_sentinel,
            config=config,
            source=source,
        )

    @classmethod
    def from_loaded_inputs(
        cls,
        loaded_inputs: LoadedDatasetInputs,
        *,
        comparisons: Sequence[ComparisonSpec] | None = None,
    ) -> PhosphoDataset:
        """Build a dataset workspace from internal trusted loader output.

        This explicit internal boundary transfers ownership of the loader-managed
        frames directly into the dataset workspace. Public callers should use
        ``PhosphoDataset(...)`` or ``PhosphoDataset.from_files(...)``.
        """
        from ..validation.requests import build_validated_dataset_inputs
        from .loaders import LoadedDatasetInputs

        if not isinstance(loaded_inputs, LoadedDatasetInputs):
            msg = (
                "from_loaded_inputs requires an internal LoadedDatasetInputs instance."
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
        from .loaders import DatasetLoader

        loader = DatasetLoader(schema=schema)
        validated_inputs = loader.load(
            total_path,
            phospho_path,
            phospho_encoding=phospho_encoding,
        )
        return cls.from_loaded_inputs(
            validated_inputs,
            comparisons=comparisons,
        )
