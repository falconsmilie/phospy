from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..errors import InputCompatibilityError
from ..internal.constants import (
    PHOSPHO_GENE_COLUMN,
    SITE_MATRIX_GENE_COLUMN,
    SITE_MATRIX_ID_COLUMN,
    ComparisonSpec,
)
from ..preprocessing.core import CorePreprocessingConfig, CoreProcessingResult
from ..preprocessing.dataset import DatasetPreprocessing
from .builders import (
    DatasetSiteMatrix,
    _build_site_metadata,
    _validate_analysis_ready_alignment,
)
from .schema import DatasetSchema

if TYPE_CHECKING:
    from ..validation.requests.dataset import DatasetInputs
    from .loaders import LoadedDatasetInputs

__all__ = [
    "AnalysisReadyPhosphoDataset",
    "AnalysisReadyPreprocessingProvenance",
    "AnalysisReadyRowCounts",
    "AnalysisReadySiteMatrixStats",
    "CoreInputs",
    "PhosphoDataset",
]

_DEFAULT_SITE_TO_PROTEIN_METADATA_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "protein",
    SITE_MATRIX_GENE_COLUMN,
    PHOSPHO_GENE_COLUMN,
)


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
    duplicate_site_strategy: str = "max_mean_signal"

    @classmethod
    def from_mapping(
        cls,
        row_drop_stats: Mapping[str, int | str],
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
            duplicate_site_strategy=str(
                row_drop_stats.get("duplicate_site_strategy", "max_mean_signal")
            ),
        )


@dataclass(frozen=True, slots=True)
class AnalysisReadyPreprocessingProvenance:
    """Preprocessing provenance for one analysis-ready phosphosite dataset."""

    source: str
    schema: DatasetSchema
    comparisons: tuple[ComparisonSpec, ...] | None
    row_counts: AnalysisReadyRowCounts
    site_matrix_stats: AnalysisReadySiteMatrixStats


@dataclass(slots=True, init=False)
class AnalysisReadyPhosphoDataset:
    """Owned analysis-ready phosphosite state between preprocessing and inference.

    This owned boundary object carries the minimum phosphosite state needed
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
        """Create an owned analysis-ready dataset from caller-supplied inputs.

        The constructor is the external ownership boundary. It defensively
        deep-copies caller-managed pandas inputs once, then delegates to the
        owned fast path used by internal preprocessing results.
        """

        owned = self.from_external(
            phospho_matrix=phospho_matrix,
            site_metadata=site_metadata,
            site_sequences=site_sequences,
            phospho_corrected=phospho_corrected,
            provenance=provenance,
        )
        self.phospho_matrix = owned.phospho_matrix
        self.site_metadata = owned.site_metadata
        self.site_sequences = owned.site_sequences
        self.phospho_corrected = owned.phospho_corrected
        self.provenance = owned.provenance

    @classmethod
    def from_external(
        cls,
        *,
        phospho_matrix: pd.DataFrame,
        site_metadata: pd.DataFrame,
        site_sequences: pd.Series,
        phospho_corrected: pd.DataFrame,
        provenance: AnalysisReadyPreprocessingProvenance,
    ) -> AnalysisReadyPhosphoDataset:
        """Create an analysis-ready dataset by taking ownership of external inputs.

        This boundary isolates stored state from later caller mutation by copying
        each pandas object once.
        """
        return cls.from_owned(
            phospho_matrix=phospho_matrix.copy(deep=True),
            site_metadata=site_metadata.copy(deep=True),
            site_sequences=site_sequences.copy(deep=True),
            phospho_corrected=phospho_corrected.copy(deep=True),
            provenance=provenance,
        )

    @classmethod
    def from_owned(
        cls,
        *,
        phospho_matrix: pd.DataFrame,
        site_metadata: pd.DataFrame,
        site_sequences: pd.Series,
        phospho_corrected: pd.DataFrame,
        provenance: AnalysisReadyPreprocessingProvenance,
    ) -> AnalysisReadyPhosphoDataset:
        """Create an analysis-ready dataset from already-owned aligned tables."""
        _validate_analysis_ready_alignment(
            phospho_matrix=phospho_matrix,
            site_metadata=site_metadata,
            site_sequences=site_sequences,
        )
        instance = cls.__new__(cls)
        instance.phospho_matrix = phospho_matrix
        instance.site_metadata = site_metadata
        instance.site_sequences = site_sequences
        instance.phospho_corrected = phospho_corrected
        instance.provenance = provenance
        return instance

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
        return cls.from_owned(
            phospho_matrix=result.site_matrix.matrix,
            site_metadata=site_metadata,
            site_sequences=result.site_matrix.sequences,
            phospho_corrected=result.phospho_corrected,
            provenance=provenance,
        )

    def resolve_site_to_protein_mapping(
        self,
        *,
        metadata_columns: Sequence[str] | None = None,
    ) -> pd.Series:
        """Resolve an aligned site-to-protein mapping from site metadata."""

        resolved_candidates = (
            _DEFAULT_SITE_TO_PROTEIN_METADATA_COLUMNS
            if metadata_columns is None
            else tuple(metadata_columns)
        )
        candidate_columns = tuple(
            str(column).strip()
            for column in resolved_candidates
            if column is not None and str(column).strip()
        )
        if not candidate_columns:
            msg = (
                "metadata_columns must contain at least one non-empty column name "
                "when resolving site-to-protein mapping."
            )
            raise ValueError(msg)

        site_index = pd.Index(
            self.phospho_matrix.index.astype("string"),
            name=SITE_MATRIX_ID_COLUMN,
        )
        metadata_index = pd.Index(
            self.site_metadata.index.astype("string"),
            name=SITE_MATRIX_ID_COLUMN,
        )
        if not site_index.equals(metadata_index):
            msg = (
                "AnalysisReadyPhosphoDataset.site_metadata must remain aligned with "
                "phospho_matrix to resolve site-to-protein mapping."
            )
            raise InputCompatibilityError(msg)

        available_columns = {
            str(column): column for column in self.site_metadata.columns
        }
        checked_columns: list[str] = []
        incomplete_column_diagnostics: list[str] = []

        for candidate in candidate_columns:
            original_column = available_columns.get(candidate)
            if original_column is None:
                continue
            checked_columns.append(candidate)

            raw_values = self.site_metadata.loc[:, original_column]
            normalized_values = raw_values.astype("string")
            stripped_values = normalized_values.str.strip()
            invalid_mask = normalized_values.isna() | stripped_values.eq("")
            if bool(invalid_mask.any()):
                invalid_site_ids = metadata_index[invalid_mask].astype(str).tolist()
                preview = ", ".join(invalid_site_ids[:3])
                diagnostic = (
                    f"column '{candidate}' has missing/empty values for: {preview}"
                )
                if len(invalid_site_ids) > 3:
                    diagnostic += ", ..."
                incomplete_column_diagnostics.append(diagnostic)
                continue

            resolved_index = pd.Index(
                site_index.astype(str),
                dtype=object,
                name=SITE_MATRIX_ID_COLUMN,
            )
            return pd.Series(
                stripped_values.astype(object).tolist(),
                index=resolved_index,
                dtype=object,
                name="protein_id",
            )

        checked_preview = ", ".join(candidate_columns)
        if not checked_columns:
            available_preview = ", ".join(
                str(column) for column in self.site_metadata.columns[:5]
            )
            if len(self.site_metadata.columns) > 5:
                available_preview += ", ..."
            msg = (
                "AnalysisReadyPhosphoDataset.site_metadata does not include a usable "
                "site-to-protein column. "
                f"Checked columns: {checked_preview}. "
            )
            if available_preview:
                msg += f"Available columns: {available_preview}. "
            else:
                msg += "site_metadata has no columns. "
            msg += "Provide site_to_protein explicitly when running signalome analysis."
            raise InputCompatibilityError(msg)

        diagnostics = "; ".join(incomplete_column_diagnostics)
        msg = (
            "AnalysisReadyPhosphoDataset.site_metadata does not contain a complete "
            "site-to-protein mapping. "
            f"Checked columns: {checked_preview}. {diagnostics}. "
            "Provide site_to_protein explicitly when running signalome analysis."
        )
        raise InputCompatibilityError(msg)


class PhosphoDataset:
    """Mutable workspace that owns owned phosphoproteomics input tables.

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
        from ..validation.requests.dataset import validate_dataset_request

        dataset_inputs = validate_dataset_request(
            total_df=total_df,
            phospho_df=phospho_df,
            schema=schema,
            comparisons=comparisons,
            context="PhosphoDataset",
        )
        self._set_state(dataset_inputs=dataset_inputs)

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

    def _set_state(self, *, dataset_inputs: DatasetInputs) -> None:
        self._inputs = CoreInputs(
            total_df=dataset_inputs.total_df,
            phospho_df=dataset_inputs.phospho_df,
        )
        self._schema = dataset_inputs.schema
        self._comparisons = dataset_inputs.comparisons

    @classmethod
    def _from_owned_dataset_inputs(
        cls,
        dataset_inputs: DatasetInputs,
    ) -> PhosphoDataset:
        """Build a dataset from already-owned input tables without copying again."""
        instance = cls.__new__(cls)
        instance._set_state(dataset_inputs=dataset_inputs)
        return instance

    @property
    def total_df_live(self) -> pd.DataFrame:
        """Return the owned total-protein workspace table.

        Mutating the returned frame mutates this dataset's owned workspace state.
        Prefer `total_df_copy` for read-oriented caller work.
        """
        return self.inputs.total_df

    @property
    def phospho_df_live(self) -> pd.DataFrame:
        """Return the owned phosphoproteomics workspace table.

        Mutating the returned frame mutates this dataset's owned workspace state.
        Prefer `phospho_df_copy` for read-oriented caller work.
        """
        return self.inputs.phospho_df

    @property
    def total_df_copy(self) -> pd.DataFrame:
        """Return a detached deep copy of the owned total-protein table."""
        return self.inputs.total_df.copy(deep=True)

    @property
    def phospho_df_copy(self) -> pd.DataFrame:
        """Return a detached deep copy of the owned phosphoproteomics table."""
        return self.inputs.phospho_df.copy(deep=True)

    def copy_inputs(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return detached deep copies suitable for caller-owned mutation and read use."""
        return self.total_df_copy, self.phospho_df_copy

    @property
    def preprocessing(self) -> DatasetPreprocessing:
        """Return the bound preprocessing facade for this dataset workspace."""
        return DatasetPreprocessing.from_owned(
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
        *,
        config: CorePreprocessingConfig,
        source: str = "dataset preprocessing",
    ) -> AnalysisReadyPhosphoDataset:
        """Run preprocessing and return an analysis-ready phosphosite dataset."""
        return self.preprocessing.run_analysis_ready(
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
        from ..validation.requests.dataset import build_dataset_inputs
        from .loaders import LoadedDatasetInputs

        if not isinstance(loaded_inputs, LoadedDatasetInputs):
            msg = (
                "from_loaded_inputs requires an internal LoadedDatasetInputs instance."
            )
            raise TypeError(msg)
        dataset_inputs = build_dataset_inputs(
            schema=loaded_inputs.schema,
            total_df=loaded_inputs.total_df,
            phospho_df=loaded_inputs.phospho_df,
            comparisons=comparisons,
            context="PhosphoDataset",
        )
        return cls._from_owned_dataset_inputs(dataset_inputs)

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
