from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..errors import InputCompatibilityError
from ..internal.constants import (
    PHOSPHO_GENE_COLUMN,
    ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY,
    ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY,
    ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY,
    ROW_DROP_DUPLICATE_SITE_STRATEGY_KEY,
    ROW_DROP_INPUT_ROWS_KEY,
    ROW_DROP_MISSING_DATA_POLICY_KEY,
    ROW_DROP_REQUIRED_OBSERVED_COUNT_KEY,
    ROW_DROP_RETAINED_ROWS_KEY,
    SITE_MATRIX_GENE_COLUMN,
    SITE_MATRIX_ID_COLUMN,
    ComparisonSpec,
)
from ..internal.pandas_copy import detached_frame_copy, detached_series_copy
from ..internal.types import (
    DUPLICATE_SITE_STRATEGY_MAX_MEAN_SIGNAL,
    SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
)
from ..preprocessing.core import CorePreprocessingConfig, CoreProcessingResult
from ..preprocessing.dataset import DatasetPreprocessing
from ..validation.values.identifiers import parse_canonical_site_id
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

_SITE_TO_PROTEIN_FALLBACK_POLICY_STRICT: str = "strict"
_SITE_TO_PROTEIN_FALLBACK_POLICY_METADATA: str = "metadata"
_SITE_TO_PROTEIN_FALLBACK_POLICIES: tuple[str, ...] = (
    _SITE_TO_PROTEIN_FALLBACK_POLICY_STRICT,
    _SITE_TO_PROTEIN_FALLBACK_POLICY_METADATA,
)
_STRICT_SITE_TO_PROTEIN_METADATA_COLUMNS: tuple[str, ...] = ("protein_id",)
_DEFAULT_FALLBACK_SITE_TO_PROTEIN_METADATA_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "protein",
    SITE_MATRIX_GENE_COLUMN,
    PHOSPHO_GENE_COLUMN,
)
_GENE_SYMBOL_METADATA_COLUMNS: tuple[str, ...] = (
    SITE_MATRIX_GENE_COLUMN,
    PHOSPHO_GENE_COLUMN,
)


@dataclass(frozen=True, slots=True)
class _SiteToProteinResolutionPolicy:
    fallback_policy: str
    candidate_columns: tuple[str, ...]


def _resolve_site_to_protein_policy(
    *,
    metadata_columns: Sequence[str] | None,
    fallback_policy: str,
) -> _SiteToProteinResolutionPolicy:
    resolved_fallback_policy = str(fallback_policy).strip().lower()
    if resolved_fallback_policy not in _SITE_TO_PROTEIN_FALLBACK_POLICIES:
        allowed = ", ".join(_SITE_TO_PROTEIN_FALLBACK_POLICIES)
        msg = (
            f"fallback_policy must be one of: {allowed}. Received: {fallback_policy!r}"
        )
        raise ValueError(msg)

    if resolved_fallback_policy == _SITE_TO_PROTEIN_FALLBACK_POLICY_STRICT:
        if metadata_columns is not None:
            msg = "metadata_columns is only supported when fallback_policy='metadata'."
            raise ValueError(msg)
        resolved_candidates = _STRICT_SITE_TO_PROTEIN_METADATA_COLUMNS
    else:
        resolved_candidates = (
            _DEFAULT_FALLBACK_SITE_TO_PROTEIN_METADATA_COLUMNS
            if metadata_columns is None
            else tuple(metadata_columns)
        )

    candidate_columns: list[str] = []
    for column in resolved_candidates:
        if column is None:
            continue
        normalized_column = str(column).strip()
        if normalized_column:
            candidate_columns.append(normalized_column)

    if not candidate_columns:
        msg = (
            "metadata_columns must contain at least one non-empty column name "
            "when resolving site-to-protein mapping."
        )
        raise ValueError(msg)

    return _SiteToProteinResolutionPolicy(
        fallback_policy=resolved_fallback_policy,
        candidate_columns=tuple(candidate_columns),
    )


def _resolve_aligned_site_metadata_indices(
    *,
    phospho_index: pd.Index,
    metadata_index: pd.Index,
) -> tuple[pd.Index, pd.Index]:
    site_index = pd.Index(
        phospho_index.astype("string"),
        name=SITE_MATRIX_ID_COLUMN,
    )
    site_metadata_index = pd.Index(
        metadata_index.astype("string"),
        name=SITE_MATRIX_ID_COLUMN,
    )
    if not site_index.equals(site_metadata_index):
        msg = (
            "AnalysisReadyPhosphoDataset.site_metadata must remain aligned with "
            "phospho_matrix to resolve site-to-protein mapping."
        )
        raise InputCompatibilityError(msg)
    return site_index, site_metadata_index


def _format_incomplete_site_to_protein_column_diagnostic(
    *,
    candidate_column: str,
    metadata_index: pd.Index,
    invalid_mask: pd.Series,
) -> str:
    invalid_site_ids = metadata_index[invalid_mask].astype(str).tolist()
    preview = ", ".join(invalid_site_ids[:3])
    diagnostic = f"column '{candidate_column}' has missing/empty values for: {preview}"
    if len(invalid_site_ids) > 3:
        diagnostic += ", ..."
    return diagnostic


def _evaluate_candidate_metadata_column_completeness(
    *,
    values: pd.Series,
    metadata_index: pd.Index,
    candidate_column: str,
) -> tuple[pd.Series | None, str | None]:
    normalized_values = values.astype("string")
    stripped_values = normalized_values.str.strip()
    invalid_mask = normalized_values.isna() | stripped_values.eq("")
    if bool(invalid_mask.any()):
        return None, _format_incomplete_site_to_protein_column_diagnostic(
            candidate_column=candidate_column,
            metadata_index=metadata_index,
            invalid_mask=invalid_mask,
        )
    return stripped_values, None


def _parse_canonical_entities_for_site_index(site_index: pd.Index) -> pd.Series:
    canonical_entities: list[object] = []
    for site_id in site_index.astype(str).tolist():
        parsed_site_id = parse_canonical_site_id(site_id)
        if parsed_site_id is None:
            canonical_entities.append(pd.NA)
            continue
        canonical_entity, _ = parsed_site_id
        canonical_entities.append(canonical_entity)
    return pd.Series(canonical_entities, index=site_index, dtype="string")


def _find_ambiguous_fallback_identifiers(
    *,
    fallback_values: pd.Series,
    canonical_entities: pd.Series,
) -> list[str]:
    ambiguity_frame = pd.DataFrame(
        {
            "fallback_value": fallback_values.astype("string"),
            "canonical_entity": canonical_entities.astype("string"),
        },
        copy=False,
    )
    parseable_rows = ambiguity_frame.loc[ambiguity_frame["canonical_entity"].notna()]
    if parseable_rows.empty:
        return []

    canonical_entity_counts = parseable_rows.groupby("fallback_value", sort=True)[
        "canonical_entity"
    ].nunique(dropna=True)
    ambiguous_identifiers = canonical_entity_counts.loc[
        canonical_entity_counts.gt(1)
    ].index.astype(str)
    return sorted(ambiguous_identifiers.tolist())


def _enforce_gene_symbol_fallback_policy(
    *,
    candidate_column: str,
    allow_gene_symbol_fallback: bool,
) -> None:
    if allow_gene_symbol_fallback:
        warnings.warn(
            (
                "Gene-symbol site-to-protein fallback is enabled for "
                f"column '{candidate_column}'. This can collapse biologically "
                "distinct proteins."
            ),
            category=UserWarning,
            stacklevel=2,
        )
        return

    msg = (
        "Gene-symbol site-to-protein fallback is disabled by default "
        "because it can collapse distinct proteins. "
        f"Resolved fallback column: '{candidate_column}'. "
        "Provide explicit site_to_protein, include a 'protein_id' "
        "metadata column, or set allow_gene_symbol_fallback=True "
        "to opt in."
    )
    raise InputCompatibilityError(msg)


def _validate_ambiguous_fallback_identifiers(
    *,
    candidate_column: str,
    ambiguous_identifiers: Sequence[str],
    allow_ambiguous_fallback: bool,
) -> None:
    if not ambiguous_identifiers:
        return

    preview = ", ".join(ambiguous_identifiers[:3])
    suffix = ", ..." if len(ambiguous_identifiers) > 3 else ""
    message = (
        "Ambiguous site-to-protein metadata mapping detected: "
        f"column '{candidate_column}' maps one fallback identifier to "
        "multiple canonical protein IDs inferred from site IDs. "
        f"Ambiguous identifiers: {preview}{suffix}"
    )
    if not allow_ambiguous_fallback:
        raise InputCompatibilityError(message)
    warnings.warn(
        f"{message}. Proceeding because allow_ambiguous_fallback=True.",
        category=UserWarning,
        stacklevel=2,
    )


def _build_site_to_protein_mapping_series(
    *,
    site_index: pd.Index,
    resolved_values: pd.Series,
) -> pd.Series:
    resolved_index = pd.Index(
        site_index.astype(str),
        dtype=object,
        name=SITE_MATRIX_ID_COLUMN,
    )
    return pd.Series(
        resolved_values.astype(object).tolist(),
        index=resolved_index,
        dtype=object,
        name="protein_id",
    )


def _raise_missing_site_to_protein_columns_error(
    *,
    fallback_policy: str,
    candidate_columns: Sequence[str],
    available_columns: pd.Index,
) -> None:
    checked_preview = ", ".join(candidate_columns)
    available_preview = ", ".join(str(column) for column in available_columns[:5])
    if len(available_columns) > 5:
        available_preview += ", ..."

    if fallback_policy == _SITE_TO_PROTEIN_FALLBACK_POLICY_STRICT:
        msg = (
            "AnalysisReadyPhosphoDataset.site_metadata does not include the "
            "required strict site-to-protein column 'protein_id'. "
            "Strict mode does not allow metadata fallback columns. "
        )
    else:
        msg = (
            "AnalysisReadyPhosphoDataset.site_metadata does not include a "
            "usable site-to-protein column. "
            f"Checked columns: {checked_preview}. "
        )

    if available_preview:
        msg += f"Available columns: {available_preview}. "
    else:
        msg += "site_metadata has no columns. "
    msg += "Provide site_to_protein explicitly when running signalome analysis."
    raise InputCompatibilityError(msg)


def _raise_incomplete_site_to_protein_mapping_error(
    *,
    candidate_columns: Sequence[str],
    diagnostics: Sequence[str],
) -> None:
    checked_preview = ", ".join(candidate_columns)
    diagnostic_message = "; ".join(diagnostics)
    msg = (
        "AnalysisReadyPhosphoDataset.site_metadata does not contain a complete "
        "site-to-protein mapping. "
        f"Checked columns: {checked_preview}. {diagnostic_message}. "
        "Provide site_to_protein explicitly when running signalome analysis."
    )
    raise InputCompatibilityError(msg)


@dataclass(slots=True)
class CoreInputs:
    """Owned mutable dataset tables used by the core preprocessing pipeline."""

    total_df: pd.DataFrame
    phospho_df: pd.DataFrame

    def copy_frames(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return detached copies suitable for caller-owned mutation."""
        return detached_frame_copy(self.total_df), detached_frame_copy(self.phospho_df)


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
    missing_data_policy: str = SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
    required_observed_count: int = 0
    duplicate_site_strategy: str = DUPLICATE_SITE_STRATEGY_MAX_MEAN_SIGNAL

    @classmethod
    def from_mapping(
        cls,
        row_drop_stats: Mapping[str, int | str],
    ) -> AnalysisReadySiteMatrixStats:
        return cls(
            input_rows=int(row_drop_stats.get(ROW_DROP_INPUT_ROWS_KEY, 0)),
            dropped_missing_sequence=int(
                row_drop_stats.get(ROW_DROP_DROPPED_MISSING_SEQUENCE_KEY, 0)
            ),
            dropped_incomplete_values=int(
                row_drop_stats.get(ROW_DROP_DROPPED_INCOMPLETE_VALUES_KEY, 0)
            ),
            missing_data_policy=str(
                row_drop_stats.get(
                    ROW_DROP_MISSING_DATA_POLICY_KEY,
                    SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
                )
            ),
            required_observed_count=int(
                row_drop_stats.get(ROW_DROP_REQUIRED_OBSERVED_COUNT_KEY, 0)
            ),
            deduplicated_site_rows=int(
                row_drop_stats.get(ROW_DROP_DEDUPLICATED_SITE_ROWS_KEY, 0)
            ),
            retained_rows=int(row_drop_stats.get(ROW_DROP_RETAINED_ROWS_KEY, 0)),
            duplicate_site_strategy=str(
                row_drop_stats.get(
                    ROW_DROP_DUPLICATE_SITE_STRATEGY_KEY,
                    DUPLICATE_SITE_STRATEGY_MAX_MEAN_SIGNAL,
                )
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
        detaches caller-managed pandas inputs once, then delegates to the owned
        fast path used by internal preprocessing results.
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
            phospho_matrix=detached_frame_copy(phospho_matrix),
            site_metadata=detached_frame_copy(site_metadata),
            site_sequences=detached_series_copy(site_sequences),
            phospho_corrected=detached_frame_copy(phospho_corrected),
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
        fallback_policy: str = _SITE_TO_PROTEIN_FALLBACK_POLICY_STRICT,
        allow_gene_symbol_fallback: bool = False,
        allow_ambiguous_fallback: bool = False,
    ) -> pd.Series:
        """Resolve an aligned site-to-protein mapping from site metadata.

        Parameters
        ----------
        metadata_columns:
            Candidate metadata columns to evaluate when ``fallback_policy`` is
            ``"metadata"``.
        fallback_policy:
            Mapping policy:
            - ``"strict"`` (default): require ``protein_id`` only.
            - ``"metadata"``: opt in to metadata fallback columns.
        allow_gene_symbol_fallback:
            Opt in to gene-symbol fallback when using ``fallback_policy="metadata"``.
            Disabled by default because gene symbols can collapse distinct proteins.
        allow_ambiguous_fallback:
            Opt in to ambiguous fallback values that map to multiple canonical
            protein IDs when those IDs are parseable from site IDs.
        """
        policy = _resolve_site_to_protein_policy(
            metadata_columns=metadata_columns,
            fallback_policy=fallback_policy,
        )
        site_index, metadata_index = _resolve_aligned_site_metadata_indices(
            phospho_index=self.phospho_matrix.index,
            metadata_index=self.site_metadata.index,
        )

        available_columns = {
            str(column): column for column in self.site_metadata.columns
        }
        checked_columns: list[str] = []
        incomplete_column_diagnostics: list[str] = []
        canonical_entities = (
            _parse_canonical_entities_for_site_index(site_index)
            if policy.fallback_policy == _SITE_TO_PROTEIN_FALLBACK_POLICY_METADATA
            else None
        )

        for candidate in policy.candidate_columns:
            original_column = available_columns.get(candidate)
            if original_column is None:
                continue
            checked_columns.append(candidate)

            raw_values = self.site_metadata.loc[:, original_column]
            stripped_values, incomplete_diagnostic = (
                _evaluate_candidate_metadata_column_completeness(
                    values=raw_values,
                    metadata_index=metadata_index,
                    candidate_column=candidate,
                )
            )
            if stripped_values is None:
                if incomplete_diagnostic is not None:
                    incomplete_column_diagnostics.append(incomplete_diagnostic)
                continue

            if (
                policy.fallback_policy == _SITE_TO_PROTEIN_FALLBACK_POLICY_METADATA
                and candidate in _GENE_SYMBOL_METADATA_COLUMNS
            ):
                _enforce_gene_symbol_fallback_policy(
                    candidate_column=candidate,
                    allow_gene_symbol_fallback=allow_gene_symbol_fallback,
                )

            if (
                policy.fallback_policy == _SITE_TO_PROTEIN_FALLBACK_POLICY_METADATA
                and candidate != "protein_id"
                and canonical_entities is not None
            ):
                ambiguous_identifiers = _find_ambiguous_fallback_identifiers(
                    fallback_values=stripped_values,
                    canonical_entities=canonical_entities,
                )
                _validate_ambiguous_fallback_identifiers(
                    candidate_column=candidate,
                    ambiguous_identifiers=ambiguous_identifiers,
                    allow_ambiguous_fallback=allow_ambiguous_fallback,
                )

            return _build_site_to_protein_mapping_series(
                site_index=site_index,
                resolved_values=stripped_values,
            )

        if not checked_columns:
            _raise_missing_site_to_protein_columns_error(
                fallback_policy=policy.fallback_policy,
                candidate_columns=policy.candidate_columns,
                available_columns=self.site_metadata.columns,
            )
        _raise_incomplete_site_to_protein_mapping_error(
            candidate_columns=policy.candidate_columns,
            diagnostics=incomplete_column_diagnostics,
        )


class PhosphoDataset:
    """Mutable workspace that owns owned phosphoproteomics input tables.

    `PhosphoDataset` is an owned in-memory processing workspace, not an immutable
    snapshot. It validates raw constructor inputs at the boundary, stores owned
    mutable pandas tables, and exposes detached copies by default.

    Ownership policy: copy once at external construction boundaries, then pass
    owned mutable frames through internal preprocessing paths without copying
    again unless a caller explicitly asks for detached copies.

    Prefer `total_df_copy` / `phospho_df_copy` or `copy_inputs()` for caller-owned
    inspection, export, and other read-oriented work. For expert workflows that
    intentionally mutate this dataset's owned state, use
    `to_mutable_frames_unsafe()`.
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
        """Return detached dataset inputs for safe caller-owned inspection.

        The returned ``CoreInputs`` bundle contains detached copies of the owned
        workspace frames. Mutating it does not change this dataset.
        """
        total_df, phospho_df = self._inputs.copy_frames()
        return CoreInputs(total_df=total_df, phospho_df=phospho_df)

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
    def total_df_copy(self) -> pd.DataFrame:
        """Return a detached copy of the owned total-protein table."""
        return detached_frame_copy(self._inputs.total_df)

    @property
    def phospho_df_copy(self) -> pd.DataFrame:
        """Return a detached copy of the owned phosphoproteomics table."""
        return detached_frame_copy(self._inputs.phospho_df)

    def copy_inputs(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return detached copies suitable for caller-owned mutation and read use."""
        return self.total_df_copy, self.phospho_df_copy

    def to_mutable_frames_unsafe(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return owned mutable workspace frames for expert in-place mutation.

        Warning: mutating these frames mutates this dataset's internal state and
        can invalidate assumptions in downstream code.
        """
        return self._inputs.total_df, self._inputs.phospho_df

    @property
    def preprocessing(self) -> DatasetPreprocessing:
        """Return the bound preprocessing facade for this dataset workspace."""
        total_df, phospho_df = self.to_mutable_frames_unsafe()
        return DatasetPreprocessing.from_owned(
            total_df=total_df,
            phospho_df=phospho_df,
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
