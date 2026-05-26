"""Dataset domain models."""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TypeVar

import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.ownership import (
    borrow_dataframe,
    borrow_optional_dataframe,
    export_dataframe,
    export_optional_dataframe,
    own_dataframe,
    own_optional_dataframe,
)
from phospy.provenance.models import RunProvenance
from phospy.science.datasets.preprocessing.report_schema import (
    COMPARISON_GROUP_STATS_COLUMNS,
    COMPARISON_PAIR_STATS_COLUMNS,
    DUPLICATE_SITE_RESOLUTION_COLUMNS,
    METADATA_CONFLICT_COLUMNS,
    OPERATIONS_COLUMNS,
    ROW_AUDIT_COLUMNS,
    ROW_COUNTS_COLUMNS,
    ComparisonGroupStatsRow,
    ComparisonPairStatsRow,
    DuplicateSiteResolutionRow,
    MetadataConflictRow,
    PreprocessingOperationRow,
    PreprocessingRowAuditRow,
    PreprocessingRowCountRow,
    dataframe_from_comparison_group_stats_rows,
    dataframe_from_comparison_pair_stats_rows,
    dataframe_from_duplicate_site_resolution_rows,
    dataframe_from_metadata_conflict_rows,
    dataframe_from_operation_rows,
    dataframe_from_row_audit_rows,
    dataframe_from_row_count_rows,
    reorder_columns,
)
from phospy.science.datasets.processing_state import (
    DatasetProcessingState,
    RuvReadinessState,
)
from phospy.science.references.models import Organism
from phospy.science.transformations.models import IntensityScaleState
from phospy.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.validation.common.dataframes import (
    require_columns,
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
)
from phospy.validation.datasets.display_site_identity import (
    enforce_unique_display_site_identity_rows,
)
from phospy.validation.transformations.state import IntensityScaleStateValidator

_INTENSITY_SCALE_STATE_VALIDATOR = IntensityScaleStateValidator()
_SITE_ID_REASON_PATTERN = re.compile(r"site identifier|site_id|site id")
_INVALID_REASON_PATTERN = re.compile(r"invalid|missing|blank|empty")
_ExpectedType = TypeVar("_ExpectedType")


@dataclass(frozen=True, slots=True)
class PreprocessingSiteAttritionSummary:
    """Compact preprocessing-owned site attrition counters."""

    input_rows: int
    output_rows: int
    rows_removed_during_preprocessing: int
    rows_removed_invalid_or_missing_site_identifiers: int
    duplicate_sites_merged_or_resolved: int


@dataclass(frozen=True, slots=True)
class SiteSequenceResolutionReport:
    """Structured provenance summary for site-sequence origin and loss."""

    total_sites: int
    provided_by_input: int
    resolved_from_fasta: int
    resolved_from_reference: int
    unresolved: int
    conflicts: int
    conflict_policy: str
    final_sequence_complete_sites: int


@dataclass(frozen=True, slots=True, init=False)
class DatasetPreprocessingReport:
    """Public provenance report for dataset preprocessing.

    Internal `_borrow_*` accessors expose mutation-isolated borrowed snapshots
    for trusted internal read paths only.
    """

    _row_counts: pd.DataFrame = field(init=False, repr=False)
    _operations: pd.DataFrame = field(init=False, repr=False)
    _row_audit: pd.DataFrame = field(init=False, repr=False)
    _duplicate_site_resolution: pd.DataFrame | None = field(init=False, repr=False)
    _metadata_conflicts: pd.DataFrame | None = field(init=False, repr=False)
    _comparison_group_stats: pd.DataFrame | None = field(init=False, repr=False)
    _comparison_pair_stats: pd.DataFrame | None = field(init=False, repr=False)
    _site_sequence_resolution: SiteSequenceResolutionReport | None = field(
        init=False,
        repr=False,
    )

    def __init__(
        self,
        row_counts: pd.DataFrame,
        operations: pd.DataFrame,
        row_audit: pd.DataFrame,
        duplicate_site_resolution: pd.DataFrame | None = None,
        metadata_conflicts: pd.DataFrame | None = None,
        comparison_group_stats: pd.DataFrame | None = None,
        comparison_pair_stats: pd.DataFrame | None = None,
        site_sequence_resolution: SiteSequenceResolutionReport | None = None,
        _assume_owned: bool = False,
    ) -> None:
        row_counts = own_dataframe(
            row_counts,
            field_name="dataset.preprocessing_report.row_counts",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        operations = own_dataframe(
            operations,
            field_name="dataset.preprocessing_report.operations",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        row_audit = own_dataframe(
            row_audit,
            field_name="dataset.preprocessing_report.row_audit",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        duplicate_site_resolution = own_optional_dataframe(
            duplicate_site_resolution,
            field_name="dataset.preprocessing_report.duplicate_site_resolution",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        metadata_conflicts = own_optional_dataframe(
            metadata_conflicts,
            field_name="dataset.preprocessing_report.metadata_conflicts",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        comparison_group_stats = own_optional_dataframe(
            comparison_group_stats,
            field_name="dataset.preprocessing_report.comparison_group_stats",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        comparison_pair_stats = own_optional_dataframe(
            comparison_pair_stats,
            field_name="dataset.preprocessing_report.comparison_pair_stats",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        require_columns(
            row_counts,
            field_name="dataset.preprocessing_report.row_counts",
            required_columns=ROW_COUNTS_COLUMNS,
            error_type=DatasetValidationError,
        )
        require_columns(
            operations,
            field_name="dataset.preprocessing_report.operations",
            required_columns=OPERATIONS_COLUMNS,
            error_type=DatasetValidationError,
        )
        require_columns(
            row_audit,
            field_name="dataset.preprocessing_report.row_audit",
            required_columns=ROW_AUDIT_COLUMNS,
            error_type=DatasetValidationError,
        )
        if duplicate_site_resolution is not None:
            require_columns(
                duplicate_site_resolution,
                field_name="dataset.preprocessing_report.duplicate_site_resolution",
                required_columns=DUPLICATE_SITE_RESOLUTION_COLUMNS,
                error_type=DatasetValidationError,
            )
        if metadata_conflicts is not None:
            require_columns(
                metadata_conflicts,
                field_name="dataset.preprocessing_report.metadata_conflicts",
                required_columns=METADATA_CONFLICT_COLUMNS,
                error_type=DatasetValidationError,
            )
        if comparison_group_stats is not None:
            require_columns(
                comparison_group_stats,
                field_name="dataset.preprocessing_report.comparison_group_stats",
                required_columns=COMPARISON_GROUP_STATS_COLUMNS,
                error_type=DatasetValidationError,
            )
        if comparison_pair_stats is not None:
            require_columns(
                comparison_pair_stats,
                field_name="dataset.preprocessing_report.comparison_pair_stats",
                required_columns=COMPARISON_PAIR_STATS_COLUMNS,
                error_type=DatasetValidationError,
            )
        row_counts = reorder_columns(row_counts, expected_columns=ROW_COUNTS_COLUMNS)
        operations = reorder_columns(operations, expected_columns=OPERATIONS_COLUMNS)
        row_audit = reorder_columns(row_audit, expected_columns=ROW_AUDIT_COLUMNS)
        if duplicate_site_resolution is not None:
            duplicate_site_resolution = reorder_columns(
                duplicate_site_resolution,
                expected_columns=DUPLICATE_SITE_RESOLUTION_COLUMNS,
            )
        if metadata_conflicts is not None:
            metadata_conflicts = reorder_columns(
                metadata_conflicts,
                expected_columns=METADATA_CONFLICT_COLUMNS,
            )
        if comparison_group_stats is not None:
            comparison_group_stats = reorder_columns(
                comparison_group_stats,
                expected_columns=COMPARISON_GROUP_STATS_COLUMNS,
            )
        if comparison_pair_stats is not None:
            comparison_pair_stats = reorder_columns(
                comparison_pair_stats,
                expected_columns=COMPARISON_PAIR_STATS_COLUMNS,
            )
        object.__setattr__(self, "_row_counts", row_counts)
        object.__setattr__(self, "_operations", operations)
        object.__setattr__(self, "_row_audit", row_audit)
        object.__setattr__(
            self, "_duplicate_site_resolution", duplicate_site_resolution
        )
        object.__setattr__(self, "_metadata_conflicts", metadata_conflicts)
        object.__setattr__(self, "_comparison_group_stats", comparison_group_stats)
        object.__setattr__(self, "_comparison_pair_stats", comparison_pair_stats)
        object.__setattr__(self, "_site_sequence_resolution", site_sequence_resolution)

    @property
    def row_counts(self) -> pd.DataFrame:
        return export_dataframe(self._row_counts)

    @property
    def operations(self) -> pd.DataFrame:
        return export_dataframe(self._operations)

    @property
    def row_audit(self) -> pd.DataFrame:
        return export_dataframe(self._row_audit)

    @property
    def duplicate_site_resolution(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._duplicate_site_resolution)

    @property
    def metadata_conflicts(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._metadata_conflicts)

    @property
    def comparison_group_stats(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._comparison_group_stats)

    @property
    def comparison_pair_stats(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._comparison_pair_stats)

    @property
    def site_sequence_resolution(self) -> SiteSequenceResolutionReport | None:
        return self._site_sequence_resolution

    def _borrow_row_counts_frame(self) -> pd.DataFrame:
        """Package-private borrowed row-count table for internal workflows."""

        return borrow_dataframe(self._row_counts)

    def _borrow_operations_frame(self) -> pd.DataFrame:
        """Package-private borrowed operations table for internal workflows."""

        return borrow_dataframe(self._operations)

    def _borrow_row_audit_frame(self) -> pd.DataFrame:
        """Package-private borrowed row-audit table for internal workflows."""

        return borrow_dataframe(self._row_audit)

    def _borrow_duplicate_site_resolution_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed duplicate-resolution table for internals."""

        return borrow_optional_dataframe(self._duplicate_site_resolution)

    def _borrow_metadata_conflicts_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed metadata-conflict table for internals."""

        return borrow_optional_dataframe(self._metadata_conflicts)

    def _borrow_comparison_group_stats_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed comparison-group stats for internals."""

        return borrow_optional_dataframe(self._comparison_group_stats)

    def _borrow_comparison_pair_stats_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed comparison-pair stats for internals."""

        return borrow_optional_dataframe(self._comparison_pair_stats)

    def row_counts_dataframe(self) -> pd.DataFrame:
        """Return a row-count snapshot; mutating it does not mutate this report."""

        return export_dataframe(self._row_counts)

    def operations_dataframe(self) -> pd.DataFrame:
        """Return an operations snapshot; mutating it does not mutate this report."""

        return export_dataframe(self._operations)

    def row_audit_dataframe(self) -> pd.DataFrame:
        """Return a row-audit snapshot; mutating it does not mutate this report."""

        return export_dataframe(self._row_audit)

    def duplicate_site_resolution_dataframe(self) -> pd.DataFrame | None:
        """Return an optional snapshot isolated from this report."""

        return export_optional_dataframe(self._duplicate_site_resolution)

    def metadata_conflicts_dataframe(self) -> pd.DataFrame | None:
        """Return an optional snapshot isolated from this report."""

        return export_optional_dataframe(self._metadata_conflicts)

    def comparison_group_stats_dataframe(self) -> pd.DataFrame | None:
        """Return an optional snapshot isolated from this report."""

        return export_optional_dataframe(self._comparison_group_stats)

    def comparison_pair_stats_dataframe(self) -> pd.DataFrame | None:
        """Return an optional snapshot isolated from this report."""

        return export_optional_dataframe(self._comparison_pair_stats)

    def site_sequence_resolution_summary(self) -> SiteSequenceResolutionReport | None:
        """Return structured site-sequence provenance summary when available."""

        return self._site_sequence_resolution

    def site_attrition_summary(self) -> PreprocessingSiteAttritionSummary:
        """Return compact preprocessing-owned site attrition counters."""

        input_rows = int(self._resolve_row_count(stage="preprocessing_input"))
        output_rows = int(self._resolve_row_count(stage="preprocessing_complete"))
        duplicate_sites_merged_or_resolved = 0
        duplicate_site_resolution = self._duplicate_site_resolution
        if (
            duplicate_site_resolution is not None
            and not duplicate_site_resolution.empty
        ):
            site_id_values = duplicate_site_resolution["site_id"]
            site_id_value_list: list[object] = list(site_id_values.tolist())
            normalized_site_ids: set[str] = set()
            for raw_site_id in site_id_value_list:
                if _is_missing_value(raw_site_id):
                    continue
                normalized_site_ids.add(str(raw_site_id).strip())
            duplicate_sites_merged_or_resolved = int(len(normalized_site_ids))
        rows_removed_invalid_or_missing_site_identifiers = (
            self._count_invalid_or_missing_identifier_drops()
        )
        return PreprocessingSiteAttritionSummary(
            input_rows=input_rows,
            output_rows=output_rows,
            rows_removed_during_preprocessing=max(input_rows - output_rows, 0),
            rows_removed_invalid_or_missing_site_identifiers=(
                rows_removed_invalid_or_missing_site_identifiers
            ),
            duplicate_sites_merged_or_resolved=duplicate_sites_merged_or_resolved,
        )

    @classmethod
    def from_rows(
        cls,
        *,
        row_count_rows: Sequence[PreprocessingRowCountRow] = (),
        operation_rows: Sequence[PreprocessingOperationRow] = (),
        row_audit_rows: Sequence[PreprocessingRowAuditRow] = (),
        duplicate_site_resolution_rows: Sequence[DuplicateSiteResolutionRow] = (),
        metadata_conflict_rows: Sequence[MetadataConflictRow] = (),
        comparison_group_stats_rows: Sequence[ComparisonGroupStatsRow] = (),
        comparison_pair_stats_rows: Sequence[ComparisonPairStatsRow] = (),
        site_sequence_resolution: SiteSequenceResolutionReport | None = None,
    ) -> DatasetPreprocessingReport:
        return cls._from_owned(
            row_counts=dataframe_from_row_count_rows(row_count_rows),
            operations=dataframe_from_operation_rows(operation_rows),
            row_audit=dataframe_from_row_audit_rows(row_audit_rows),
            duplicate_site_resolution=dataframe_from_duplicate_site_resolution_rows(
                duplicate_site_resolution_rows
            ),
            metadata_conflicts=dataframe_from_metadata_conflict_rows(
                metadata_conflict_rows
            ),
            comparison_group_stats=dataframe_from_comparison_group_stats_rows(
                comparison_group_stats_rows
            ),
            comparison_pair_stats=dataframe_from_comparison_pair_stats_rows(
                comparison_pair_stats_rows
            ),
            site_sequence_resolution=site_sequence_resolution,
        )

    @classmethod
    def _from_owned(
        cls,
        *,
        row_counts: pd.DataFrame,
        operations: pd.DataFrame,
        row_audit: pd.DataFrame,
        duplicate_site_resolution: pd.DataFrame | None = None,
        metadata_conflicts: pd.DataFrame | None = None,
        comparison_group_stats: pd.DataFrame | None = None,
        comparison_pair_stats: pd.DataFrame | None = None,
        site_sequence_resolution: SiteSequenceResolutionReport | None = None,
    ) -> DatasetPreprocessingReport:
        return cls(
            row_counts=row_counts,
            operations=operations,
            row_audit=row_audit,
            duplicate_site_resolution=duplicate_site_resolution,
            metadata_conflicts=metadata_conflicts,
            comparison_group_stats=comparison_group_stats,
            comparison_pair_stats=comparison_pair_stats,
            site_sequence_resolution=site_sequence_resolution,
            _assume_owned=True,
        )

    def _resolve_row_count(self, *, stage: str) -> int:
        row_counts = self._row_counts
        if row_counts.empty:
            return 0
        stage_values: list[object] = list(row_counts["stage"].tolist())
        output_values: list[object] = list(row_counts["output_rows"].tolist())
        target_stage = str(stage)
        matched_values: list[object] = []
        for stage_value, output_value in zip(stage_values, output_values, strict=True):
            if str(stage_value) == target_stage:
                matched_values.append(output_value)
        if not matched_values:
            return 0
        last_value = matched_values[-1]
        if isinstance(last_value, bool):
            return int(last_value)
        if isinstance(last_value, (int, float, str)):
            return int(last_value)
        raise DatasetValidationError(
            "dataset.preprocessing_report.row_counts.output_rows must be "
            "int-like values"
        )

    def _count_invalid_or_missing_identifier_drops(self) -> int:
        row_audit = self._row_audit
        if row_audit.empty:
            return 0
        action_values: list[object] = list(row_audit["action"].tolist())
        reason_values: list[object] = list(row_audit["reason"].tolist())
        drop_count = 0
        for action_value, reason_value in zip(
            action_values, reason_values, strict=True
        ):
            if str(action_value) != "dropped":
                continue
            if _is_missing_value(reason_value):
                continue
            normalized_reason = str(reason_value).lower()
            if _SITE_ID_REASON_PATTERN.search(
                normalized_reason
            ) and _INVALID_REASON_PATTERN.search(normalized_reason):
                drop_count += 1
        return drop_count


@dataclass(frozen=True, slots=True, init=False)
class AnalysisReadyPhosphoDataset:
    """Public analysis-ready dataset contract.

    `phospho` stores the quantitative matrix after builder preprocessing policy
    has been applied. When total/protein correction is enabled in the builder
    lane, corrected values are represented directly in this matrix. When
    site-matrix construction is enabled in the builder lane, this matrix already
    reflects the constructed site-matrix-ready rows. Intermediate site-matrix
    artefacts remain private to preprocessing internals.
    Optional `comparisons` can carry builder-constructed dataset-level pairwise
    columns aligned to `phospho.index`.
    Site identity is strict at this boundary: standard phospho row IDs must be
    coherent with `site_metadata.gene_symbol` / `site_metadata.site`.

    Provenance in this object describes owned internal state at creation time.
    Public export helpers return defensive snapshots; mutating exports does not
    mutate this owning dataset.
    Internal `_borrow_*` accessors are reserved for trusted internal paths and
    return mutation-isolated borrowed snapshots.
    """

    intensity_scale_state: IntensityScaleState
    processing_state: DatasetProcessingState
    organism: Organism | None = None
    preprocessing_report: DatasetPreprocessingReport | None = None
    provenance: RunProvenance | None = None
    allow_opaque_site_values: bool = False
    _phospho: pd.DataFrame = field(init=False, repr=False)
    _site_metadata: pd.DataFrame = field(init=False, repr=False)
    _sample_metadata: pd.DataFrame | None = field(init=False, repr=False)
    _total: pd.DataFrame | None = field(init=False, repr=False)
    _comparisons: pd.DataFrame | None = field(init=False, repr=False)
    _allow_opaque_site_values: bool = field(init=False, repr=False, default=False)
    _init_payload: (
        tuple[
            pd.DataFrame,
            pd.DataFrame,
            pd.DataFrame | None,
            pd.DataFrame | None,
            pd.DataFrame | None,
            bool,
            bool,
        ]
        | None
    ) = field(init=False, repr=False, default=None)

    def __init__(
        self,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        provenance: RunProvenance | None = None,
        allow_opaque_site_values: bool = False,
        _assume_owned: bool = False,
    ) -> None:
        _require_instance(
            allow_opaque_site_values,
            expected_type=bool,
            error_message="dataset.allow_opaque_site_values must be a bool",
        )
        object.__setattr__(self, "intensity_scale_state", intensity_scale_state)
        object.__setattr__(self, "processing_state", processing_state)
        object.__setattr__(self, "organism", organism)
        object.__setattr__(self, "preprocessing_report", preprocessing_report)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "allow_opaque_site_values", allow_opaque_site_values)
        object.__setattr__(
            self,
            "_init_payload",
            (
                phospho,
                site_metadata,
                sample_metadata,
                total,
                comparisons,
                allow_opaque_site_values,
                _assume_owned,
            ),
        )
        self.__post_init__()

    def __post_init__(self) -> None:
        payload = self._init_payload
        if payload is None:
            raise DatasetValidationError(
                "dataset internal initialization payload missing"
            )
        (
            phospho,
            site_metadata,
            sample_metadata,
            total,
            comparisons,
            allow_opaque_site_values,
            _assume_owned,
        ) = payload
        phospho = own_dataframe(
            phospho,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        site_metadata = own_dataframe(
            site_metadata,
            field_name="dataset.site_metadata",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        sample_metadata = own_optional_dataframe(
            sample_metadata,
            field_name="dataset.sample_metadata",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        total = own_optional_dataframe(
            total,
            field_name="dataset.total",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        comparisons = own_optional_dataframe(
            comparisons,
            field_name="dataset.comparisons",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        enforce_unique_display_site_identity_rows(
            site_metadata=site_metadata,
            display_site_ids=pd.Series(
                phospho.index.tolist(),
                index=pd.Index(site_metadata.index),
                name="display_site_id",
                dtype="object",
            ),
            field_name="dataset.display_site_identity",
            error_type=DatasetValidationError,
        )
        phospho_table = PhosphoIntensityMatrix(
            frame=phospho,
            _assume_owned=True,
        )
        site_metadata_table = SiteMetadataTable(
            frame=site_metadata,
            expected_index=phospho_table.frame.index,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=True,
        )
        sample_metadata_table = (
            None
            if sample_metadata is None
            else SampleMetadataTable(
                frame=sample_metadata,
                expected_index=phospho_table.frame.columns,
                _assume_owned=True,
            )
        )
        total_table = (
            None
            if total is None
            else TotalProteinMatrix(
                frame=total,
                expected_sample_index=phospho_table.frame.columns,
                _assume_owned=True,
            )
        )
        if comparisons is not None:
            comparisons_frame = require_dataframe(
                comparisons,
                field_name="dataset.comparisons",
                allow_empty=True,
                error_type=DatasetValidationError,
            )
            require_non_empty_dataframe(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
            )
            require_numeric_dataframe(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
            )
            require_finite_numeric_dataframe(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
                allow_missing=False,
            )
            require_unique_columns(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
            )
            require_exact_index_match(
                left=comparisons_frame.index,
                right=phospho_table.frame.index,
                left_name="dataset.comparisons.index",
                right_name="dataset.phospho.index",
                error_type=DatasetValidationError,
            )
            comparisons = comparisons_frame
        validated_intensity_scale_state = _INTENSITY_SCALE_STATE_VALIDATOR.run(
            intensity_scale_state=self.intensity_scale_state,
            has_total_matrix=total_table is not None,
            require_established=True,
        )
        _require_optional_instance(
            self.organism,
            expected_type=Organism,
            error_message="dataset.organism must be an Organism enum value or None",
        )
        _require_instance(
            self.processing_state,
            expected_type=DatasetProcessingState,
            error_message=(
                "dataset.processing_state must be a DatasetProcessingState instance"
            ),
        )
        _require_instance(
            self.processing_state.ruv_readiness,
            expected_type=RuvReadinessState,
            error_message=(
                "dataset.processing_state.ruv_readiness must be a "
                "RuvReadinessState instance"
            ),
        )
        if (
            not self.processing_state.ruv_readiness.enabled
            and self.processing_state.ruv_readiness.ready
        ):
            raise DatasetValidationError(
                "dataset.processing_state.ruv_readiness.ready must be False when "
                "ruv_readiness.enabled is False"
            )
        if self.processing_state.intensity_scale != validated_intensity_scale_state:
            raise DatasetValidationError(
                "dataset.processing_state.intensity_scale must match "
                "dataset.intensity_scale_state"
            )
        if not bool(self.processing_state.missing_data.complete_matrix):
            raise DatasetValidationError(
                "dataset.processing_state.missing_data.complete_matrix must be True "
                "at AnalysisReadyPhosphoDataset boundary"
            )
        _require_optional_instance(
            self.preprocessing_report,
            expected_type=DatasetPreprocessingReport,
            error_message=(
                "dataset.preprocessing_report must be DatasetPreprocessingReport "
                "or None"
            ),
        )
        _require_optional_instance(
            self.provenance,
            expected_type=RunProvenance,
            error_message="dataset.provenance must be RunProvenance or None",
        )
        object.__setattr__(self, "_phospho", phospho_table.frame)
        object.__setattr__(self, "_site_metadata", site_metadata_table.frame)
        object.__setattr__(
            self,
            "_sample_metadata",
            None if sample_metadata_table is None else sample_metadata_table.frame,
        )
        object.__setattr__(
            self, "_total", None if total_table is None else total_table.frame
        )
        object.__setattr__(self, "_comparisons", comparisons)
        object.__setattr__(self, "_allow_opaque_site_values", allow_opaque_site_values)
        object.__setattr__(self, "allow_opaque_site_values", allow_opaque_site_values)
        object.__setattr__(
            self, "intensity_scale_state", validated_intensity_scale_state
        )
        object.__setattr__(self, "_init_payload", None)

    @property
    def phospho(self) -> pd.DataFrame:
        return export_dataframe(self._phospho)

    @property
    def site_metadata(self) -> pd.DataFrame:
        return export_dataframe(self._site_metadata)

    @property
    def sample_metadata(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._sample_metadata)

    @property
    def total(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._total)

    @property
    def comparisons(self) -> pd.DataFrame | None:
        return export_optional_dataframe(self._comparisons)

    @property
    def opaque_site_values_allowed(self) -> bool:
        return bool(self._allow_opaque_site_values)

    def _borrow_phospho_frame(self) -> pd.DataFrame:
        """Package-private borrowed phospho matrix for internal workflows."""

        return borrow_dataframe(self._phospho)

    def _borrow_site_metadata_frame(self) -> pd.DataFrame:
        """Package-private borrowed site-metadata table for internals."""

        return borrow_dataframe(self._site_metadata)

    def _borrow_sample_metadata_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed sample-metadata table for internals."""

        return borrow_optional_dataframe(self._sample_metadata)

    def _borrow_total_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed total-protein table for internals."""

        return borrow_optional_dataframe(self._total)

    def _borrow_comparisons_frame(self) -> pd.DataFrame | None:
        """Package-private borrowed comparisons table for internals."""

        return borrow_optional_dataframe(self._comparisons)

    @classmethod
    def _from_owned(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
        provenance: RunProvenance | None = None,
        allow_opaque_site_values: bool = False,
    ) -> AnalysisReadyPhosphoDataset:
        return cls(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            organism=organism,
            preprocessing_report=preprocessing_report,
            provenance=provenance,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=True,
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Return a phospho snapshot; mutating it does not mutate this dataset."""

        return export_dataframe(self._phospho)

    def site_metadata_dataframe(self) -> pd.DataFrame:
        """Return a site-metadata snapshot isolated from this dataset."""

        return export_dataframe(self._site_metadata)

    def sample_metadata_dataframe(self) -> pd.DataFrame | None:
        """Return an optional sample-metadata snapshot isolated from this dataset."""

        return export_optional_dataframe(self._sample_metadata)

    def total_dataframe(self) -> pd.DataFrame | None:
        """Return an optional total-protein snapshot isolated from this dataset."""

        return export_optional_dataframe(self._total)

    def comparisons_dataframe(self) -> pd.DataFrame | None:
        """Return an optional comparisons snapshot isolated from this dataset."""

        return export_optional_dataframe(self._comparisons)


def _is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _require_instance(
    value: object,
    *,
    expected_type: type[_ExpectedType],
    error_message: str,
) -> None:
    if not isinstance(value, expected_type):
        raise DatasetValidationError(error_message)


def _require_optional_instance(
    value: object | None,
    *,
    expected_type: type[_ExpectedType],
    error_message: str,
) -> None:
    if value is None:
        return
    _require_instance(
        value,
        expected_type=expected_type,
        error_message=error_message,
    )
