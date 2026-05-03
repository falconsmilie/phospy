"""Dataset domain models."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import pandas as pd

from phospy._frame_ownership import (
    export_dataframe,
    export_optional_dataframe,
    own_dataframe,
    own_optional_dataframe,
)
from phospy.datasets.preprocessing.report_schema import (
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
from phospy.datasets.processing_state import DatasetProcessingState, RuvReadinessState
from phospy.errors.validation import DatasetValidationError
from phospy.provenance.models import RunProvenance
from phospy.references.models import Organism
from phospy.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.transformations.models import IntensityScaleState
from phospy.validation.common.dataframes import (
    require_columns,
    require_dataframe,
    require_exact_index_match,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_unique_columns,
)
from phospy.validation.transformations.state import IntensityScaleStateValidator

_INTENSITY_SCALE_STATE_VALIDATOR = IntensityScaleStateValidator()


@dataclass(frozen=True, slots=True)
class PreprocessingSiteAttritionSummary:
    """Compact preprocessing-owned site attrition counters."""

    input_rows: int
    output_rows: int
    rows_removed_during_preprocessing: int
    rows_removed_invalid_or_missing_site_identifiers: int
    duplicate_sites_merged_or_resolved: int


@dataclass(frozen=True, slots=True, init=False)
class DatasetPreprocessingReport:
    """Public provenance report for dataset preprocessing."""

    _row_counts: pd.DataFrame = field(init=False, repr=False)
    _operations: pd.DataFrame = field(init=False, repr=False)
    _row_audit: pd.DataFrame = field(init=False, repr=False)
    _duplicate_site_resolution: pd.DataFrame | None = field(init=False, repr=False)
    _metadata_conflicts: pd.DataFrame | None = field(init=False, repr=False)
    _comparison_group_stats: pd.DataFrame | None = field(init=False, repr=False)
    _comparison_pair_stats: pd.DataFrame | None = field(init=False, repr=False)

    def __init__(
        self,
        row_counts: pd.DataFrame,
        operations: pd.DataFrame,
        row_audit: pd.DataFrame,
        duplicate_site_resolution: pd.DataFrame | None = None,
        metadata_conflicts: pd.DataFrame | None = None,
        comparison_group_stats: pd.DataFrame | None = None,
        comparison_pair_stats: pd.DataFrame | None = None,
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
            duplicate_sites_merged_or_resolved = int(
                duplicate_site_resolution.loc[:, "site_id"]
                .astype("string")
                .str.strip()
                .dropna()
                .nunique()
            )
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
    ) -> DatasetPreprocessingReport:
        return cls(
            row_counts=row_counts,
            operations=operations,
            row_audit=row_audit,
            duplicate_site_resolution=duplicate_site_resolution,
            metadata_conflicts=metadata_conflicts,
            comparison_group_stats=comparison_group_stats,
            comparison_pair_stats=comparison_pair_stats,
            _assume_owned=True,
        )

    def _resolve_row_count(self, *, stage: str) -> int:
        row_counts = self._row_counts
        if row_counts.empty:
            return 0
        stage_mask = row_counts.loc[:, "stage"].astype(str) == str(stage)
        if not bool(stage_mask.any()):
            return 0
        return int(row_counts.loc[stage_mask, "output_rows"].iloc[-1])

    def _count_invalid_or_missing_identifier_drops(self) -> int:
        row_audit = self._row_audit
        if row_audit.empty:
            return 0
        dropped = row_audit.loc[row_audit.loc[:, "action"].astype(str) == "dropped", :]
        if dropped.empty:
            return 0
        reasons = dropped.loc[:, "reason"].astype("string").str.lower()
        site_id_reason = reasons.str.contains(
            "site identifier|site_id|site id",
            regex=True,
            na=False,
        )
        invalid_reason = reasons.str.contains(
            "invalid|missing|blank|empty",
            regex=True,
            na=False,
        )
        return int((site_id_reason & invalid_reason).sum())


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
    Site identity is strict at this boundary: canonical phospho row IDs must be
    coherent with `site_metadata.gene_symbol` / `site_metadata.site`.

    Provenance in this object describes owned internal state at creation time.
    Public export helpers return defensive snapshots; mutating exports does not
    mutate this owning dataset.
    """

    intensity_scale_state: IntensityScaleState
    processing_state: DatasetProcessingState
    organism: Organism | None = None
    preprocessing_report: DatasetPreprocessingReport | None = None
    provenance: RunProvenance | None = None
    _phospho: pd.DataFrame = field(init=False, repr=False)
    _site_metadata: pd.DataFrame = field(init=False, repr=False)
    _sample_metadata: pd.DataFrame | None = field(init=False, repr=False)
    _total: pd.DataFrame | None = field(init=False, repr=False)
    _comparisons: pd.DataFrame | None = field(init=False, repr=False)
    _init_payload: (
        tuple[
            pd.DataFrame,
            pd.DataFrame,
            pd.DataFrame | None,
            pd.DataFrame | None,
            pd.DataFrame | None,
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
        _assume_owned: bool = False,
    ) -> None:
        object.__setattr__(self, "intensity_scale_state", intensity_scale_state)
        object.__setattr__(self, "processing_state", processing_state)
        object.__setattr__(self, "organism", organism)
        object.__setattr__(self, "preprocessing_report", preprocessing_report)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(
            self,
            "_init_payload",
            (
                phospho,
                site_metadata,
                sample_metadata,
                total,
                comparisons,
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
        phospho_table = PhosphoIntensityMatrix(
            frame=phospho,
            _assume_owned=True,
        )
        site_metadata_table = SiteMetadataTable(
            frame=site_metadata,
            expected_index=phospho_table.frame.index,
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
        if self.organism is not None and not isinstance(self.organism, Organism):
            raise DatasetValidationError(
                "dataset.organism must be an Organism enum value or None"
            )
        validated_intensity_scale_state = _INTENSITY_SCALE_STATE_VALIDATOR.run(
            intensity_scale_state=self.intensity_scale_state,
            has_total_matrix=total_table is not None,
            require_established=True,
        )
        if not isinstance(self.processing_state, DatasetProcessingState):
            raise DatasetValidationError(
                "dataset.processing_state must be a DatasetProcessingState instance"
            )
        if not isinstance(self.processing_state.ruv_readiness, RuvReadinessState):
            raise DatasetValidationError(
                "dataset.processing_state.ruv_readiness must be a "
                "RuvReadinessState instance"
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
        if self.preprocessing_report is not None and not isinstance(
            self.preprocessing_report,
            DatasetPreprocessingReport,
        ):
            raise DatasetValidationError(
                "dataset.preprocessing_report must be DatasetPreprocessingReport "
                "or None"
            )
        if self.provenance is not None and not isinstance(
            self.provenance, RunProvenance
        ):
            raise DatasetValidationError(
                "dataset.provenance must be RunProvenance or None"
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
