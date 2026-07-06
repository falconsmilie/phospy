"""Dataset preprocessing-state summary models."""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TypeVar, cast

import numpy as np
import numpy.typing as npt
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
from phospy.science.datasets._processing_state.json_contracts import (
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
    TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1,
    JsonPrimitive,
    JsonValue,
)
from phospy.science.datasets._processing_state.missing_data import (
    MissingDataDiagnostics,
    MissingDataDiagnosticsV1,
)
from phospy.science.datasets._processing_state.models import (
    ComparisonState,
    DatasetProcessingState,
    MissingDataState,
    NormalisationState,
    RuvReadinessState,
    SiteMatrixState,
    SiteSequenceResolutionRowDiagnostic,
    SiteSequenceResolutionState,
    TotalProteinCorrectionState,
    default_ruv_readiness_state,
)
from phospy.science.datasets._processing_state.total_protein import (
    TotalProteinCorrectionDiagnostics,
    TotalProteinCorrectionDiagnosticsV1,
)
from phospy.science.datasets.preprocessing.batch_correction import (
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationReport,
)
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
from phospy.validation.common.dataframes import (
    require_columns,
    require_dataframe,
    require_exact_index_match,
)

_ExpectedType = TypeVar("_ExpectedType")
_SITE_ID_REASON_PATTERN = re.compile(r"site identifier|site_id|site id")
_INVALID_REASON_PATTERN = re.compile(r"invalid|missing|blank|empty")
IMPUTATION_FEATURE_METADATA_COLUMNS = (
    "imputed_cell_count",
    "observed_cell_count",
    "imputed_fraction",
)
IMPUTATION_OBSERVATION_SUMMARY_COLUMNS = (
    "feature_id",
    "observed_cell_count",
    "imputed_cell_count",
    "total_analysed_cell_count",
    "imputed_fraction",
)


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

    Internal `_borrow_*` accessors return mutation-isolated internal table
    snapshots for trusted internal read paths only.
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
    _batch_correction: BatchCorrectionReport | None = field(init=False, repr=False)
    _protein_aware_preparation: ProteinAwarePreparationReport | None = field(
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
        batch_correction: BatchCorrectionReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationReport | None = None,
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
        object.__setattr__(self, "_batch_correction", batch_correction)
        _require_optional_instance(
            protein_aware_preparation,
            expected_type=ProteinAwarePreparationReport,
            error_message=(
                "dataset.preprocessing_report.protein_aware_preparation must be "
                "ProteinAwarePreparationReport or None"
            ),
        )
        object.__setattr__(
            self,
            "_protein_aware_preparation",
            protein_aware_preparation,
        )

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

    @property
    def batch_correction(self) -> BatchCorrectionReport | None:
        return self._batch_correction

    @property
    def protein_aware_preparation(self) -> ProteinAwarePreparationReport | None:
        return self._protein_aware_preparation

    def _borrow_row_counts_frame(self) -> pd.DataFrame:
        """Package-private row-count snapshot for internal read paths."""

        return borrow_dataframe(self._row_counts)

    def _borrow_operations_frame(self) -> pd.DataFrame:
        """Package-private operations snapshot for internal read paths."""

        return borrow_dataframe(self._operations)

    def _borrow_row_audit_frame(self) -> pd.DataFrame:
        """Package-private row-audit snapshot for internal read paths."""

        return borrow_dataframe(self._row_audit)

    def _borrow_duplicate_site_resolution_frame(self) -> pd.DataFrame | None:
        """Package-private duplicate-resolution snapshot for internals."""

        return borrow_optional_dataframe(self._duplicate_site_resolution)

    def _borrow_metadata_conflicts_frame(self) -> pd.DataFrame | None:
        """Package-private metadata-conflict snapshot for internals."""

        return borrow_optional_dataframe(self._metadata_conflicts)

    def _borrow_comparison_group_stats_frame(self) -> pd.DataFrame | None:
        """Package-private comparison-group stats snapshot for internals."""

        return borrow_optional_dataframe(self._comparison_group_stats)

    def _borrow_comparison_pair_stats_frame(self) -> pd.DataFrame | None:
        """Package-private comparison-pair stats snapshot for internals."""

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

    def batch_correction_summary(self) -> BatchCorrectionReport | None:
        """Return structured batch-correction provenance when available."""

        return self._batch_correction

    def protein_aware_preparation_summary(
        self,
    ) -> ProteinAwarePreparationReport | None:
        """Return protein-aware preparation provenance when available."""

        return self._protein_aware_preparation

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
        batch_correction: BatchCorrectionReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationReport | None = None,
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
            batch_correction=batch_correction,
            protein_aware_preparation=protein_aware_preparation,
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
        batch_correction: BatchCorrectionReport | None = None,
        protein_aware_preparation: ProteinAwarePreparationReport | None = None,
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
            batch_correction=batch_correction,
            protein_aware_preparation=protein_aware_preparation,
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
class ImputationObservationMetadata:
    """Dataset-owned originally-observed vs imputed-cell metadata.

    The internal mask is aligned to `dataset.phospho`: rows are phosphosite
    features, columns are samples, and True means the value was originally
    observed rather than imputed. Public accessors return defensive snapshots.
    """

    _observed_mask: pd.DataFrame = field(init=False, repr=False)
    _feature_summary: pd.DataFrame = field(init=False, repr=False)

    def __init__(
        self,
        *,
        observed_mask: pd.DataFrame,
        phospho_index: pd.Index,
        sample_index: pd.Index,
        _assume_owned: bool = False,
    ) -> None:
        mask = own_dataframe(
            observed_mask,
            field_name="dataset.imputation_observation_mask",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        require_dataframe(
            mask,
            field_name="dataset.imputation_observation_mask",
            allow_empty=False,
            error_type=DatasetValidationError,
        )
        require_exact_index_match(
            left=mask.index,
            right=phospho_index,
            left_name="dataset.imputation_observation_mask.index",
            right_name="dataset.phospho.index",
            error_type=DatasetValidationError,
        )
        require_exact_index_match(
            left=mask.columns,
            right=sample_index,
            left_name="dataset.imputation_observation_mask.columns",
            right_name="dataset.phospho.columns",
            error_type=DatasetValidationError,
        )
        _require_boolean_observation_mask(mask)

        observed = pd.DataFrame(
            mask.to_numpy(dtype=bool),
            index=mask.index,
            columns=mask.columns,
        )
        observed.index.name = mask.index.name
        observed.columns.name = mask.columns.name
        feature_summary = _build_imputation_feature_summary(observed)
        object.__setattr__(self, "_observed_mask", observed)
        object.__setattr__(self, "_feature_summary", feature_summary)

    @property
    def feature_summary(self) -> pd.DataFrame:
        return export_dataframe(self._feature_summary)

    @property
    def observed_mask(self) -> pd.DataFrame:
        return export_dataframe(self._observed_mask)

    def feature_summary_dataframe(self) -> pd.DataFrame:
        """Return per-feature imputation counts isolated from this metadata."""

        return export_dataframe(self._feature_summary)

    def observed_mask_dataframe(self) -> pd.DataFrame:
        """Return a defensive observed-cell mask snapshot."""

        return export_dataframe(self._observed_mask)

    def feature_observation_summary_dataframe(
        self,
        *,
        feature_ids: Sequence[object],
        sample_ids: Sequence[object],
    ) -> pd.DataFrame:
        """Return feature-level observation counts for a requested subset."""

        requested_feature_ids = _requested_label_list(
            feature_ids,
            field_name="dataset.imputation_observation_summary.feature_ids",
        )
        requested_sample_ids = _requested_label_list(
            sample_ids,
            field_name="dataset.imputation_observation_summary.sample_ids",
        )
        _require_requested_labels_present(
            requested_feature_ids,
            available_labels=self._observed_mask.index,
            field_name="dataset.imputation_observation_summary.feature_ids",
            available_field_name="dataset.imputation_observation_mask.index",
        )
        _require_requested_labels_present(
            requested_sample_ids,
            available_labels=self._observed_mask.columns,
            field_name="dataset.imputation_observation_summary.sample_ids",
            available_field_name="dataset.imputation_observation_mask.columns",
        )
        observed_values = self._observed_mask.to_numpy(dtype=bool)
        feature_positions = _label_positions(
            requested_feature_ids,
            available_labels=self._observed_mask.index,
        )
        sample_positions = _label_positions(
            requested_sample_ids,
            available_labels=self._observed_mask.columns,
        )
        observed_subset = pd.DataFrame(
            observed_values[np.ix_(feature_positions, sample_positions)],
            index=pd.Index(
                requested_feature_ids,
                name=self._observed_mask.index.name,
            ),
            columns=pd.Index(
                requested_sample_ids,
                name=self._observed_mask.columns.name,
            ),
        )
        summary = _build_imputation_observation_summary(observed_subset)
        return export_dataframe(summary)

    def aggregated_observed_mask_dataframe(
        self,
        *,
        sample_groups: Sequence[tuple[object, Sequence[object]]],
    ) -> pd.DataFrame:
        """Return an observed-cell mask collapsed to requested sample groups."""

        if not sample_groups:
            raise DatasetValidationError(
                "dataset.imputation_observation_mask sample_groups must contain "
                "at least one sample group"
            )
        observed_values = self._observed_mask.to_numpy(dtype=bool)
        aggregated_column_values: list[np.ndarray] = []
        output_labels: list[object] = []
        for output_label, input_sample_ids in sample_groups:
            requested_sample_ids = _requested_label_list(
                input_sample_ids,
                field_name=(
                    "dataset.imputation_observation_mask sample_groups."
                    f"{output_label!r}.sample_ids"
                ),
            )
            _require_requested_labels_present(
                requested_sample_ids,
                available_labels=self._observed_mask.columns,
                field_name=(
                    "dataset.imputation_observation_mask sample_groups."
                    f"{output_label!r}.sample_ids"
                ),
                available_field_name="dataset.imputation_observation_mask.columns",
            )
            sample_positions = _label_positions(
                requested_sample_ids,
                available_labels=self._observed_mask.columns,
            )
            collapsed = np.all(observed_values[:, sample_positions], axis=1)
            aggregated_column_values.append(collapsed.astype(bool))
            output_labels.append(output_label)
        aggregated_values = np.column_stack(aggregated_column_values).astype(bool)
        aggregated = pd.DataFrame(
            aggregated_values,
            index=pd.Index(
                self._observed_mask.index.tolist(),
                name=self._observed_mask.index.name,
            ),
            columns=pd.Index(output_labels, name=self._observed_mask.columns.name),
        )
        return export_dataframe(aggregated)

    def _borrow_observed_mask_frame(self) -> pd.DataFrame:
        """Package-private defensive mask snapshot for internal read paths."""

        return borrow_dataframe(self._observed_mask)


def _build_imputation_observation_metadata_or_none(
    *,
    imputation_observation_mask: pd.DataFrame | None,
    phospho_index: pd.Index,
    sample_index: pd.Index,
) -> ImputationObservationMetadata | None:
    if imputation_observation_mask is None:
        return None
    return ImputationObservationMetadata(
        observed_mask=imputation_observation_mask,
        phospho_index=phospho_index,
        sample_index=sample_index,
        _assume_owned=True,
    )


def _require_boolean_observation_mask(mask: pd.DataFrame) -> None:
    values = mask.to_numpy(dtype="object", copy=False)
    missing_values: npt.NDArray[np.bool_] = np.asarray(
        pd.isna(values),
        dtype=bool,
    )
    if bool(missing_values.any()):
        raise DatasetValidationError(
            "dataset.imputation_observation_mask must not contain missing values"
        )
    if all(pd.api.types.is_bool_dtype(dtype) for dtype in mask.dtypes):
        return

    boolean_result = np.frompyfunc(
        lambda value: isinstance(value, (bool, np.bool_)),
        1,
        1,
    )(values)
    boolean_cells: npt.NDArray[np.bool_] = np.asarray(boolean_result, dtype=bool)
    invalid_locations = np.argwhere(~boolean_cells)
    if invalid_locations.size == 0:
        return
    row_index, column_index = invalid_locations[0]
    raise DatasetValidationError(
        "dataset.imputation_observation_mask must contain only boolean "
        "values; "
        f"invalid_cell=({mask.index[int(row_index)]!r}, "
        f"{mask.columns[int(column_index)]!r})"
    )


def _build_imputation_feature_summary(observed_mask: pd.DataFrame) -> pd.DataFrame:
    sample_count = int(observed_mask.shape[1])
    observed_values = observed_mask.to_numpy(dtype=bool)
    observed_counts = observed_values.sum(axis=1).astype(np.int64)
    imputed_counts = (sample_count - observed_counts).astype(np.int64)
    imputed_fraction = imputed_counts.astype(float) / float(sample_count)
    summary = pd.DataFrame(
        {
            "imputed_cell_count": imputed_counts,
            "observed_cell_count": observed_counts,
            "imputed_fraction": imputed_fraction,
        },
        index=observed_mask.index,
    )
    summary.index.name = observed_mask.index.name
    return summary


def _build_imputation_observation_summary(observed_mask: pd.DataFrame) -> pd.DataFrame:
    sample_count = int(observed_mask.shape[1])
    observed_values = observed_mask.to_numpy(dtype=bool)
    observed_counts = observed_values.sum(axis=1).astype(np.int64)
    imputed_counts = (sample_count - observed_counts).astype(np.int64)
    imputed_fraction = imputed_counts.astype(float) / float(sample_count)
    feature_ids = observed_mask.index.tolist()
    summary = pd.DataFrame(
        {
            "feature_id": feature_ids,
            "observed_cell_count": observed_counts,
            "imputed_cell_count": imputed_counts,
            "total_analysed_cell_count": np.full(
                shape=observed_counts.shape,
                fill_value=sample_count,
                dtype=np.int64,
            ),
            "imputed_fraction": imputed_fraction,
        },
        index=pd.Index(feature_ids, name=observed_mask.index.name),
        columns=list(IMPUTATION_OBSERVATION_SUMMARY_COLUMNS),
    )
    return summary


def _requested_label_list(
    labels: Sequence[object],
    *,
    field_name: str,
) -> list[object]:
    requested: list[object]
    if isinstance(labels, str | bytes):
        requested = [labels]
    else:
        try:
            requested = list(labels)
        except TypeError as exc:
            raise DatasetValidationError(
                f"{field_name} must be a sequence of labels"
            ) from exc
    if not requested:
        raise DatasetValidationError(f"{field_name} must contain at least one label")
    return requested


def _label_positions(
    requested_labels: Sequence[object],
    *,
    available_labels: pd.Index,
) -> list[int]:
    positions_by_label = {
        label: int(position) for position, label in enumerate(available_labels.tolist())
    }
    return [positions_by_label[label] for label in requested_labels]


def _require_requested_labels_present(
    requested_labels: Sequence[object],
    *,
    available_labels: pd.Index,
    field_name: str,
    available_field_name: str,
) -> None:
    available = set(available_labels.tolist())
    missing: list[object] = []
    for label in requested_labels:
        if label in available:
            continue
        missing.append(label)
    if missing:
        raise DatasetValidationError(
            f"{field_name} contains labels absent from {available_field_name}: "
            + _format_label_preview(missing)
        )


def _format_label_preview(labels: Sequence[object]) -> str:
    preview = [repr(label) for label in labels[:5]]
    if len(labels) > 5:
        preview.append("...")
    return ", ".join(preview)


def _missing_data_state_claims_no_missing_values(
    processing_state: DatasetProcessingState,
) -> bool:
    missing_data = processing_state.missing_data
    if bool(missing_data.complete_matrix):
        return True
    if missing_data.has_missing_values is False:
        return True
    if missing_data.missing_value_count == 0:
        return True
    diagnostics = missing_data.diagnostics
    if diagnostics is None:
        return False
    return diagnostics.get("output_missing_cell_count") == 0


def _is_missing_value(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value = cast(object, value)
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value = cast(object, value)
        return str(temporal_value) == "NaT"
    return False


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


build_imputation_observation_metadata_or_none = (
    _build_imputation_observation_metadata_or_none
)
is_missing_value = _is_missing_value
missing_data_state_claims_no_missing_values = (
    _missing_data_state_claims_no_missing_values
)
require_boolean_observation_mask = _require_boolean_observation_mask
require_instance = _require_instance
require_optional_instance = _require_optional_instance


__all__ = [
    "ComparisonState",
    "DatasetProcessingState",
    "DatasetPreprocessingReport",
    "IMPUTATION_FEATURE_METADATA_COLUMNS",
    "IMPUTATION_OBSERVATION_SUMMARY_COLUMNS",
    "ImputationObservationMetadata",
    "MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1",
    "MissingDataDiagnostics",
    "MissingDataDiagnosticsV1",
    "MissingDataState",
    "NormalisationState",
    "PreprocessingSiteAttritionSummary",
    "RuvReadinessState",
    "SiteMatrixState",
    "SiteSequenceResolutionReport",
    "SiteSequenceResolutionRowDiagnostic",
    "SiteSequenceResolutionState",
    "TOTAL_PROTEIN_CORRECTION_DIAGNOSTICS_SCHEMA_VERSION_V1",
    "TotalProteinCorrectionDiagnostics",
    "TotalProteinCorrectionDiagnosticsV1",
    "TotalProteinCorrectionState",
    "default_ruv_readiness_state",
    "JsonPrimitive",
    "JsonValue",
]
