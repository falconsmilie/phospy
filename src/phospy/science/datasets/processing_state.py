"""Dataset preprocessing-state summary models."""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TypeVar, cast

import numpy as np
import pandas as pd

from phospy.errors.validation import DatasetValidationError
from phospy.frames.comparison import dataframe_equals, optional_dataframe_equals
from phospy.frames.ownership import (
    borrow_dataframe,
    borrow_optional_dataframe,
    export_dataframe,
    export_optional_dataframe,
    own_dataframe,
    own_optional_dataframe,
)
from phospy.frames.validation import (
    require_columns,
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
from phospy.science.datasets.imputation_metadata import (
    IMPUTATION_FEATURE_METADATA_COLUMNS,
    IMPUTATION_OBSERVATION_SUMMARY_COLUMNS,
    ImputationObservationMetadata,
    _build_imputation_observation_metadata_or_none,
    _require_boolean_observation_mask,
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
from phospy.science.transformations.models import IntensityTransformationEvent

_ExpectedType = TypeVar("_ExpectedType")
_SITE_ID_REASON_PATTERN = re.compile(r"site identifier|site_id|site id")
_INVALID_REASON_PATTERN = re.compile(r"invalid|missing|blank|empty")


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


@dataclass(frozen=True, slots=True, init=False, eq=False)
class DatasetPreprocessingReport:
    """Public provenance report for dataset preprocessing.

    Internal `_borrow_*` accessors return mutation-isolated internal table
    snapshots for trusted internal read paths only.
    """

    __hash__ = object.__hash__

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
    _intensity_transformation_event: IntensityTransformationEvent | None = field(
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
        intensity_transformation_event: IntensityTransformationEvent | None = None,
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
        _require_optional_instance(
            intensity_transformation_event,
            expected_type=IntensityTransformationEvent,
            error_message=(
                "dataset.preprocessing_report.intensity_transformation_event must "
                "be IntensityTransformationEvent or None"
            ),
        )
        object.__setattr__(
            self,
            "_intensity_transformation_event",
            intensity_transformation_event,
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

    @property
    def intensity_transformation_event(
        self,
    ) -> IntensityTransformationEvent | None:
        return self._intensity_transformation_event

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

    def intensity_transformation_summary(
        self,
    ) -> IntensityTransformationEvent | None:
        """Return typed intensity-scale transition evidence when available."""

        return self._intensity_transformation_event

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another report owns the same report content."""

        if not isinstance(other, DatasetPreprocessingReport):
            return False
        return (
            dataframe_equals(self._row_counts, other._row_counts)
            and dataframe_equals(self._operations, other._operations)
            and dataframe_equals(self._row_audit, other._row_audit)
            and optional_dataframe_equals(
                self._duplicate_site_resolution,
                other._duplicate_site_resolution,
            )
            and optional_dataframe_equals(
                self._metadata_conflicts,
                other._metadata_conflicts,
            )
            and optional_dataframe_equals(
                self._comparison_group_stats,
                other._comparison_group_stats,
            )
            and optional_dataframe_equals(
                self._comparison_pair_stats,
                other._comparison_pair_stats,
            )
            and self._site_sequence_resolution == other._site_sequence_resolution
            and self._batch_correction == other._batch_correction
            and _optional_scientific_equals(
                self._protein_aware_preparation,
                other._protein_aware_preparation,
            )
            and self._intensity_transformation_event
            == other._intensity_transformation_event
        )

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
        intensity_transformation_event: IntensityTransformationEvent | None = None,
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
            intensity_transformation_event=intensity_transformation_event,
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
        intensity_transformation_event: IntensityTransformationEvent | None = None,
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
            intensity_transformation_event=intensity_transformation_event,
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


def _optional_scientific_equals(left: object | None, right: object | None) -> bool:
    if left is None or right is None:
        return left is right
    method = getattr(left, "scientifically_equals", None)
    if callable(method):
        return bool(method(right))
    return left == right


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
