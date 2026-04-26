"""Dataset domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import own_dataframe, own_optional_dataframe
from phospy.datasets.processing_state import DatasetProcessingState
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
    require_dataframe,
    require_exact_index_match,
    require_numeric_dataframe,
    require_unique_columns,
)
from phospy.validation.common.missing_values import (
    MissingValuePolicy,
    require_missing_value_policy,
)
from phospy.validation.transformations.state import IntensityScaleStateValidator

_INTENSITY_SCALE_STATE_VALIDATOR = IntensityScaleStateValidator()
_PREPROCESSING_REPORT_ROW_COUNT_COLUMNS = (
    "stage",
    "input_rows",
    "output_rows",
    "dropped_rows",
)
_PREPROCESSING_REPORT_OPERATION_COLUMNS = (
    "step_order",
    "stage",
    "operation",
    "parameters",
    "input_rows",
    "output_rows",
    "notes",
)
_PREPROCESSING_REPORT_ROW_AUDIT_COLUMNS = (
    "stage",
    "action",
    "reason",
    "source_row_id",
    "site_id",
    "retained",
    "retained_row_id",
    "source_rows",
    "retained_row",
    "parameter_snapshot",
)
_PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_COLUMNS = (
    "site_id",
    "source_row_id",
    "retained",
    "resolution_strategy",
    "retained_reason",
    "dropped_reason",
    "observed_values",
    "mean_signal",
    "n_source_rows",
    "n_aggregated_rows",
    "source_protein_id",
    "source_gene_symbol",
    "source_site",
    "source_site_sequence",
    "metadata_conflict_detected",
)
_PREPROCESSING_REPORT_METADATA_CONFLICT_COLUMNS = (
    "site_id",
    "field",
    "values",
    "n_distinct_values",
    "source_row_ids",
)
_PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_COLUMNS = (
    "site_id",
    "group",
    "n",
    "mean",
    "sd",
    "sem",
)
_PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_COLUMNS = (
    "site_id",
    "comparison",
    "left_group",
    "right_group",
    "left_n",
    "right_n",
    "left_mean",
    "right_mean",
    "left_sd",
    "right_sd",
    "left_sem",
    "right_sem",
    "effect_size",
)


@dataclass(frozen=True, slots=True)
class DatasetPreprocessingReport:
    """Public provenance report for dataset preprocessing."""

    row_counts: pd.DataFrame
    operations: pd.DataFrame
    row_audit: pd.DataFrame
    duplicate_site_resolution: pd.DataFrame | None = None
    metadata_conflicts: pd.DataFrame | None = None
    comparison_group_stats: pd.DataFrame | None = None
    comparison_pair_stats: pd.DataFrame | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        row_counts = own_dataframe(
            self.row_counts,
            field_name="dataset.preprocessing_report.row_counts",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        operations = own_dataframe(
            self.operations,
            field_name="dataset.preprocessing_report.operations",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        row_audit = own_dataframe(
            self.row_audit,
            field_name="dataset.preprocessing_report.row_audit",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        duplicate_site_resolution = own_optional_dataframe(
            self.duplicate_site_resolution,
            field_name="dataset.preprocessing_report.duplicate_site_resolution",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        metadata_conflicts = own_optional_dataframe(
            self.metadata_conflicts,
            field_name="dataset.preprocessing_report.metadata_conflicts",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        comparison_group_stats = own_optional_dataframe(
            self.comparison_group_stats,
            field_name="dataset.preprocessing_report.comparison_group_stats",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        comparison_pair_stats = own_optional_dataframe(
            self.comparison_pair_stats,
            field_name="dataset.preprocessing_report.comparison_pair_stats",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        missing_row_columns = [
            column
            for column in _PREPROCESSING_REPORT_ROW_COUNT_COLUMNS
            if column not in row_counts.columns
        ]
        if missing_row_columns:
            missing = ", ".join(missing_row_columns)
            raise DatasetValidationError(
                "dataset.preprocessing_report.row_counts is missing required "
                f"columns: {missing}"
            )
        missing_operation_columns = [
            column
            for column in _PREPROCESSING_REPORT_OPERATION_COLUMNS
            if column not in operations.columns
        ]
        if missing_operation_columns:
            missing = ", ".join(missing_operation_columns)
            raise DatasetValidationError(
                "dataset.preprocessing_report.operations is missing required "
                f"columns: {missing}"
            )
        missing_row_audit_columns = [
            column
            for column in _PREPROCESSING_REPORT_ROW_AUDIT_COLUMNS
            if column not in row_audit.columns
        ]
        if missing_row_audit_columns:
            missing = ", ".join(missing_row_audit_columns)
            raise DatasetValidationError(
                "dataset.preprocessing_report.row_audit is missing required "
                f"columns: {missing}"
            )
        if duplicate_site_resolution is not None:
            missing_duplicate_columns = [
                column
                for column in _PREPROCESSING_REPORT_DUPLICATE_SITE_RESOLUTION_COLUMNS
                if column not in duplicate_site_resolution.columns
            ]
            if missing_duplicate_columns:
                missing = ", ".join(missing_duplicate_columns)
                raise DatasetValidationError(
                    "dataset.preprocessing_report.duplicate_site_resolution is "
                    f"missing required columns: {missing}"
                )
        if metadata_conflicts is not None:
            missing_conflict_columns = [
                column
                for column in _PREPROCESSING_REPORT_METADATA_CONFLICT_COLUMNS
                if column not in metadata_conflicts.columns
            ]
            if missing_conflict_columns:
                missing = ", ".join(missing_conflict_columns)
                raise DatasetValidationError(
                    "dataset.preprocessing_report.metadata_conflicts is missing "
                    f"required columns: {missing}"
                )
        if comparison_group_stats is not None:
            missing_group_stats_columns = [
                column
                for column in _PREPROCESSING_REPORT_COMPARISON_GROUP_STATS_COLUMNS
                if column not in comparison_group_stats.columns
            ]
            if missing_group_stats_columns:
                missing = ", ".join(missing_group_stats_columns)
                raise DatasetValidationError(
                    "dataset.preprocessing_report.comparison_group_stats is "
                    f"missing required columns: {missing}"
                )
        if comparison_pair_stats is not None:
            missing_pair_stats_columns = [
                column
                for column in _PREPROCESSING_REPORT_COMPARISON_PAIR_STATS_COLUMNS
                if column not in comparison_pair_stats.columns
            ]
            if missing_pair_stats_columns:
                missing = ", ".join(missing_pair_stats_columns)
                raise DatasetValidationError(
                    "dataset.preprocessing_report.comparison_pair_stats is "
                    f"missing required columns: {missing}"
                )
        object.__setattr__(self, "row_counts", row_counts)
        object.__setattr__(self, "operations", operations)
        object.__setattr__(self, "row_audit", row_audit)
        object.__setattr__(self, "duplicate_site_resolution", duplicate_site_resolution)
        object.__setattr__(self, "metadata_conflicts", metadata_conflicts)
        object.__setattr__(self, "comparison_group_stats", comparison_group_stats)
        object.__setattr__(self, "comparison_pair_stats", comparison_pair_stats)

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


@dataclass(frozen=True, slots=True)
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
    """

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    intensity_scale_state: IntensityScaleState
    processing_state: DatasetProcessingState
    sample_metadata: pd.DataFrame | None = None
    total: pd.DataFrame | None = None
    comparisons: pd.DataFrame | None = None
    organism: Organism | None = None
    preprocessing_report: DatasetPreprocessingReport | None = None
    provenance: RunProvenance | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        phospho = own_dataframe(
            self.phospho,
            field_name="dataset.phospho",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        site_metadata = own_dataframe(
            self.site_metadata,
            field_name="dataset.site_metadata",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        sample_metadata = own_optional_dataframe(
            self.sample_metadata,
            field_name="dataset.sample_metadata",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        total = own_optional_dataframe(
            self.total,
            field_name="dataset.total",
            error_type=DatasetValidationError,
            assume_owned=_assume_owned,
        )
        comparisons = own_optional_dataframe(
            self.comparisons,
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
                allow_empty=False,
                error_type=DatasetValidationError,
            )
            require_numeric_dataframe(
                comparisons_frame,
                field_name="dataset.comparisons",
                error_type=DatasetValidationError,
            )
            require_missing_value_policy(
                comparisons_frame,
                field_name="dataset.comparisons",
                policy=MissingValuePolicy.FORBID,
                error_type=DatasetValidationError,
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
        intensity_scale_state = _INTENSITY_SCALE_STATE_VALIDATOR.run(
            intensity_scale_state=self.intensity_scale_state,
            has_total_matrix=total_table is not None,
            require_established=True,
        )
        if not isinstance(self.processing_state, DatasetProcessingState):
            raise DatasetValidationError(
                "dataset.processing_state must be a DatasetProcessingState instance"
            )
        if self.processing_state.intensity_scale != intensity_scale_state:
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
        object.__setattr__(self, "phospho", phospho_table.frame)
        object.__setattr__(self, "site_metadata", site_metadata_table.frame)
        object.__setattr__(
            self,
            "sample_metadata",
            None if sample_metadata_table is None else sample_metadata_table.frame,
        )
        object.__setattr__(
            self, "total", None if total_table is None else total_table.frame
        )
        object.__setattr__(self, "comparisons", comparisons)
        object.__setattr__(self, "intensity_scale_state", intensity_scale_state)
        object.__setattr__(self, "processing_state", self.processing_state)

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
