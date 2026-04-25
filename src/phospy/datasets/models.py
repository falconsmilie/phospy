"""Dataset domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass

import pandas as pd

from phospy._frame_ownership import own_dataframe, own_optional_dataframe
from phospy.errors.validation import DatasetValidationError
from phospy.references.models import Organism
from phospy.transformations.models import TransformationState
from phospy.validation.datasets.analysis_ready import AnalysisReadyDatasetValidator
from phospy.validation.transformations.state import TransformationStateValidator

_DATASET_VALIDATOR = AnalysisReadyDatasetValidator()
_TRANSFORMATION_STATE_VALIDATOR = TransformationStateValidator()
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
        duplicate_site_resolution: pd.DataFrame | None = None,
        metadata_conflicts: pd.DataFrame | None = None,
        comparison_group_stats: pd.DataFrame | None = None,
        comparison_pair_stats: pd.DataFrame | None = None,
    ) -> DatasetPreprocessingReport:
        return cls(
            row_counts=row_counts,
            operations=operations,
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
    transformation_state: TransformationState
    sample_metadata: pd.DataFrame | None = None
    total: pd.DataFrame | None = None
    comparisons: pd.DataFrame | None = None
    organism: Organism | None = None
    preprocessing_report: DatasetPreprocessingReport | None = None
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
        _DATASET_VALIDATOR.run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            organism=self.organism,
        )
        transformation_state = _TRANSFORMATION_STATE_VALIDATOR.run(
            transformation_state=self.transformation_state,
            has_total_matrix=total is not None,
            require_established=True,
        )
        if self.preprocessing_report is not None and not isinstance(
            self.preprocessing_report,
            DatasetPreprocessingReport,
        ):
            raise DatasetValidationError(
                "dataset.preprocessing_report must be DatasetPreprocessingReport "
                "or None"
            )
        object.__setattr__(self, "phospho", phospho)
        object.__setattr__(self, "site_metadata", site_metadata)
        object.__setattr__(self, "sample_metadata", sample_metadata)
        object.__setattr__(self, "total", total)
        object.__setattr__(self, "comparisons", comparisons)
        object.__setattr__(self, "transformation_state", transformation_state)

    @classmethod
    def _from_owned(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        transformation_state: TransformationState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        organism: Organism | None = None,
        preprocessing_report: DatasetPreprocessingReport | None = None,
    ) -> AnalysisReadyPhosphoDataset:
        return cls(
            phospho=phospho,
            site_metadata=site_metadata,
            transformation_state=transformation_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            organism=organism,
            preprocessing_report=preprocessing_report,
            _assume_owned=True,
        )
