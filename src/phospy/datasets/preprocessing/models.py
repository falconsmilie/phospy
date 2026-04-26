"""Internal dataset preprocessing planning and state models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Protocol

import pandas as pd

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN,
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY,
    DATASET_MISSING_DATA_POLICY_FORBID,
    DATASET_NORMALISATION_POLICY_NONE,
    DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    DATASET_SITE_MATRIX_POLICY_AS_INPUT,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
    DatasetComparisonBuildingPolicy,
    DatasetComparisonPair,
    DatasetIntensityTransformPolicy,
    DatasetMissingDataPolicy,
    DatasetNormalisationPolicy,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixDuplicateSitePolicy,
    DatasetSiteMatrixMissingDataPolicy,
    DatasetSiteMatrixPolicy,
    DatasetTotalProteinCorrectionPolicy,
)
from phospy.datasets.preprocessing.report_schema import (
    ROW_AUDIT_COLUMNS,
    PreprocessingRowAuditRow,
    dataframe_from_row_audit_rows,
    reorder_columns,
)

DATASET_PREPROCESSING_STAGE_MISSING_DATA = "missing_data"
DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION = "total_protein_correction"
DATASET_PREPROCESSING_STAGE_SITE_MATRIX = "site_matrix"
DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM = "intensity_transform"
DATASET_PREPROCESSING_STAGE_NORMALISATION = "normalisation"
DATASET_PREPROCESSING_STAGE_COMPARISONS = "comparisons"
DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT = (DATASET_PREPROCESSING_STAGE_MISSING_DATA,)


@dataclass(frozen=True, slots=True)
class PreprocessingPlan:
    """Execution-ready internal preprocessing plan derived from public config."""

    intensity_transform_policy: DatasetIntensityTransformPolicy = (
        DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY
    )
    intensity_transform_pseudocount: float = 1.0
    normalisation_policy: DatasetNormalisationPolicy = DATASET_NORMALISATION_POLICY_NONE
    missing_data_policy: DatasetMissingDataPolicy = DATASET_MISSING_DATA_POLICY_FORBID
    missing_data_min_observed_values: int | None = None
    total_protein_correction_policy: DatasetTotalProteinCorrectionPolicy = (
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )
    site_matrix_policy: DatasetSiteMatrixPolicy = DATASET_SITE_MATRIX_POLICY_AS_INPUT
    comparison_building_policy: DatasetComparisonBuildingPolicy = (
        DATASET_COMPARISON_BUILDING_POLICY_NONE
    )
    site_matrix_duplicate_site_policy: DatasetSiteMatrixDuplicateSitePolicy = (
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL
    )
    site_matrix_missing_data_policy: DatasetSiteMatrixMissingDataPolicy = (
        DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
    )
    site_matrix_minimum_observed_values: int | None = None
    comparison_sample_group_column: str = (
        DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN
    )
    comparison_pairs: tuple[DatasetComparisonPair, ...] | None = None
    stage_order: tuple[str, ...] = DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT

    @classmethod
    def from_config(cls, config: DatasetPreprocessingConfig) -> PreprocessingPlan:
        stage_order: list[str] = [DATASET_PREPROCESSING_STAGE_MISSING_DATA]
        if (
            config.intensity_transform.policy
            != DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY
        ):
            stage_order.append(DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM)
        if (
            config.total_protein_correction.policy
            != DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
        ):
            stage_order.append(DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION)
        if config.site_matrix.policy != DATASET_SITE_MATRIX_POLICY_AS_INPUT:
            stage_order.append(DATASET_PREPROCESSING_STAGE_SITE_MATRIX)
        if config.normalisation.policy != DATASET_NORMALISATION_POLICY_NONE:
            stage_order.append(DATASET_PREPROCESSING_STAGE_NORMALISATION)
        if config.comparisons.policy != DATASET_COMPARISON_BUILDING_POLICY_NONE:
            stage_order.append(DATASET_PREPROCESSING_STAGE_COMPARISONS)
        return cls(
            intensity_transform_policy=config.intensity_transform.policy,
            intensity_transform_pseudocount=float(
                config.intensity_transform.pseudocount
            ),
            normalisation_policy=config.normalisation.policy,
            missing_data_policy=config.missing_data.policy,
            missing_data_min_observed_values=config.missing_data.min_observed_values,
            total_protein_correction_policy=config.total_protein_correction.policy,
            site_matrix_policy=config.site_matrix.policy,
            site_matrix_duplicate_site_policy=config.site_matrix.duplicate_site_policy,
            site_matrix_missing_data_policy=config.site_matrix.missing_data_policy,
            site_matrix_minimum_observed_values=config.site_matrix.minimum_observed_values,
            comparison_building_policy=config.comparisons.policy,
            comparison_sample_group_column=config.comparisons.sample_group_column,
            comparison_pairs=(
                None
                if config.comparisons.pairs is None
                else tuple(config.comparisons.pairs)
            ),
            stage_order=tuple(stage_order),
        )

    @classmethod
    def default(cls) -> PreprocessingPlan:
        return cls.from_config(DatasetPreprocessingConfig())


@dataclass(frozen=True, slots=True)
class PreprocessingState:
    """Internal preprocessing state carried between ordered stages."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    plan: PreprocessingPlan
    comparisons: pd.DataFrame | None = None
    comparison_group_stats: pd.DataFrame | None = None
    comparison_pair_stats: pd.DataFrame | None = None
    duplicate_site_resolution: pd.DataFrame | None = None
    metadata_conflicts: pd.DataFrame | None = None
    row_audit: pd.DataFrame | None = None
    report_rows: tuple[PreprocessingReportRow, ...] = ()


def empty_preprocessing_row_audit() -> pd.DataFrame:
    """Return an empty stable-schema preprocessing row-audit table."""

    return dataframe_from_row_audit_rows(())


def append_row_audit_records(
    state: PreprocessingState,
    records: Sequence[PreprocessingRowAuditRow],
) -> PreprocessingState:
    """Append row-audit records without mutating existing state frames."""

    if not records:
        return state

    existing = (
        empty_preprocessing_row_audit()
        if state.row_audit is None
        else reorder_columns(
            state.row_audit,
            expected_columns=ROW_AUDIT_COLUMNS,
        ).copy(deep=True)
    )
    appended = dataframe_from_row_audit_rows(records)
    combined = pd.concat([existing, appended], axis=0, ignore_index=True)
    return replace(state, row_audit=combined)


@dataclass(frozen=True, slots=True)
class ComparisonBuildResult:
    """Structured comparison-building output with provenance sidecars."""

    comparisons: pd.DataFrame
    comparison_group_stats: pd.DataFrame
    comparison_pair_stats: pd.DataFrame


@dataclass(frozen=True, slots=True)
class DuplicateSiteResolutionResult:
    """Duplicate-site policy output with structured provenance tables."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    dropped_row_count: int
    duplicate_site_resolution: pd.DataFrame
    metadata_conflicts: pd.DataFrame


@dataclass(frozen=True, slots=True)
class PreprocessingReportRow:
    """Structured stage-owned contribution to preprocessing report assembly."""

    table: str
    values: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class PreprocessingStageResult:
    """Structured output for a single preprocessing stage execution."""

    state: PreprocessingState
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    report_rows: Sequence[PreprocessingReportRow] = ()


@dataclass(frozen=True, slots=True)
class PreprocessingStageExecution:
    """Executed preprocessing stage provenance trace."""

    stage: str
    operation: str
    parameters: dict[str, object]
    input_shape: tuple[int, int]
    output_shape: tuple[int, int]
    input_hash: str
    output_hash: str
    dropped_row_ids: tuple[str, ...] = ()
    dropped_row_count: int = 0
    imputed_cell_count: int = 0
    imputed_row_ids: tuple[str, ...] = ()
    notes: str | None = None
    diagnostics: dict[str, object] = field(default_factory=dict)

    @property
    def input_rows(self) -> int:
        return int(self.input_shape[0])

    @property
    def output_rows(self) -> int:
        return int(self.output_shape[0])


class PreprocessingStage(Protocol):
    """Single internal preprocessing stage contract."""

    stage_key: str

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        """Apply a preprocessing stage and return its structured result."""


__all__ = [
    "DATASET_PREPROCESSING_STAGE_COMPARISONS",
    "DATASET_PREPROCESSING_STAGE_INTENSITY_TRANSFORM",
    "DATASET_PREPROCESSING_STAGE_MISSING_DATA",
    "DATASET_PREPROCESSING_STAGE_NORMALISATION",
    "DATASET_PREPROCESSING_STAGE_ORDER_DEFAULT",
    "DATASET_PREPROCESSING_STAGE_SITE_MATRIX",
    "DATASET_PREPROCESSING_STAGE_TOTAL_PROTEIN_CORRECTION",
    "ComparisonBuildResult",
    "DuplicateSiteResolutionResult",
    "append_row_audit_records",
    "empty_preprocessing_row_audit",
    "PreprocessingPlan",
    "PreprocessingReportRow",
    "PreprocessingStageResult",
    "PreprocessingStageExecution",
    "PreprocessingStage",
    "PreprocessingState",
]
