"""Typed policy outcomes and shared missing-data stage inputs."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum

import pandas as pd

from phospy.science.datasets.preprocessing.sample_group_metadata import (
    ResolvedSampleGroups,
)
from phospy.science.datasets.processing_state import JsonValue


class GroupMissingnessClassification(str, Enum):
    """Observed-pattern classification for one phosphosite/group block."""

    COMPLETE = "complete"
    SUPPORTED_PARTIAL = "supported_partial"
    UNSUPPORTED_PARTIAL = "unsupported_partial"
    SUPPORTED_FULLY_MISSING = "supported_fully_missing"
    UNSUPPORTED_FULLY_MISSING = "unsupported_fully_missing"


class GroupMissingnessRoute(str, Enum):
    """Allowed numerical route assigned from original missingness only."""

    NONE = "none"
    KNN = "knn"
    MINPROB = "minprob"


class GroupRoutingAssumption(str, Enum):
    """Scientific assumption attached to an eligible routing decision."""

    NONE = "none"
    SIMILARITY_BASED_ELIGIBILITY = "similarity_based_eligibility"
    ASYMMETRIC_LEFT_CENSORED = "asymmetric_left_censored"


@dataclass(frozen=True, slots=True)
class GroupRoutingFact:
    """Compact classification facts for one phosphosite/group block."""

    row_id: str
    group_label: str
    group_sample_count: int
    observed_finite_count: int
    missing_count: int
    observed_fraction: float
    classification: GroupMissingnessClassification
    route: GroupMissingnessRoute
    routing_assumption: GroupRoutingAssumption
    qualifying_reference_group_labels: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class DroppedRowRoutingRecord:
    """Unsupported group states that make one row ineligible for retention."""

    row_id: str
    reasons_by_group: Mapping[str, GroupMissingnessClassification]


@dataclass(frozen=True, slots=True)
class GroupAwareRoutingOutcome:
    """Typed, pre-imputation result of group-aware missingness routing."""

    retained_row_mask: pd.Series
    retained_row_ids: tuple[str, ...]
    dropped_row_ids: tuple[str, ...]
    dropped_row_reasons: tuple[DroppedRowRoutingRecord, ...]
    resolved_groups: ResolvedSampleGroups
    group_facts: tuple[GroupRoutingFact, ...]
    knn_target_mask: pd.DataFrame
    minprob_target_mask: pd.DataFrame
    knn_target_cell_count: int
    minprob_target_cell_count: int
    unsupported_group_count: int
    unsupported_partial_group_count: int
    unsupported_fully_missing_group_count: int
    min_partial_observed_fraction: float
    min_reference_observed_fraction: float
    original_missingness_mask_hash: str

    def retain_rows(self, frame: pd.DataFrame) -> pd.DataFrame:
        """Return a copy containing only rows eligible for numerical imputation."""

        if not frame.index.equals(self.retained_row_mask.index):
            raise ValueError(
                "group-aware routing can retain rows only from a frame aligned "
                "to the classified phospho index"
            )
        return frame.loc[self.retained_row_mask].copy(deep=True)


@dataclass(frozen=True, slots=True)
class MissingDataInputProfile:
    """Input missingness profile computed once by the stage coordinator."""

    input_missing_cell_count: int
    affected_row_ids: tuple[str, ...]
    affected_column_ids: tuple[str, ...]
    missingness_mask_hash: str

    @property
    def affected_row_count(self) -> int:
        return len(self.affected_row_ids)

    @property
    def affected_column_count(self) -> int:
        return len(self.affected_column_ids)


@dataclass(frozen=True, slots=True)
class RowImputationRecord:
    """Per-row imputation summary used to build row-audit records."""

    row_id: str
    imputed_columns: tuple[str, ...]
    imputed_cell_count: int
    nearest_neighbour_imputed_columns: tuple[str, ...] = ()
    column_mean_fallback_columns: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class RowMedianPolicyOutcome:
    """Numerical output for row-median missing-data policy."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    imputed_mask: pd.DataFrame
    min_observed_values: int
    dropped_row_ids: tuple[str, ...]
    dropped_row_observed_values: tuple[tuple[str, int], ...]
    imputed_cell_count: int
    imputed_row_ids: tuple[str, ...]
    imputed_column_ids: tuple[str, ...]
    output_missing_cell_count: int
    rows_not_imputable: tuple[str, ...]
    row_medians_used: dict[str, float]
    imputed_rows: tuple[RowImputationRecord, ...]


@dataclass(frozen=True, slots=True)
class KnnPolicyOutcome:
    """Numerical output for KNN missing-data policy."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    imputed_mask: pd.DataFrame
    nearest_neighbour_imputed_mask: pd.DataFrame
    column_mean_fallback_imputed_mask: pd.DataFrame
    k: int
    distance: str
    max_missing_fraction_per_row: float
    no_overlap_policy: str
    no_overlap_policy_version: int
    dropped_row_ids: tuple[str, ...]
    dropped_rows_missing_fraction: tuple[tuple[str, float], ...]
    imputed_cell_count: int
    imputed_row_ids: tuple[str, ...]
    imputed_column_ids: tuple[str, ...]
    nearest_neighbour_imputed_cell_count: int
    nearest_neighbour_imputed_row_ids: tuple[str, ...]
    nearest_neighbour_imputed_column_ids: tuple[str, ...]
    column_mean_fallback_imputed_cell_count: int
    column_mean_fallback_row_ids: tuple[str, ...]
    column_mean_fallback_column_ids: tuple[str, ...]
    nearest_neighbour_imputation_mask_hash: str
    column_mean_fallback_imputation_mask_hash: str
    fully_column_mean_fallback_row_ids: tuple[str, ...]
    output_missing_cell_count: int
    rows_not_imputable: tuple[str, ...]
    imputed_rows: tuple[RowImputationRecord, ...]


@dataclass(frozen=True, slots=True)
class MinProbPolicyOutcome:
    """Numerical output for minprob missing-data policy."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    imputed_mask: pd.DataFrame
    q: float
    width: float
    seed: int
    max_missing_fraction_per_row: float
    dropped_row_ids: tuple[str, ...]
    dropped_rows_missing_fraction: tuple[tuple[str, float], ...]
    imputed_cell_count: int
    imputed_row_ids: tuple[str, ...]
    imputed_column_ids: tuple[str, ...]
    output_missing_cell_count: int
    rows_not_imputable: tuple[str, ...]
    per_column_distribution_parameters: dict[str, dict[str, JsonValue]]
    imputed_rows: tuple[RowImputationRecord, ...]


@dataclass(frozen=True, slots=True)
class GroupAwarePolicyOutcome:
    """Merged numerical output for group-aware mixed-mechanism imputation."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    imputed_mask: pd.DataFrame
    knn_target_mask: pd.DataFrame
    minprob_target_mask: pd.DataFrame
    knn_imputed_mask: pd.DataFrame
    minprob_imputed_mask: pd.DataFrame
    routing: GroupAwareRoutingOutcome
    q: float
    width: float
    seed: int
    k: int
    distance: str
    no_overlap_policy: str
    per_column_distribution_parameters: dict[str, dict[str, JsonValue]]
    dropped_row_ids: tuple[str, ...]
    imputed_cell_count: int
    imputed_row_ids: tuple[str, ...]
    imputed_column_ids: tuple[str, ...]
    output_missing_cell_count: int
    rows_not_imputable: tuple[str, ...]
    imputed_rows: tuple[RowImputationRecord, ...]
