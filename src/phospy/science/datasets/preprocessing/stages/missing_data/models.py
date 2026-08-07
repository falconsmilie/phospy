"""Typed policy outcomes and shared missing-data stage inputs."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.science.datasets.processing_state import JsonValue


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
