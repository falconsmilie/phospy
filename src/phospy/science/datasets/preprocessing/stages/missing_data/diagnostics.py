"""Diagnostics and provenance helpers for missing-data stage."""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from typing import Protocol

import numpy as np
import pandas as pd

from phospy.errors.provenance import ProvenanceFingerprintError
from phospy.provenance.hashing import (
    fingerprint_table_normalized_axes,
    hash_table_tolerance,
)
from phospy.science.datasets.processing_state import (
    MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
    JsonValue,
)

from .models import MissingDataInputProfile


class _HashUpdater(Protocol):
    def update(self, data: bytes, /) -> object: ...


def build_input_profile(phospho: pd.DataFrame) -> MissingDataInputProfile:
    """Return missingness summary for stage-level diagnostics and audit."""

    input_missing_mask = phospho.isna()
    return MissingDataInputProfile(
        input_missing_cell_count=int(input_missing_mask.to_numpy().sum()),
        affected_row_ids=tuple(
            str(row_id) for row_id in phospho.index[input_missing_mask.any(axis=1)]
        ),
        affected_column_ids=tuple(
            str(column_id)
            for column_id in phospho.columns[input_missing_mask.any(axis=0)]
        ),
        missingness_mask_hash=hash_missingness_mask(input_missing_mask),
    )


def hash_missingness_mask(mask: pd.DataFrame) -> str:
    """Return stable fingerprint for input missingness structure."""

    return hash_table_tolerance(
        mask.astype("int8"),
        name="missing_data.input_missingness_mask",
    )


def hash_imputation_mask(mask: pd.DataFrame) -> str:
    """Return stable fingerprint for policy-owned imputed-cell mask."""

    return hash_table_tolerance(
        mask.astype("int8"),
        name="missing_data.imputation_mask",
    )


def hash_knn_mechanism_mask(mask: pd.DataFrame, *, mechanism_name: str) -> str:
    """Return row/column-order-stable fingerprint for KNN mechanism masks."""

    hash_name = f"missing_data.knn_{mechanism_name}_imputation_mask"
    fast_hash = _hash_unique_string_labeled_boolean_mask_normalized_axes(
        mask,
        name=hash_name,
    )
    if fast_hash is not None:
        return fast_hash

    normalized_mask = mask.astype("int8")
    try:
        return fingerprint_table_normalized_axes(
            normalized_mask,
            name=hash_name,
        ).tolerance_hash_value
    except ProvenanceFingerprintError:
        return hash_table_tolerance(normalized_mask, name=hash_name)


def _hash_unique_string_labeled_boolean_mask_normalized_axes(
    mask: pd.DataFrame,
    *,
    name: str,
) -> str | None:
    row_labels = _unique_string_labels(mask.index)
    column_labels = _unique_string_labels(mask.columns)
    if row_labels is None or column_labels is None:
        return None

    row_positions = np.asarray(
        sorted(range(len(row_labels)), key=row_labels.__getitem__),
        dtype=np.intp,
    )
    column_positions = np.asarray(
        sorted(range(len(column_labels)), key=column_labels.__getitem__),
        dtype=np.intp,
    )
    values = mask.to_numpy(dtype=bool, copy=False)
    normalized_values = values[np.ix_(row_positions, column_positions)]

    digest = hashlib.sha256()
    digest.update(b"phospy:boolean-mask-normalized-axes:v1\n")
    _update_hash_string(digest, name)
    _update_hash_string_axis(
        digest,
        (row_labels[int(position)] for position in row_positions),
    )
    _update_hash_string_axis(
        digest,
        (column_labels[int(position)] for position in column_positions),
    )
    digest.update(str(int(normalized_values.shape[0])).encode("ascii"))
    digest.update(b"x")
    digest.update(str(int(normalized_values.shape[1])).encode("ascii"))
    digest.update(b"\n")
    digest.update(np.packbits(normalized_values.ravel(order="C")).tobytes())
    return digest.hexdigest()


def _unique_string_labels(index: pd.Index) -> list[str] | None:
    labels = index.tolist()
    if not all(isinstance(label, str) for label in labels):
        return None
    if len(set(labels)) != len(labels):
        return None
    return labels


def _update_hash_string_axis(
    digest: _HashUpdater,
    labels: Iterable[str],
) -> None:
    materialized_labels = tuple(labels)
    digest.update(str(len(materialized_labels)).encode("ascii"))
    digest.update(b"\n")
    for label in materialized_labels:
        _update_hash_string(digest, label)


def _update_hash_string(digest: _HashUpdater, value: str) -> None:
    encoded = value.encode("utf-8")
    digest.update(str(len(encoded)).encode("ascii"))
    digest.update(b":")
    digest.update(encoded)
    digest.update(b"\n")


def label_preview(values: list[object], *, max_items: int = 3) -> str:
    """Render a bounded preview of labels for human-facing error messages."""

    if not values:
        return "none"
    rendered = [repr(str(value)) for value in values[:max_items]]
    remaining_count = len(values) - len(rendered)
    if remaining_count > 0:
        rendered.append(f"+{remaining_count} more")
    return ", ".join(rendered)


def build_missing_data_diagnostics(
    *,
    missing_data_policy: str,
    imputation_method_id: str | None,
    imputation_method_family: str | None,
    input_missing_cell_count: int,
    output_missing_cell_count: int,
    imputed_cell_count: int,
    affected_row_ids: tuple[str, ...],
    affected_column_ids: tuple[str, ...],
    imputed_row_ids: tuple[str, ...],
    imputed_column_ids: tuple[str, ...],
    dropped_row_ids: tuple[str, ...],
    random_seed: int | None,
    method_parameters: dict[str, JsonValue],
    matrix_scale_requirement: str | None,
    imputation_input_scale: str | None,
    imputation_input_scale_source: str | None,
    imputation_operation_order: str | None,
    stage_order: tuple[str, ...],
    missingness_mask_hash: str,
    left_censored_assumption: bool | None,
    rows_not_imputable: tuple[str, ...],
    row_medians_used: dict[str, float],
    per_column_distribution_parameters: dict[str, dict[str, JsonValue]] | None,
    dropped_rows_above_max_missing_fraction: tuple[str, ...],
    neighbour_count: int | None,
    distance_metric: str | None,
    imputation_mask_hash: str | None,
    knn_no_overlap_policy: str | None = None,
    knn_no_overlap_policy_version: int | None = None,
    knn_nearest_neighbour_imputed_cell_count: int | None = None,
    knn_nearest_neighbour_imputed_row_ids: tuple[str, ...] = (),
    knn_nearest_neighbour_imputed_column_ids: tuple[str, ...] = (),
    knn_column_mean_fallback_imputed_cell_count: int | None = None,
    knn_column_mean_fallback_row_ids: tuple[str, ...] = (),
    knn_column_mean_fallback_column_ids: tuple[str, ...] = (),
    knn_nearest_neighbour_imputation_mask_hash: str | None = None,
    knn_column_mean_fallback_imputation_mask_hash: str | None = None,
    knn_fully_column_mean_fallback_row_ids: tuple[str, ...] = (),
    diagnostic_caveat_codes: tuple[str, ...] = (),
) -> dict[str, JsonValue]:
    diagnostics: dict[str, JsonValue] = {
        "diagnostics_schema_version": MISSING_DATA_DIAGNOSTICS_SCHEMA_VERSION_V1,
        "missing_data_policy": missing_data_policy,
        "imputation_method_id": imputation_method_id,
        "imputation_method_family": imputation_method_family,
        "input_missing_cell_count": int(input_missing_cell_count),
        "output_missing_cell_count": int(output_missing_cell_count),
        "imputed_cell_count": int(imputed_cell_count),
        "affected_row_count": int(len(affected_row_ids)),
        "affected_column_count": int(len(affected_column_ids)),
        "affected_row_ids": list(affected_row_ids),
        "affected_column_ids": list(affected_column_ids),
        "imputed_row_ids": list(imputed_row_ids),
        "imputed_column_ids": list(imputed_column_ids),
        "dropped_row_ids": list(dropped_row_ids),
        "imputed_row_count": int(len(imputed_row_ids)),
        "imputed_column_count": int(len(imputed_column_ids)),
        "dropped_row_count": int(len(dropped_row_ids)),
        "random_seed": random_seed,
        "method_parameters": dict(method_parameters),
        "matrix_scale_requirement": matrix_scale_requirement,
        "imputation_input_scale": imputation_input_scale,
        "imputation_input_scale_source": imputation_input_scale_source,
        "imputation_operation_order": imputation_operation_order,
        "stage_order": list(stage_order),
        "missingness_mask_hash": missingness_mask_hash,
        "left_censored_assumption": left_censored_assumption,
        "rows_not_imputable": list(rows_not_imputable),
        "row_medians_used": {
            str(row_id): float(row_median)
            for row_id, row_median in row_medians_used.items()
        },
        "dropped_rows_above_max_missing_fraction": list(
            dropped_rows_above_max_missing_fraction
        ),
        "neighbour_count": neighbour_count,
        "distance_metric": distance_metric,
    }
    if imputation_mask_hash is not None:
        diagnostics["imputation_mask_hash"] = str(imputation_mask_hash)
    if knn_no_overlap_policy is not None:
        diagnostics["knn_no_overlap_policy"] = str(knn_no_overlap_policy)
    if knn_no_overlap_policy_version is not None:
        diagnostics["knn_no_overlap_policy_version"] = int(
            knn_no_overlap_policy_version
        )
    if knn_nearest_neighbour_imputed_cell_count is not None:
        diagnostics["knn_nearest_neighbour_imputed_cell_count"] = int(
            knn_nearest_neighbour_imputed_cell_count
        )
    if knn_column_mean_fallback_imputed_cell_count is not None:
        diagnostics["knn_column_mean_fallback_imputed_cell_count"] = int(
            knn_column_mean_fallback_imputed_cell_count
        )
    if knn_nearest_neighbour_imputed_row_ids:
        diagnostics["knn_nearest_neighbour_imputed_row_ids"] = list(
            knn_nearest_neighbour_imputed_row_ids
        )
    if knn_nearest_neighbour_imputed_column_ids:
        diagnostics["knn_nearest_neighbour_imputed_column_ids"] = list(
            knn_nearest_neighbour_imputed_column_ids
        )
    if knn_column_mean_fallback_row_ids:
        diagnostics["knn_column_mean_fallback_row_ids"] = list(
            knn_column_mean_fallback_row_ids
        )
    if knn_column_mean_fallback_column_ids:
        diagnostics["knn_column_mean_fallback_column_ids"] = list(
            knn_column_mean_fallback_column_ids
        )
    if knn_nearest_neighbour_imputation_mask_hash is not None:
        diagnostics["knn_nearest_neighbour_imputation_mask_hash"] = str(
            knn_nearest_neighbour_imputation_mask_hash
        )
    if knn_column_mean_fallback_imputation_mask_hash is not None:
        diagnostics["knn_column_mean_fallback_imputation_mask_hash"] = str(
            knn_column_mean_fallback_imputation_mask_hash
        )
    if knn_fully_column_mean_fallback_row_ids:
        diagnostics["knn_fully_column_mean_fallback_row_ids"] = list(
            knn_fully_column_mean_fallback_row_ids
        )
    if diagnostic_caveat_codes:
        diagnostics["diagnostic_caveat_codes"] = list(diagnostic_caveat_codes)
    if per_column_distribution_parameters is not None:
        per_column_distribution_payload: dict[str, JsonValue] = {
            str(column_name): dict(parameters)
            for column_name, parameters in per_column_distribution_parameters.items()
        }
        diagnostics["per_column_distribution_parameters"] = (
            per_column_distribution_payload
        )
    return diagnostics
