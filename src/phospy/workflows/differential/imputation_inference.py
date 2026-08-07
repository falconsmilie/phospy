"""Imputation-aware inferential status summaries for differential results."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.differential.models import (
    DIFFERENTIAL_RESULT_STATUS_TESTED,
)
from phospy.workflows.differential.models import (
    DifferentialFeatureEligibilityInputs,
    DifferentialImputationPolicyInputs,
)

DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_NO_TESTED_IMPUTED_VALUES = (
    "no_tested_imputed_values"
)
DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_RETAINED_IMPUTED_VALUES = (
    "retained_imputed_values_without_observed_only_fit"
)
DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_FULLY_OBSERVED = "tested_fully_observed"
DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_RETAINED_IMPUTED_VALUES = (
    "tested_retained_imputed_values"
)
DIFFERENTIAL_ROW_INFERENCE_STATUS_WITHHELD_NOT_TESTED = "withheld_not_tested"
DIFFERENTIAL_OBSERVED_ONLY_FIT = False
DIFFERENTIAL_RESIDUAL_DF_ADJUSTED_FOR_IMPUTATION = False


@dataclass(frozen=True, slots=True)
class DifferentialImputationInferenceSummary:
    """Tested-row imputation status for downstream inference interpretation."""

    total_feature_count: int
    tested_feature_count: int
    withheld_feature_count: int
    tested_imputed_feature_count: int
    tested_imputed_cell_count: int
    tested_imputed_feature_ids: tuple[str, ...]
    observed_only_fit: bool
    residual_df_adjusted_for_imputation: bool
    inferential_status: str
    adjusted_p_value_denominator_feature_count: int
    status_counts: dict[str, int]


def imputation_inference_columns(
    *,
    feature_metadata: pd.DataFrame,
    result_status: pd.Series,
) -> dict[str, NDArray[np.bool_] | NDArray[np.object_]]:
    """Build row-level inference columns from analysed imputation facts."""

    _require_imputation_inference_alignment(
        feature_metadata=feature_metadata,
        result_status=result_status,
    )
    imputed_cell_counts = _imputed_cell_counts(feature_metadata)
    status_values = result_status.astype(str).to_numpy(dtype=str)
    tested_mask = status_values == DIFFERENTIAL_RESULT_STATUS_TESTED
    contains_imputed_cells = imputed_cell_counts > 0
    inferential_status = np.full(
        int(status_values.size),
        DIFFERENTIAL_ROW_INFERENCE_STATUS_WITHHELD_NOT_TESTED,
        dtype=object,
    )
    inferential_status[tested_mask & ~contains_imputed_cells] = (
        DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_FULLY_OBSERVED
    )
    inferential_status[tested_mask & contains_imputed_cells] = (
        DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_RETAINED_IMPUTED_VALUES
    )
    return {
        "contains_imputed_cells": contains_imputed_cells,
        "observed_only_fit": np.full(
            int(status_values.size),
            DIFFERENTIAL_OBSERVED_ONLY_FIT,
            dtype=bool,
        ),
        "residual_df_adjusted_for_imputation": np.full(
            int(status_values.size),
            DIFFERENTIAL_RESIDUAL_DF_ADJUSTED_FOR_IMPUTATION,
            dtype=bool,
        ),
        "inferential_status": inferential_status,
    }


def summarize_differential_imputation_inference(
    *,
    imputation_policy_inputs: DifferentialImputationPolicyInputs,
    feature_eligibility_inputs: DifferentialFeatureEligibilityInputs | None = None,
) -> DifferentialImputationInferenceSummary:
    """Summarize retained imputed-row inference from final row eligibility."""

    feature_metadata = imputation_policy_inputs.feature_metadata
    result_status = imputation_policy_inputs.result_status
    if (
        feature_eligibility_inputs is not None
        and "imputed_cell_count" in feature_eligibility_inputs.feature_metadata.columns
    ):
        feature_metadata = feature_eligibility_inputs.feature_metadata
        result_status = feature_eligibility_inputs.result_status
    _require_imputation_inference_alignment(
        feature_metadata=feature_metadata,
        result_status=result_status,
    )
    imputed_cell_counts = _imputed_cell_counts(feature_metadata)
    status_values = result_status.astype(str).to_numpy(dtype=str)
    tested_mask = status_values == DIFFERENTIAL_RESULT_STATUS_TESTED
    tested_imputed_mask = tested_mask & (imputed_cell_counts > 0)
    tested_imputed_feature_ids = tuple(
        str(feature_id)
        for feature_id, retained in zip(
            feature_metadata.index,
            tested_imputed_mask,
            strict=True,
        )
        if bool(retained)
    )
    tested_imputed_cell_count = int(imputed_cell_counts[tested_imputed_mask].sum())
    tested_feature_count = int(np.count_nonzero(tested_mask))
    total_feature_count = int(status_values.size)
    status_counts = _status_counts(result_status)
    inferential_status = (
        DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_RETAINED_IMPUTED_VALUES
        if tested_imputed_cell_count > 0
        else DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_NO_TESTED_IMPUTED_VALUES
    )
    return DifferentialImputationInferenceSummary(
        total_feature_count=total_feature_count,
        tested_feature_count=tested_feature_count,
        withheld_feature_count=total_feature_count - tested_feature_count,
        tested_imputed_feature_count=len(tested_imputed_feature_ids),
        tested_imputed_cell_count=tested_imputed_cell_count,
        tested_imputed_feature_ids=tested_imputed_feature_ids,
        observed_only_fit=DIFFERENTIAL_OBSERVED_ONLY_FIT,
        residual_df_adjusted_for_imputation=(
            DIFFERENTIAL_RESIDUAL_DF_ADJUSTED_FOR_IMPUTATION
        ),
        inferential_status=inferential_status,
        adjusted_p_value_denominator_feature_count=tested_feature_count,
        status_counts=status_counts,
    )


def imputation_inference_summary_payload(
    summary: DifferentialImputationInferenceSummary,
    *,
    max_feature_examples: int = 5,
) -> dict[str, object]:
    """Return a compact JSON-compatible summary payload."""

    return {
        "total_feature_count": int(summary.total_feature_count),
        "tested_feature_count": int(summary.tested_feature_count),
        "withheld_feature_count": int(summary.withheld_feature_count),
        "tested_imputed_feature_count": int(summary.tested_imputed_feature_count),
        "tested_imputed_cell_count": int(summary.tested_imputed_cell_count),
        "tested_imputed_feature_examples": list(
            summary.tested_imputed_feature_ids[:max_feature_examples]
        ),
        "tested_imputed_feature_example_count": int(
            min(len(summary.tested_imputed_feature_ids), max_feature_examples)
        ),
        "observed_only_fit": bool(summary.observed_only_fit),
        "residual_df_adjusted_for_imputation": bool(
            summary.residual_df_adjusted_for_imputation
        ),
        "inferential_status": summary.inferential_status,
        "adjusted_p_value_denominator_feature_count": int(
            summary.adjusted_p_value_denominator_feature_count
        ),
        "status_counts": dict(summary.status_counts),
        "withholding_counts": {
            status: count
            for status, count in summary.status_counts.items()
            if status != DIFFERENTIAL_RESULT_STATUS_TESTED
        },
    }


def _imputed_cell_counts(feature_metadata: pd.DataFrame) -> np.ndarray:
    if "imputed_cell_count" not in feature_metadata.columns:
        raise WorkflowBoundaryError(
            seam="differential.imputation_inference.metadata_contract",
            next_action=(
                "ensure imputation-aware differential inference metadata is built "
                "from feature metadata that includes imputed_cell_count"
            ),
            message_prefix="differential workflow boundary validation failed",
        )
    return np.asarray(
        feature_metadata["imputed_cell_count"].to_numpy(dtype=np.int64),
        dtype=np.int64,
    )


def _status_counts(result_status: pd.Series) -> dict[str, int]:
    counts = result_status.astype(str).value_counts(sort=False)
    return {str(status): int(count) for status, count in counts.items()}


def _require_imputation_inference_alignment(
    *,
    feature_metadata: pd.DataFrame,
    result_status: pd.Series,
) -> None:
    if not feature_metadata.index.equals(result_status.index):
        raise WorkflowBoundaryError(
            seam="differential.imputation_inference.result_status_alignment",
            next_action=(
                "ensure imputation-aware differential inference metadata uses the "
                "same feature index as result_status"
            ),
            message_prefix="differential workflow boundary validation failed",
        )


__all__ = [
    "DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_NO_TESTED_IMPUTED_VALUES",
    "DIFFERENTIAL_IMPUTATION_INFERENCE_STATUS_RETAINED_IMPUTED_VALUES",
    "DIFFERENTIAL_OBSERVED_ONLY_FIT",
    "DIFFERENTIAL_RESIDUAL_DF_ADJUSTED_FOR_IMPUTATION",
    "DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_FULLY_OBSERVED",
    "DIFFERENTIAL_ROW_INFERENCE_STATUS_TESTED_RETAINED_IMPUTED_VALUES",
    "DIFFERENTIAL_ROW_INFERENCE_STATUS_WITHHELD_NOT_TESTED",
    "DifferentialImputationInferenceSummary",
    "imputation_inference_columns",
    "imputation_inference_summary_payload",
    "summarize_differential_imputation_inference",
]
