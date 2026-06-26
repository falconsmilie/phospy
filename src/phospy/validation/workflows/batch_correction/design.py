"""Design validation wrapper for batch-correction workflows."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from phospy.contracts.configs.preprocessing import (
    CorrectionMissingnessPolicy,
    InternalBatchCorrectionMethod,
    TemporaryImputationMethod,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.control_sites import ControlSiteMapping
from phospy.validation.datasets.batch_correction import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
    ResolvedBatchDesignMetadata,
)
from phospy.workflows.batch_correction.contracts import BatchCorrectionWorkflowRequest

_MIN_NUMERICAL_RANK_TOLERANCE = 0.0


class BatchCorrectionWorkflowDesignValidator:
    """Resolve sample metadata and validate batch/condition design adequacy."""

    def __init__(
        self,
        *,
        metadata_validator: BatchDesignMetadataValidator | None = None,
        adequacy_validator: BatchCorrectionAdequacyValidator | None = None,
    ) -> None:
        self._metadata_validator = metadata_validator or BatchDesignMetadataValidator()
        self._adequacy_validator = (
            adequacy_validator or BatchCorrectionAdequacyValidator()
        )

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
    ) -> ResolvedBatchDesignMetadata:
        config = request.config
        metadata = self._metadata_validator.run(
            phospho=request.phospho,
            sample_metadata=request.sample_metadata,
            batch_column=config.batch_column,
            condition_columns=config.condition_columns,
            replicate_column=config.replicate_column,
            require_replicate_column=(
                config.method is InternalBatchCorrectionMethod.RUV_III_STYLE
            ),
            context="batch-correction workflow",
        )
        self._adequacy_validator.run(
            batch_by_sample=metadata.batch_by_sample,
            condition_by_sample=metadata.condition_by_sample,
            sample_order=metadata.sample_order,
            preserve_condition_effects=True,
            context="SPS/RUV-style batch correction design validation",
        )
        return metadata


class BatchCorrectionWorkflowFactorFeasibilityValidator:
    """Validate requested unwanted-factor count before numerical execution."""

    def run(
        self,
        *,
        request: BatchCorrectionWorkflowRequest,
        dataset_metadata: ResolvedBatchDesignMetadata,
        control_site_mapping: ControlSiteMapping,
        missingness_policy: CorrectionMissingnessPolicy,
    ) -> None:
        requested = request.config.n_unwanted_factors
        if requested is None:
            requested = 1
        if isinstance(requested, bool) or not isinstance(requested, int):
            raise PhosPyInputError(
                "batch-correction factor feasibility validation requires "
                "n_unwanted_factors to be an int when provided"
            )
        if requested < 1:
            raise PhosPyInputError(
                "batch-correction factor feasibility validation requires "
                "n_unwanted_factors >= 1"
            )

        controls = tuple(
            row for row in control_site_mapping.row_eligibility if row.is_control
        )
        control_count = len(controls)
        if control_count <= requested:
            raise PhosPyInputError(
                "batch-correction factor feasibility validation failed: eligible "
                f"controls are too few for n_unwanted_factors={requested}; "
                f"observed {control_count}, required at least {requested + 1}"
            )

        sample_order = tuple(dataset_metadata.sample_order)
        protected = _condition_design(dataset_metadata)
        protected_rank = _matrix_rank(protected)
        sample_residual_df = len(sample_order) - protected_rank
        if sample_residual_df < requested:
            raise PhosPyInputError(
                "batch-correction factor feasibility validation failed: requested "
                f"n_unwanted_factors={requested} exceeds sample/design rank "
                f"capacity after protected terms (samples={len(sample_order)}, "
                f"protected_design_rank={protected_rank}, "
                f"residual_degrees_of_freedom={sample_residual_df})"
            )

        working = _prepared_rank_matrix(
            phospho=request.phospho,
            missingness_policy=missingness_policy,
        )
        control_positions = _control_positions(
            phospho=request.phospho,
            controls=controls,
        )
        response = working.to_numpy(dtype="float64", copy=True).T
        protected_coefficients = np.linalg.pinv(protected) @ response
        control_residual = (response - protected @ protected_coefficients)[
            :, list(control_positions)
        ]
        control_residual_rank = _matrix_rank(control_residual)
        if control_residual_rank < requested:
            raise PhosPyInputError(
                "batch-correction factor feasibility validation failed: requested "
                f"n_unwanted_factors={requested} exceeds eligible control "
                f"residual rank {control_residual_rank}; eligible controls do "
                "not contain enough independent variation after protected design "
                "terms"
            )


def _condition_design(metadata: ResolvedBatchDesignMetadata) -> np.ndarray:
    return _treatment_coded_design(metadata.condition_labels, include_intercept=True)


def _treatment_coded_design(
    labels: Sequence[str],
    *,
    include_intercept: bool,
) -> np.ndarray:
    levels = _levels_in_order(labels)
    row_width = (1 if include_intercept else 0) + max(len(levels) - 1, 0)
    if row_width == 0:
        return np.empty((len(labels), 0), dtype="float64")
    rows: list[list[float]] = []
    for label in labels:
        row: list[float] = []
        if include_intercept:
            row.append(1.0)
        row.extend(1.0 if label == level else 0.0 for level in levels[1:])
        rows.append(row)
    return np.asarray(rows, dtype="float64")


def _levels_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    levels: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        levels.append(label)
    return tuple(levels)


def _prepared_rank_matrix(
    *,
    phospho: pd.DataFrame,
    missingness_policy: CorrectionMissingnessPolicy,
) -> pd.DataFrame:
    matrix = phospho.astype("float64").copy(deep=True)
    if not bool(matrix.isna().to_numpy().any()):
        return matrix
    policy = missingness_policy.temporary_imputation
    method = TemporaryImputationMethod.parse(
        policy.method,
        field_name="temporary imputation policy.method",
    )
    if (
        method is not TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY
        or not policy.allowed
    ):
        raise PhosPyInputError(
            "batch-correction factor feasibility validation requires executable "
            "row_median_temporary imputation before rank checks when missing "
            "values are present"
        )
    min_observed_values = dict(policy.method_parameters).get("min_observed_values", 1)
    if isinstance(min_observed_values, bool) or not isinstance(
        min_observed_values, int
    ):
        raise PhosPyInputError(
            "batch-correction factor feasibility validation requires integer "
            "row_median_temporary min_observed_values"
        )
    for row_id in matrix.index.tolist():
        row = matrix.loc[row_id, :]
        observed = row.dropna()
        if int(observed.shape[0]) < int(min_observed_values):
            raise PhosPyInputError(
                "batch-correction factor feasibility validation cannot temporarily "
                f"impute row {str(row_id)!r}; observed values are below "
                "min_observed_values"
            )
        missing = row.isna()
        if bool(missing.any()):
            matrix.loc[row_id, missing] = float(observed.median())
    if not np.isfinite(matrix.to_numpy(dtype="float64", copy=True)).all():
        raise PhosPyInputError(
            "batch-correction factor feasibility validation requires finite values "
            "after temporary imputation"
        )
    return matrix


def _control_positions(
    *,
    phospho: pd.DataFrame,
    controls: Sequence[object],
) -> tuple[int, ...]:
    row_count = int(phospho.shape[0])
    positions: list[int] = []
    for row in controls:
        position = int(getattr(row, "row_position", -1))
        if position < 0 or position >= row_count:
            raise PhosPyInputError(
                "batch-correction factor feasibility validation found a control "
                "row_position outside the phospho matrix"
            )
        expected = str(getattr(row, "site_key", ""))
        observed = str(phospho.index[position])
        if observed != expected:
            raise PhosPyInputError(
                "batch-correction factor feasibility validation found a control "
                "row_position that does not match the phospho index "
                f"(expected {expected!r}, observed {observed!r})"
            )
        positions.append(position)
    if len(set(positions)) != len(positions):
        raise PhosPyInputError(
            "batch-correction factor feasibility validation requires unique "
            "eligible control rows"
        )
    return tuple(positions)


def _matrix_rank(matrix: np.ndarray) -> int:
    if matrix.size == 0:
        return 0
    if _MIN_NUMERICAL_RANK_TOLERANCE > 0.0:
        return int(np.linalg.matrix_rank(matrix, tol=_MIN_NUMERICAL_RANK_TOLERANCE))
    return int(np.linalg.matrix_rank(matrix))


__all__ = [
    "BatchCorrectionWorkflowDesignValidator",
    "BatchCorrectionWorkflowFactorFeasibilityValidator",
]
