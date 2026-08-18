"""Feature-wise REML duplicate-correlation estimation.

The implementation in this module belongs to the differential science domain.
It estimates feature-level compound-symmetry correlations from the fixed-effects
design and block identities, then aggregates finite estimates on Fisher's
``atanh`` scale using the internal duplicate-correlation contract.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import cast

import numpy as np
import numpy.typing as npt
from scipy.linalg import solve_triangular

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.compound_symmetry_gls import (
    CompoundSymmetryGLSFit,
    fit_compound_symmetry_gls,
    fit_duplicate_correlation_gls,
)
from phospy.science.differential.models.duplicate_correlation import (
    DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
    DUPLICATE_CORRELATION_MINIMUM_RESIDUAL_DEGREES_OF_FREEDOM,
    DUPLICATE_CORRELATION_TRIM_FRACTION,
    DuplicateCorrelationBlockStructureSummary,
    DuplicateCorrelationBoundary,
    DuplicateCorrelationBoundarySummary,
    DuplicateCorrelationConsensusResult,
    DuplicateCorrelationConvergenceSummary,
    DuplicateCorrelationFailureReason,
    DuplicateCorrelationFeatureEstimate,
    DuplicateCorrelationFeatureStatus,
    DuplicateCorrelationReasonCount,
)

DUPLICATE_CORRELATION_POSITIVE_DEFINITE_TOLERANCE = 1.0e-10
DUPLICATE_CORRELATION_FISHER_BOUNDARY_TOLERANCE = (
    DUPLICATE_CORRELATION_POSITIVE_DEFINITE_TOLERANCE
)
DUPLICATE_CORRELATION_UPPER_CORRELATION_CAP = 0.99
DUPLICATE_CORRELATION_OPTIMIZER_ABSOLUTE_TOLERANCE = 1.0e-10
DUPLICATE_CORRELATION_OPTIMIZER_MAX_ITERATIONS = 500
DUPLICATE_CORRELATION_BOUNDARY_DETECTION_TOLERANCE = 1.0e-7

_ZERO_RESIDUAL_VARIATION_RELATIVE_TOLERANCE = np.finfo(np.float64).eps ** 2


@dataclass(frozen=True, slots=True)
class _FeatureBounds:
    lower: float
    upper: float


@dataclass(frozen=True, slots=True)
class _OptimizerResult:
    converged: bool
    estimate: float
    objective_value: float
    iteration_count: int


@dataclass(frozen=True, slots=True)
class _TrimmedConsensus:
    correlation: float
    trim_count_each_tail: int
    retained_count: int


def estimate_duplicate_correlation_reml_consensus(
    matrix: npt.ArrayLike,
    design: npt.ArrayLike,
    block_ids: Sequence[object],
    *,
    feature_ids: Sequence[object] | None = None,
    design_column_names: Sequence[object] | None = None,
    block_id_field_name: str = "block_id",
    optimizer_absolute_tolerance: float = (
        DUPLICATE_CORRELATION_OPTIMIZER_ABSOLUTE_TOLERANCE
    ),
    optimizer_max_iterations: int = DUPLICATE_CORRELATION_OPTIMIZER_MAX_ITERATIONS,
) -> DuplicateCorrelationConsensusResult:
    """Estimate feature-wise REML correlations and one Fisher-trimmed consensus.

    ``matrix`` is feature-by-sample. ``design`` is sample-by-fixed-effect.
    ``block_ids`` is aligned to the design rows and matrix columns. Missing
    feature observations are represented by non-finite matrix values and are
    removed feature-by-feature before fitting.
    """

    matrix_values = _as_float_matrix(
        matrix,
        field_name="duplicate_correlation.matrix",
    )
    design_values = _as_float_matrix(
        design,
        field_name="duplicate_correlation.design",
    )
    if matrix_values.shape[1] != design_values.shape[0]:
        raise PhosPyInputError(
            "duplicate_correlation.matrix column count must match "
            "duplicate_correlation.design row count; "
            f"matrix_columns={int(matrix_values.shape[1])}, "
            f"design_rows={int(design_values.shape[0])}"
        )
    if matrix_values.shape[0] < 1:
        raise PhosPyInputError(
            "duplicate_correlation.matrix must contain at least one feature"
        )
    if design_values.shape[1] < 1:
        raise PhosPyInputError(
            "duplicate_correlation.design must contain at least one fixed effect"
        )
    if not np.isfinite(design_values).all():
        raise PhosPyInputError(
            "duplicate_correlation.design must contain only finite numeric values"
        )

    feature_id_values = _coerce_feature_ids(
        feature_ids,
        feature_count=int(matrix_values.shape[0]),
    )
    block_id_values = _coerce_block_ids(
        block_ids,
        expected_count=int(design_values.shape[0]),
    )
    _reject_block_design_columns(
        design_column_names,
        coefficient_count=int(design_values.shape[1]),
    )
    block_structure = _summarize_blocks(
        block_id_values,
        block_id_field_name=block_id_field_name,
    )
    design_rank = _scaled_design_rank(design_values)
    if design_rank < int(design_values.shape[1]):
        raise PhosPyInputError(
            "duplicate_correlation.design is rank deficient under scaled SVD "
            "rank assessment; "
            f"rank={design_rank}, columns={int(design_values.shape[1])}"
        )
    residual_degrees_of_freedom = float(int(design_values.shape[0]) - design_rank)
    if (
        residual_degrees_of_freedom
        < DUPLICATE_CORRELATION_MINIMUM_RESIDUAL_DEGREES_OF_FREEDOM
    ):
        raise PhosPyInputError(
            "duplicate_correlation.design requires at least two residual degrees "
            "of freedom for REML correlation estimation; "
            f"samples={int(design_values.shape[0])}, "
            f"design_rank={design_rank}, "
            f"residual_degrees_of_freedom={residual_degrees_of_freedom}"
        )

    optimizer_absolute_tolerance = _require_positive_float(
        optimizer_absolute_tolerance,
        field_name="duplicate_correlation.optimizer_absolute_tolerance",
    )
    optimizer_max_iterations = _require_positive_int(
        optimizer_max_iterations,
        field_name="duplicate_correlation.optimizer_max_iterations",
    )

    estimates = tuple(
        _estimate_one_feature(
            feature_id=feature_id,
            values=matrix_values[feature_position, :],
            design=design_values,
            block_ids=block_id_values,
            optimizer_absolute_tolerance=optimizer_absolute_tolerance,
            optimizer_max_iterations=optimizer_max_iterations,
        )
        for feature_position, feature_id in enumerate(feature_id_values)
    )
    return _build_consensus_result(
        estimates=estimates,
        attempted_feature_count=int(matrix_values.shape[0]),
        block_structure=block_structure,
        design_rank=design_rank,
        sample_count=int(design_values.shape[0]),
    )


def _estimate_one_feature(
    *,
    feature_id: str,
    values: npt.NDArray[np.float64],
    design: npt.NDArray[np.float64],
    block_ids: tuple[str, ...],
    optimizer_absolute_tolerance: float,
    optimizer_max_iterations: int,
) -> DuplicateCorrelationFeatureEstimate:
    finite_mask = np.isfinite(values)
    observed_count = int(np.count_nonzero(finite_mask))
    coefficient_count = int(design.shape[1])
    if observed_count < coefficient_count:
        return _failed_feature_estimate(
            feature_id=feature_id,
            status=DuplicateCorrelationFeatureStatus.INSUFFICIENT_FINITE_OBSERVATIONS,
            failure_reason=(
                DuplicateCorrelationFailureReason.INSUFFICIENT_FINITE_OBSERVATIONS
            ),
            observed_value_count=observed_count,
            design_rank=None,
            residual_degrees_of_freedom=None,
        )

    observed_values = values[finite_mask]
    observed_design = design[finite_mask, :]
    observed_blocks = tuple(
        block_id
        for block_id, observed in zip(block_ids, finite_mask, strict=True)
        if observed
    )
    observed_design_rank = _scaled_design_rank(observed_design)
    if observed_design_rank < coefficient_count:
        return _failed_feature_estimate(
            feature_id=feature_id,
            status=DuplicateCorrelationFeatureStatus.LOST_FIXED_EFFECT_ESTIMABILITY,
            failure_reason=(
                DuplicateCorrelationFailureReason.LOST_FIXED_EFFECT_ESTIMABILITY
            ),
            observed_value_count=observed_count,
            design_rank=observed_design_rank,
            residual_degrees_of_freedom=None,
        )

    residual_dof = float(observed_count - observed_design_rank)
    if residual_dof < DUPLICATE_CORRELATION_MINIMUM_RESIDUAL_DEGREES_OF_FREEDOM:
        return _failed_feature_estimate(
            feature_id=feature_id,
            status=(
                DuplicateCorrelationFeatureStatus.INSUFFICIENT_RESIDUAL_DEGREES_OF_FREEDOM
            ),
            failure_reason=(
                DuplicateCorrelationFailureReason.INSUFFICIENT_RESIDUAL_DEGREES_OF_FREEDOM
            ),
            observed_value_count=observed_count,
            design_rank=observed_design_rank,
            residual_degrees_of_freedom=residual_dof,
        )

    max_observed_block_size = _max_block_size(observed_blocks)
    if max_observed_block_size < 2:
        return _failed_feature_estimate(
            feature_id=feature_id,
            status=DuplicateCorrelationFeatureStatus.NO_REPEATED_OBSERVATIONS,
            failure_reason=(
                DuplicateCorrelationFailureReason.NO_REPEATED_OBSERVATIONS_AFTER_SUBSETTING
            ),
            observed_value_count=observed_count,
            design_rank=observed_design_rank,
            residual_degrees_of_freedom=residual_dof,
        )

    residual_sum_of_squares = _ols_residual_sum_of_squares(
        observed_values,
        observed_design,
    )
    residual_tolerance = _residual_variation_tolerance(
        observed_values,
        sample_count=observed_count,
        coefficient_count=coefficient_count,
    )
    if (
        not math.isfinite(residual_sum_of_squares)
        or residual_sum_of_squares <= residual_tolerance
    ):
        return _failed_feature_estimate(
            feature_id=feature_id,
            status=(
                DuplicateCorrelationFeatureStatus.ZERO_OR_UNUSABLE_RESIDUAL_VARIATION
            ),
            failure_reason=(
                DuplicateCorrelationFailureReason.ZERO_OR_UNUSABLE_RESIDUAL_VARIATION
            ),
            observed_value_count=observed_count,
            design_rank=observed_design_rank,
            residual_degrees_of_freedom=residual_dof,
        )

    bounds = _feature_bounds(max_observed_block_size)

    def objective(correlation: float) -> float:
        return _restricted_reml_objective(
            correlation=correlation,
            observed_values=observed_values,
            observed_design=observed_design,
            observed_blocks=observed_blocks,
            residual_degrees_of_freedom=residual_dof,
        )

    optimizer_result = _bounded_golden_section_minimize(
        objective,
        lower_bound=bounds.lower,
        upper_bound=bounds.upper,
        absolute_tolerance=optimizer_absolute_tolerance,
        max_iterations=optimizer_max_iterations,
    )
    if (
        not math.isfinite(optimizer_result.estimate)
        or not math.isfinite(optimizer_result.objective_value)
        or not bounds.lower <= optimizer_result.estimate <= bounds.upper
    ):
        return _failed_feature_estimate(
            feature_id=feature_id,
            status=(DuplicateCorrelationFeatureStatus.NON_FINITE_OBJECTIVE_OR_ESTIMATE),
            failure_reason=(
                DuplicateCorrelationFailureReason.NON_FINITE_OBJECTIVE_OR_ESTIMATE
            ),
            observed_value_count=observed_count,
            design_rank=observed_design_rank,
            residual_degrees_of_freedom=residual_dof,
            lower_correlation_bound=bounds.lower,
            upper_correlation_bound=bounds.upper,
        )
    if not optimizer_result.converged:
        return _failed_feature_estimate(
            feature_id=feature_id,
            status=DuplicateCorrelationFeatureStatus.OPTIMISATION_FAILED,
            failure_reason=DuplicateCorrelationFailureReason.OPTIMISATION_FAILED,
            observed_value_count=observed_count,
            design_rank=observed_design_rank,
            residual_degrees_of_freedom=residual_dof,
            lower_correlation_bound=bounds.lower,
            upper_correlation_bound=bounds.upper,
        )

    boundary = _boundary_kind(
        optimizer_result.estimate,
        lower_bound=bounds.lower,
        upper_bound=bounds.upper,
    )
    return DuplicateCorrelationFeatureEstimate(
        feature_id=feature_id,
        status=(
            DuplicateCorrelationFeatureStatus.ESTIMATED
            if boundary is None
            else DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED
        ),
        correlation=float(optimizer_result.estimate),
        observed_value_count=observed_count,
        design_rank=observed_design_rank,
        residual_degrees_of_freedom=residual_dof,
        objective_value=float(optimizer_result.objective_value),
        lower_correlation_bound=bounds.lower,
        upper_correlation_bound=bounds.upper,
        boundary=boundary,
    )


def _build_consensus_result(
    *,
    estimates: tuple[DuplicateCorrelationFeatureEstimate, ...],
    attempted_feature_count: int,
    block_structure: DuplicateCorrelationBlockStructureSummary,
    design_rank: int,
    sample_count: int,
) -> DuplicateCorrelationConsensusResult:
    estimated_correlations = tuple(
        cast(float, estimate.correlation)
        for estimate in estimates
        if estimate.status
        in {
            DuplicateCorrelationFeatureStatus.ESTIMATED,
            DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED,
        }
    )
    estimated_feature_count = len(estimated_correlations)
    failed_feature_count = len(estimates) - estimated_feature_count
    non_finite_feature_count = sum(
        estimate.status
        is DuplicateCorrelationFeatureStatus.NON_FINITE_OBJECTIVE_OR_ESTIMATE
        for estimate in estimates
    )
    reason_counts = _failure_reason_counts(estimates)
    convergence_summary = _convergence_summary(estimates)
    boundary_summary = _boundary_summary(
        block_structure=block_structure,
        estimates=estimates,
    )

    if estimated_feature_count == 0:
        return DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=False,
            consensus_correlation=None,
            attempted_feature_count=attempted_feature_count,
            eligible_feature_count=len(estimates),
            estimated_feature_count=0,
            failed_feature_count=failed_feature_count,
            non_finite_feature_count=non_finite_feature_count,
            failure_reason=_empty_consensus_failure_reason(
                estimates,
                block_structure=block_structure,
            ),
            feature_estimates=estimates,
            trimmed_feature_count_each_tail=0,
            retained_feature_count_after_trimming=0,
            failure_reason_counts=reason_counts,
            convergence_summary=convergence_summary,
            boundary_summary=boundary_summary,
            block_structure=block_structure,
            design_rank=design_rank,
            sample_count=sample_count,
        )

    consensus = _fisher_trimmed_mean_consensus(estimated_correlations)
    if not _correlation_within_workflow_bounds(
        consensus.correlation,
        block_structure=block_structure,
    ):
        return DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=False,
            consensus_correlation=None,
            attempted_feature_count=attempted_feature_count,
            eligible_feature_count=len(estimates),
            estimated_feature_count=estimated_feature_count,
            failed_feature_count=failed_feature_count,
            non_finite_feature_count=non_finite_feature_count,
            failure_reason=(
                DuplicateCorrelationFailureReason.INVALID_OR_NON_POSITIVE_DEFINITE_COVARIANCE
            ),
            feature_estimates=estimates,
            trimmed_feature_count_each_tail=consensus.trim_count_each_tail,
            retained_feature_count_after_trimming=consensus.retained_count,
            failure_reason_counts=reason_counts,
            convergence_summary=convergence_summary,
            boundary_summary=boundary_summary,
            block_structure=block_structure,
            design_rank=design_rank,
            sample_count=sample_count,
        )

    return DuplicateCorrelationConsensusResult(
        method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
        trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
        success=True,
        consensus_correlation=consensus.correlation,
        attempted_feature_count=attempted_feature_count,
        eligible_feature_count=len(estimates),
        estimated_feature_count=estimated_feature_count,
        failed_feature_count=failed_feature_count,
        non_finite_feature_count=non_finite_feature_count,
        failure_reason=None,
        feature_estimates=estimates,
        trimmed_feature_count_each_tail=consensus.trim_count_each_tail,
        retained_feature_count_after_trimming=consensus.retained_count,
        failure_reason_counts=reason_counts,
        convergence_summary=convergence_summary,
        boundary_summary=boundary_summary,
        block_structure=block_structure,
        design_rank=design_rank,
        sample_count=sample_count,
    )


def _restricted_reml_objective(
    *,
    correlation: float,
    observed_values: npt.NDArray[np.float64],
    observed_design: npt.NDArray[np.float64],
    observed_blocks: tuple[str, ...],
    residual_degrees_of_freedom: float,
) -> float:
    covariance = _compound_symmetry_correlation_matrix(
        correlation,
        block_ids=observed_blocks,
    )
    try:
        cholesky_covariance = np.asarray(
            np.linalg.cholesky(covariance),
            dtype=np.float64,
        )
        logdet_covariance = _cholesky_logdet(cholesky_covariance)
        covariance_solved_design = _cholesky_solve(
            cholesky_covariance,
            observed_design,
        )
        gls_gram = np.asarray(
            observed_design.T @ covariance_solved_design,
            dtype=np.float64,
        )
        cholesky_gls_gram = np.asarray(
            np.linalg.cholesky(gls_gram),
            dtype=np.float64,
        )
        logdet_gls_gram = _cholesky_logdet(cholesky_gls_gram)
        covariance_solved_values = _cholesky_solve(
            cholesky_covariance,
            observed_values,
        )
        gls_rhs = np.asarray(
            observed_design.T @ covariance_solved_values,
            dtype=np.float64,
        )
        coefficients = _cholesky_solve(cholesky_gls_gram, gls_rhs)
        residuals = np.asarray(
            observed_values - (observed_design @ coefficients),
            dtype=np.float64,
        )
        whitened_residuals = solve_triangular(
            cholesky_covariance,
            residuals,
            lower=True,
            check_finite=False,
        )
        residual_quadratic_form = float(whitened_residuals @ whitened_residuals)
    except (ValueError, np.linalg.LinAlgError):
        return math.inf

    if (
        not math.isfinite(logdet_covariance)
        or not math.isfinite(logdet_gls_gram)
        or not math.isfinite(residual_quadratic_form)
        or residual_quadratic_form <= 0.0
        or residual_degrees_of_freedom <= 0.0
    ):
        return math.inf
    objective = (
        logdet_covariance
        + logdet_gls_gram
        + residual_degrees_of_freedom
        * math.log(residual_quadratic_form / residual_degrees_of_freedom)
    )
    if not math.isfinite(objective):
        return math.inf
    return float(objective)


def _bounded_golden_section_minimize(
    objective: Callable[[float], float],
    *,
    lower_bound: float,
    upper_bound: float,
    absolute_tolerance: float,
    max_iterations: int,
) -> _OptimizerResult:
    if not lower_bound < upper_bound:
        return _OptimizerResult(
            converged=False,
            estimate=math.nan,
            objective_value=math.inf,
            iteration_count=0,
        )
    inverse_phi = (math.sqrt(5.0) - 1.0) / 2.0
    lower = float(lower_bound)
    upper = float(upper_bound)
    left = upper - inverse_phi * (upper - lower)
    right = lower + inverse_phi * (upper - lower)
    left_value = objective(left)
    right_value = objective(right)
    best_estimate = lower_bound
    best_objective = objective(lower_bound)
    best_estimate, best_objective = _prefer_lower_objective(
        current_estimate=best_estimate,
        current_objective=best_objective,
        candidate_estimate=upper_bound,
        candidate_objective=objective(upper_bound),
    )
    best_estimate, best_objective = _prefer_lower_objective(
        current_estimate=best_estimate,
        current_objective=best_objective,
        candidate_estimate=left,
        candidate_objective=left_value,
    )
    best_estimate, best_objective = _prefer_lower_objective(
        current_estimate=best_estimate,
        current_objective=best_objective,
        candidate_estimate=right,
        candidate_objective=right_value,
    )

    iteration_count = 0
    for iteration in range(1, max_iterations + 1):
        iteration_count = iteration
        if abs(upper - lower) <= absolute_tolerance:
            break
        if left_value <= right_value:
            upper = right
            right = left
            right_value = left_value
            left = upper - inverse_phi * (upper - lower)
            left_value = objective(left)
            candidate_estimate = left
            candidate_objective = left_value
        else:
            lower = left
            left = right
            left_value = right_value
            right = lower + inverse_phi * (upper - lower)
            right_value = objective(right)
            candidate_estimate = right
            candidate_objective = right_value
        best_estimate, best_objective = _prefer_lower_objective(
            current_estimate=best_estimate,
            current_objective=best_objective,
            candidate_estimate=candidate_estimate,
            candidate_objective=candidate_objective,
        )

    converged = abs(upper - lower) <= absolute_tolerance and math.isfinite(
        best_objective
    )
    return _OptimizerResult(
        converged=converged,
        estimate=float(best_estimate),
        objective_value=float(best_objective),
        iteration_count=iteration_count,
    )


def _prefer_lower_objective(
    *,
    current_estimate: float,
    current_objective: float,
    candidate_estimate: float,
    candidate_objective: float,
) -> tuple[float, float]:
    if not math.isfinite(candidate_objective):
        return current_estimate, current_objective
    if not math.isfinite(current_objective) or candidate_objective < current_objective:
        return candidate_estimate, candidate_objective
    return current_estimate, current_objective


def _compound_symmetry_correlation_matrix(
    correlation: float,
    *,
    block_ids: tuple[str, ...],
) -> npt.NDArray[np.float64]:
    sample_count = len(block_ids)
    covariance = np.eye(sample_count, dtype=np.float64)
    for row in range(sample_count):
        for column in range(row + 1, sample_count):
            if block_ids[row] == block_ids[column]:
                covariance[row, column] = correlation
                covariance[column, row] = correlation
    return covariance


def _cholesky_solve(
    cholesky_factor: npt.NDArray[np.float64],
    rhs: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    forward = cast(
        npt.NDArray[np.float64],
        solve_triangular(
            cholesky_factor,
            rhs,
            lower=True,
            check_finite=False,
        ),
    )
    return cast(
        npt.NDArray[np.float64],
        solve_triangular(
            cholesky_factor.T,
            forward,
            lower=False,
            check_finite=False,
        ),
    )


def _cholesky_logdet(cholesky_factor: npt.NDArray[np.float64]) -> float:
    diagonal = np.diagonal(cholesky_factor)
    if not np.isfinite(diagonal).all() or np.any(diagonal <= 0.0):
        return math.inf
    return float(2.0 * np.sum(np.log(diagonal)))


def _fisher_trimmed_mean_consensus(
    correlations: Sequence[float],
) -> _TrimmedConsensus:
    transformed = sorted(
        math.atanh(_require_correlation(value)) for value in correlations
    )
    estimate_count = len(transformed)
    if estimate_count == 0:
        raise ValueError("at least one finite feature correlation is required")
    trim_count = int(
        math.floor(float(estimate_count) * DUPLICATE_CORRELATION_TRIM_FRACTION)
    )
    retained = transformed[trim_count : estimate_count - trim_count]
    if not retained:
        raise ValueError("duplicate-correlation trimming retained no estimates")
    mean_transformed = math.fsum(retained) / float(len(retained))
    correlation = math.tanh(mean_transformed)
    if not math.isfinite(correlation) or not -1.0 < correlation < 1.0:
        raise ValueError("duplicate-correlation consensus is not finite")
    return _TrimmedConsensus(
        correlation=float(correlation),
        trim_count_each_tail=trim_count,
        retained_count=len(retained),
    )


def _failure_reason_counts(
    estimates: tuple[DuplicateCorrelationFeatureEstimate, ...],
) -> tuple[DuplicateCorrelationReasonCount, ...]:
    counts: Counter[DuplicateCorrelationFailureReason] = Counter(
        estimate.failure_reason
        for estimate in estimates
        if estimate.failure_reason is not None
    )
    return tuple(
        DuplicateCorrelationReasonCount(reason=reason, count=counts[reason])
        for reason in sorted(counts, key=lambda value: value.value)
    )


def _convergence_summary(
    estimates: tuple[DuplicateCorrelationFeatureEstimate, ...],
) -> DuplicateCorrelationConvergenceSummary:
    converged = sum(
        estimate.status
        in {
            DuplicateCorrelationFeatureStatus.ESTIMATED,
            DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED,
        }
        for estimate in estimates
    )
    boundary = sum(
        estimate.status is DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED
        for estimate in estimates
    )
    failed_optimisation = sum(
        estimate.status is DuplicateCorrelationFeatureStatus.OPTIMISATION_FAILED
        for estimate in estimates
    )
    non_finite = sum(
        estimate.status
        is DuplicateCorrelationFeatureStatus.NON_FINITE_OBJECTIVE_OR_ESTIMATE
        for estimate in estimates
    )
    return DuplicateCorrelationConvergenceSummary(
        converged_feature_count=converged,
        boundary_feature_count=boundary,
        failed_optimisation_feature_count=failed_optimisation,
        non_finite_objective_or_estimate_feature_count=non_finite,
    )


def _boundary_summary(
    *,
    block_structure: DuplicateCorrelationBlockStructureSummary,
    estimates: tuple[DuplicateCorrelationFeatureEstimate, ...],
) -> DuplicateCorrelationBoundarySummary:
    bounds = _workflow_bounds(block_structure)
    lower_count = sum(
        estimate.boundary is DuplicateCorrelationBoundary.LOWER
        for estimate in estimates
    )
    upper_count = sum(
        estimate.boundary is DuplicateCorrelationBoundary.UPPER
        for estimate in estimates
    )
    return DuplicateCorrelationBoundarySummary(
        lower_correlation_bound=bounds.lower,
        upper_correlation_bound=bounds.upper,
        lower_boundary_feature_count=lower_count,
        upper_boundary_feature_count=upper_count,
        positive_definite_tolerance=(DUPLICATE_CORRELATION_POSITIVE_DEFINITE_TOLERANCE),
        fisher_boundary_tolerance=DUPLICATE_CORRELATION_FISHER_BOUNDARY_TOLERANCE,
    )


def _empty_consensus_failure_reason(
    estimates: tuple[DuplicateCorrelationFeatureEstimate, ...],
    *,
    block_structure: DuplicateCorrelationBlockStructureSummary,
) -> DuplicateCorrelationFailureReason:
    if block_structure.repeated_block_count == 0:
        return DuplicateCorrelationFailureReason.NO_REPEATED_BLOCKS
    statuses = {estimate.status for estimate in estimates}
    if statuses == {DuplicateCorrelationFeatureStatus.OPTIMISATION_FAILED}:
        return DuplicateCorrelationFailureReason.NUMERICAL_NON_CONVERGENCE
    if statuses == {DuplicateCorrelationFeatureStatus.NON_FINITE_OBJECTIVE_OR_ESTIMATE}:
        return (
            DuplicateCorrelationFailureReason.ALL_ELIGIBLE_FEATURE_ESTIMATES_NON_FINITE
        )
    return DuplicateCorrelationFailureReason.NO_FEATURE_WITH_ESTIMABLE_CORRELATION


def _correlation_within_workflow_bounds(
    correlation: float,
    *,
    block_structure: DuplicateCorrelationBlockStructureSummary,
) -> bool:
    bounds = _workflow_bounds(block_structure)
    return bounds.lower <= correlation <= bounds.upper


def _workflow_bounds(
    block_structure: DuplicateCorrelationBlockStructureSummary,
) -> _FeatureBounds:
    if block_structure.repeated_block_count == 0:
        return _FeatureBounds(
            lower=-1.0 + DUPLICATE_CORRELATION_FISHER_BOUNDARY_TOLERANCE,
            upper=1.0 - DUPLICATE_CORRELATION_FISHER_BOUNDARY_TOLERANCE,
        )
    if block_structure.maximum_block_size is None:
        raise PhosPyInputError(
            "duplicate_correlation.block_structure.maximum_block_size is required "
            "for consensus covariance validation"
        )
    maximum_block_size = block_structure.maximum_block_size
    return _feature_bounds(maximum_block_size)


def _feature_bounds(maximum_observed_block_size: int) -> _FeatureBounds:
    if maximum_observed_block_size < 2:
        return _FeatureBounds(
            lower=-1.0 + DUPLICATE_CORRELATION_FISHER_BOUNDARY_TOLERANCE,
            upper=1.0 - DUPLICATE_CORRELATION_FISHER_BOUNDARY_TOLERANCE,
        )
    mathematical_lower = -1.0 / float(maximum_observed_block_size - 1)
    lower = max(
        mathematical_lower + DUPLICATE_CORRELATION_POSITIVE_DEFINITE_TOLERANCE,
        -1.0 + DUPLICATE_CORRELATION_POSITIVE_DEFINITE_TOLERANCE,
    )
    upper = min(
        1.0 - DUPLICATE_CORRELATION_POSITIVE_DEFINITE_TOLERANCE,
        DUPLICATE_CORRELATION_UPPER_CORRELATION_CAP,
    )
    if not lower < upper:
        lower = math.nextafter(mathematical_lower, upper)
    return _FeatureBounds(lower=float(lower), upper=float(upper))


def _boundary_kind(
    estimate: float,
    *,
    lower_bound: float,
    upper_bound: float,
) -> DuplicateCorrelationBoundary | None:
    if estimate - lower_bound <= DUPLICATE_CORRELATION_BOUNDARY_DETECTION_TOLERANCE:
        return DuplicateCorrelationBoundary.LOWER
    if upper_bound - estimate <= DUPLICATE_CORRELATION_BOUNDARY_DETECTION_TOLERANCE:
        return DuplicateCorrelationBoundary.UPPER
    return None


def _failed_feature_estimate(
    *,
    feature_id: str,
    status: DuplicateCorrelationFeatureStatus,
    failure_reason: DuplicateCorrelationFailureReason,
    observed_value_count: int,
    design_rank: int | None,
    residual_degrees_of_freedom: float | None,
    lower_correlation_bound: float | None = None,
    upper_correlation_bound: float | None = None,
) -> DuplicateCorrelationFeatureEstimate:
    return DuplicateCorrelationFeatureEstimate(
        feature_id=feature_id,
        status=status,
        failure_reason=failure_reason,
        observed_value_count=observed_value_count,
        design_rank=design_rank,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
        lower_correlation_bound=lower_correlation_bound,
        upper_correlation_bound=upper_correlation_bound,
    )


def _as_float_matrix(
    values: npt.ArrayLike, *, field_name: str
) -> npt.NDArray[np.float64]:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise PhosPyInputError(f"{field_name} must be numeric") from error
    if array.ndim != 2:
        raise PhosPyInputError(f"{field_name} must be two-dimensional")
    return array


def _coerce_feature_ids(
    feature_ids: Sequence[object] | None,
    *,
    feature_count: int,
) -> tuple[str, ...]:
    if feature_ids is None:
        return tuple(f"feature_{index + 1}" for index in range(feature_count))
    if isinstance(feature_ids, (str, bytes, bytearray)):
        raise PhosPyInputError("duplicate_correlation.feature_ids must be a sequence")
    coerced = tuple(
        _require_non_empty_string(
            value,
            field_name="duplicate_correlation.feature_ids[]",
        )
        for value in feature_ids
    )
    if len(coerced) != feature_count:
        raise PhosPyInputError(
            "duplicate_correlation.feature_ids length must match matrix rows; "
            f"feature_ids={len(coerced)}, matrix_rows={feature_count}"
        )
    if len(set(coerced)) != len(coerced):
        raise PhosPyInputError("duplicate_correlation.feature_ids must be unique")
    return coerced


def _coerce_block_ids(
    block_ids: Sequence[object],
    *,
    expected_count: int,
) -> tuple[str, ...]:
    if isinstance(block_ids, (str, bytes, bytearray)):
        raise PhosPyInputError("duplicate_correlation.block_ids must be a sequence")
    coerced = tuple(
        _require_non_empty_string(
            value,
            field_name="duplicate_correlation.block_ids[]",
        )
        for value in block_ids
    )
    if len(coerced) != expected_count:
        raise PhosPyInputError(
            "duplicate_correlation.block_ids length must match design rows; "
            f"block_ids={len(coerced)}, design_rows={expected_count}"
        )
    return coerced


def _reject_block_design_columns(
    design_column_names: Sequence[object] | None,
    *,
    coefficient_count: int,
) -> None:
    if design_column_names is None:
        return
    if isinstance(design_column_names, (str, bytes, bytearray)):
        raise PhosPyInputError(
            "duplicate_correlation.design_column_names must be a sequence"
        )
    column_names = tuple(
        _require_non_empty_string(
            value,
            field_name="duplicate_correlation.design_column_names[]",
        )
        for value in design_column_names
    )
    if len(column_names) != coefficient_count:
        raise PhosPyInputError(
            "duplicate_correlation.design_column_names length must match design "
            "columns; "
            f"design_column_names={len(column_names)}, "
            f"design_columns={coefficient_count}"
        )
    block_columns = tuple(
        column_name
        for column_name in column_names
        if column_name.lower().startswith("block[")
        or column_name.lower().startswith("block_")
    )
    if block_columns:
        raise PhosPyInputError(
            "duplicate_correlation.design must exclude fixed block dummy columns "
            "when using a block-correlation structure; "
            f"block_design_columns={list(block_columns)}"
        )


def _summarize_blocks(
    block_ids: tuple[str, ...],
    *,
    block_id_field_name: str,
) -> DuplicateCorrelationBlockStructureSummary:
    counts = Counter(block_ids)
    levels = tuple(sorted(counts))
    repeated_count = sum(count > 1 for count in counts.values())
    singleton_count = sum(count == 1 for count in counts.values())
    correlated_pair_count = sum(
        (count * (count - 1)) // 2 for count in counts.values() if count > 1
    )
    return DuplicateCorrelationBlockStructureSummary(
        block_id_field_name=block_id_field_name,
        sample_count=len(block_ids),
        block_count=len(levels),
        repeated_block_count=repeated_count,
        singleton_block_count=singleton_count,
        correlated_pair_count=correlated_pair_count,
        block_levels=levels,
        minimum_block_size=min(counts.values()),
        maximum_block_size=max(counts.values()),
    )


def _max_block_size(block_ids: tuple[str, ...]) -> int:
    if not block_ids:
        return 0
    return max(Counter(block_ids).values())


def _scaled_design_rank(design: npt.NDArray[np.float64]) -> int:
    if design.ndim != 2:
        raise PhosPyInputError("duplicate_correlation.design must be two-dimensional")
    sample_count = int(design.shape[0])
    coefficient_count = int(design.shape[1])
    if sample_count < 1 or coefficient_count < 1:
        return 0
    if not np.isfinite(design).all():
        raise PhosPyInputError(
            "duplicate_correlation.design must contain only finite numeric values"
        )
    column_scales = np.asarray(np.linalg.norm(design, axis=0), dtype=np.float64)
    usable_columns = np.isfinite(column_scales) & (column_scales > 0.0)
    if not np.any(usable_columns):
        return 0
    scaled_design = np.asarray(
        design[:, usable_columns] / column_scales[usable_columns][np.newaxis, :],
        dtype=np.float64,
    )
    try:
        singular_values = np.asarray(
            np.linalg.svd(scaled_design, compute_uv=False),
            dtype=np.float64,
        )
    except np.linalg.LinAlgError as error:
        raise PhosPyInputError(
            "duplicate_correlation.design scaled SVD rank assessment failed"
        ) from error
    if singular_values.size == 0 or not np.isfinite(singular_values).all():
        return 0
    largest = float(singular_values[0])
    tolerance = (
        np.finfo(np.float64).eps * float(max(sample_count, coefficient_count)) * largest
    )
    return int(np.count_nonzero(singular_values > tolerance))


def _ols_residual_sum_of_squares(
    values: npt.NDArray[np.float64],
    design: npt.NDArray[np.float64],
) -> float:
    try:
        coefficients = np.asarray(
            np.linalg.lstsq(design, values, rcond=None)[0],
            dtype=np.float64,
        )
    except np.linalg.LinAlgError:
        return math.inf
    residuals = values - (design @ coefficients)
    return float(residuals @ residuals)


def _residual_variation_tolerance(
    values: npt.NDArray[np.float64],
    *,
    sample_count: int,
    coefficient_count: int,
) -> float:
    scale = max(float(np.max(np.abs(values))), 1.0)
    return (
        _ZERO_RESIDUAL_VARIATION_RELATIVE_TOLERANCE
        * float(max(sample_count, coefficient_count))
        * scale
        * scale
    )


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if value is None:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    if isinstance(value, float) and math.isnan(value):
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    text = str(value).strip()
    if not text or text.lower() == "nan":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return text


def _require_positive_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric")
    coerced = float(value)
    if not math.isfinite(coerced) or coerced <= 0.0:
        raise PhosPyInputError(f"{field_name} must be finite and > 0.0")
    return coerced


def _require_positive_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise PhosPyInputError(f"{field_name} must be a positive integer")
    return int(value)


def _require_correlation(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValueError("duplicate-correlation feature estimate must be numeric")
    correlation = float(value)
    if not math.isfinite(correlation) or not -1.0 < correlation < 1.0:
        raise ValueError(
            "duplicate-correlation feature estimate must be finite and in (-1, 1)"
        )
    return correlation


__all__ = [
    "CompoundSymmetryGLSFit",
    "DUPLICATE_CORRELATION_BOUNDARY_DETECTION_TOLERANCE",
    "DUPLICATE_CORRELATION_FISHER_BOUNDARY_TOLERANCE",
    "DUPLICATE_CORRELATION_OPTIMIZER_ABSOLUTE_TOLERANCE",
    "DUPLICATE_CORRELATION_OPTIMIZER_MAX_ITERATIONS",
    "DUPLICATE_CORRELATION_POSITIVE_DEFINITE_TOLERANCE",
    "estimate_duplicate_correlation_reml_consensus",
    "fit_compound_symmetry_gls",
    "fit_duplicate_correlation_gls",
]
