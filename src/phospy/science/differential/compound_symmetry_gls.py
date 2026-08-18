"""Compound-symmetry GLS fitting for duplicate-correlation workflows.

The functions in this module are internal differential-science utilities. They
consume a fixed-effects design, block identities, and a known consensus
correlation, then return the fit quantities needed before contrast fitting and
empirical-Bayes moderation.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast

import numpy as np
import numpy.typing as npt
from scipy.linalg import solve_triangular

from phospy.errors.input import PhosPyInputError

COMPOUND_SYMMETRY_GLS_SOLVER = "compound_symmetry_cholesky_gls"
COMPOUND_SYMMETRY_GLS_COVARIANCE_MODEL = "block_compound_symmetry_correlation"
COMPOUND_SYMMETRY_GLS_RANK_TOLERANCE_POLICY = (
    "rank = count(singular_value > eps * max(n_observations, n_coefficients) * "
    "largest_singular_value) after L2 column scaling of the whitened design"
)

COMPOUND_SYMMETRY_GLS_STATUS_FIT = "fit"
COMPOUND_SYMMETRY_GLS_STATUS_PARTIAL_RANK_LOSS = "partial_rank_loss"
COMPOUND_SYMMETRY_GLS_STATUS_NO_OBSERVATIONS = "no_observations"
COMPOUND_SYMMETRY_GLS_STATUS_NO_ESTIMABLE_COEFFICIENTS = "no_estimable_coefficients"

_NEGATIVE_VARIANCE_ROUNDOFF_TOLERANCE = 1.0e-12

_FloatArray = npt.NDArray[np.float64]
_BoolArray = npt.NDArray[np.bool_]
_IntArray = npt.NDArray[np.int64]


@dataclass(frozen=True, slots=True, eq=False)
class CompoundSymmetryGLSFit:
    """Internal GLS fit artifact aligned to differential model semantics.

    Numeric arrays use the same orientation as ``DifferentialLinearFit`` where
    possible: coefficient-by-feature for coefficient-level quantities and
    sample-by-feature for fitted values/residuals. Missing observations and
    non-estimable coefficients are represented by ``NaN`` in numeric outputs and
    by explicit estimability/status fields.
    """

    __hash__ = object.__hash__

    feature_ids: tuple[str, ...]
    coefficient_names: tuple[str, ...]
    sample_count: int
    coefficient_count: int
    consensus_correlation: float
    block_ids: tuple[str, ...]
    solver: str
    covariance_model: str
    rank_tolerance_policy: str
    coefficients: _FloatArray
    coefficient_estimability: _BoolArray
    stdev_unscaled: _FloatArray
    residual_standard_deviation: _FloatArray
    residual_variance: _FloatArray
    residual_degrees_of_freedom: _FloatArray
    residual_sum_of_squares: _FloatArray
    fitted_values: _FloatArray
    residuals: _FloatArray
    average_expression: _FloatArray
    observed_value_counts: _IntArray
    feature_ranks: _IntArray
    feature_fit_statuses: tuple[str, ...]
    coefficient_covariance: _FloatArray
    feature_coefficient_covariances: _FloatArray
    factorization_cache_size: int
    factorization_cache_hit_count: int
    contrast_names: tuple[str, ...] = ()
    contrast_coefficients: _FloatArray | None = None
    contrast_stdev_unscaled: _FloatArray | None = None
    contrast_covariance: _FloatArray | None = None
    feature_contrast_covariances: _FloatArray | None = None

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another GLS fit has identical scientific content."""

        if not isinstance(other, CompoundSymmetryGLSFit):
            return False
        return (
            self.feature_ids == other.feature_ids
            and self.coefficient_names == other.coefficient_names
            and self.sample_count == other.sample_count
            and self.coefficient_count == other.coefficient_count
            and self.consensus_correlation == other.consensus_correlation
            and self.block_ids == other.block_ids
            and self.solver == other.solver
            and self.covariance_model == other.covariance_model
            and self.rank_tolerance_policy == other.rank_tolerance_policy
            and np.array_equal(self.coefficients, other.coefficients, equal_nan=True)
            and np.array_equal(
                self.coefficient_estimability,
                other.coefficient_estimability,
            )
            and np.array_equal(
                self.stdev_unscaled,
                other.stdev_unscaled,
                equal_nan=True,
            )
            and np.array_equal(
                self.residual_standard_deviation,
                other.residual_standard_deviation,
                equal_nan=True,
            )
            and np.array_equal(
                self.residual_variance,
                other.residual_variance,
                equal_nan=True,
            )
            and np.array_equal(
                self.residual_degrees_of_freedom,
                other.residual_degrees_of_freedom,
                equal_nan=True,
            )
            and np.array_equal(
                self.residual_sum_of_squares,
                other.residual_sum_of_squares,
                equal_nan=True,
            )
            and np.array_equal(self.fitted_values, other.fitted_values, equal_nan=True)
            and np.array_equal(self.residuals, other.residuals, equal_nan=True)
            and np.array_equal(
                self.average_expression,
                other.average_expression,
                equal_nan=True,
            )
            and np.array_equal(self.observed_value_counts, other.observed_value_counts)
            and np.array_equal(self.feature_ranks, other.feature_ranks)
            and self.feature_fit_statuses == other.feature_fit_statuses
            and np.array_equal(
                self.coefficient_covariance,
                other.coefficient_covariance,
                equal_nan=True,
            )
            and np.array_equal(
                self.feature_coefficient_covariances,
                other.feature_coefficient_covariances,
                equal_nan=True,
            )
            and self.factorization_cache_size == other.factorization_cache_size
            and self.factorization_cache_hit_count
            == other.factorization_cache_hit_count
            and self.contrast_names == other.contrast_names
            and _optional_array_equal(
                self.contrast_coefficients,
                other.contrast_coefficients,
            )
            and _optional_array_equal(
                self.contrast_stdev_unscaled,
                other.contrast_stdev_unscaled,
            )
            and _optional_array_equal(
                self.contrast_covariance,
                other.contrast_covariance,
            )
            and _optional_array_equal(
                self.feature_contrast_covariances,
                other.feature_contrast_covariances,
            )
        )

    def coefficient_covariance_for_feature(
        self,
        feature_position: int,
    ) -> _FloatArray:
        """Return a detached coefficient covariance matrix for one feature."""

        position = _require_position(
            feature_position,
            upper_bound=len(self.feature_ids),
            field_name="compound_symmetry_gls.feature_position",
        )
        return np.array(
            self.feature_coefficient_covariances[position, :, :],
            dtype=np.float64,
            copy=True,
        )

    def contrast_covariance_for_feature(
        self,
        feature_position: int,
    ) -> _FloatArray:
        """Return a detached contrast covariance matrix for one feature."""

        if self.feature_contrast_covariances is None:
            raise PhosPyInputError(
                "compound_symmetry_gls.contrasts were not supplied for this fit"
            )
        position = _require_position(
            feature_position,
            upper_bound=len(self.feature_ids),
            field_name="compound_symmetry_gls.feature_position",
        )
        return np.array(
            self.feature_contrast_covariances[position, :, :],
            dtype=np.float64,
            copy=True,
        )


@dataclass(frozen=True, slots=True)
class _ObservationFactorization:
    mask: tuple[bool, ...]
    positions: _IntArray
    observed_count: int
    rank: int
    residual_degrees_of_freedom: float
    selected_column_positions: tuple[int, ...]
    whitened_design: _FloatArray
    covariance_cholesky: _FloatArray | None
    coefficient_covariance: _FloatArray
    stdev_unscaled: _FloatArray
    status: str


def fit_compound_symmetry_gls(
    matrix: npt.ArrayLike,
    design: npt.ArrayLike,
    block_ids: Sequence[object],
    consensus_correlation: object,
    *,
    feature_ids: Sequence[object] | None = None,
    coefficient_names: Sequence[object] | None = None,
    observation_mask: npt.ArrayLike | None = None,
    contrasts: npt.ArrayLike | None = None,
    contrast_names: Sequence[object] | None = None,
) -> CompoundSymmetryGLSFit:
    """Fit feature-wise GLS under a block compound-symmetry correlation.

    ``matrix`` is feature-by-sample. ``design`` is sample-by-fixed-effect.
    ``block_ids`` is aligned to design rows and matrix columns. Non-finite
    feature observations are omitted feature-by-feature and are never replaced.
    The supplied ``consensus_correlation`` is treated as known.
    """

    matrix_values = _as_float_matrix(
        matrix,
        field_name="compound_symmetry_gls.matrix",
    )
    design_values = _as_float_matrix(
        design,
        field_name="compound_symmetry_gls.design",
    )
    if matrix_values.shape[0] < 1:
        raise PhosPyInputError(
            "compound_symmetry_gls.matrix must contain at least one feature"
        )
    if matrix_values.shape[1] != design_values.shape[0]:
        raise PhosPyInputError(
            "compound_symmetry_gls.matrix column count must match "
            "compound_symmetry_gls.design row count; "
            f"matrix_columns={int(matrix_values.shape[1])}, "
            f"design_rows={int(design_values.shape[0])}"
        )
    if design_values.shape[1] < 1:
        raise PhosPyInputError(
            "compound_symmetry_gls.design must contain at least one fixed effect"
        )
    if not np.isfinite(design_values).all():
        raise PhosPyInputError(
            "compound_symmetry_gls.design must contain only finite numeric values"
        )
    observation_masks = _coerce_observation_mask(
        observation_mask,
        matrix=matrix_values,
    )

    feature_id_values = _coerce_unique_labels(
        feature_ids,
        expected_count=int(matrix_values.shape[0]),
        field_name="compound_symmetry_gls.feature_ids",
        default_prefix="feature",
    )
    coefficient_name_values = _coerce_unique_labels(
        coefficient_names,
        expected_count=int(design_values.shape[1]),
        field_name="compound_symmetry_gls.coefficient_names",
        default_prefix="coefficient",
    )
    _reject_block_design_columns(coefficient_name_values)
    block_id_values = _coerce_block_ids(
        block_ids,
        expected_count=int(design_values.shape[0]),
    )
    correlation = _require_valid_consensus_correlation(
        consensus_correlation,
        block_ids=block_id_values,
    )
    _require_full_design_rank(design_values)

    contrast_values: _FloatArray | None = None
    contrast_name_values: tuple[str, ...] = ()
    if contrasts is not None:
        contrast_values = _coerce_contrasts(
            contrasts,
            coefficient_count=int(design_values.shape[1]),
        )
        contrast_name_values = _coerce_unique_labels(
            contrast_names,
            expected_count=int(contrast_values.shape[1]),
            field_name="compound_symmetry_gls.contrast_names",
            default_prefix="contrast",
        )

    full_factorization = _build_observation_factorization(
        tuple(True for _ in range(int(design_values.shape[0]))),
        design=design_values,
        block_ids=block_id_values,
        consensus_correlation=correlation,
    )
    if full_factorization.rank < int(design_values.shape[1]):
        raise PhosPyInputError(
            "compound_symmetry_gls.design lost rank under the full "
            "compound-symmetry covariance; "
            f"rank={full_factorization.rank}, columns={int(design_values.shape[1])}"
        )
    common_coefficient_correlation = _correlation_from_covariance(
        full_factorization.coefficient_covariance
    )

    feature_count = int(matrix_values.shape[0])
    sample_count = int(matrix_values.shape[1])
    coefficient_count = int(design_values.shape[1])
    contrast_count = 0 if contrast_values is None else int(contrast_values.shape[1])

    coefficients = np.full(
        (coefficient_count, feature_count),
        np.nan,
        dtype=np.float64,
    )
    coefficient_estimability = np.zeros(
        (coefficient_count, feature_count),
        dtype=np.bool_,
    )
    stdev_unscaled = np.full_like(coefficients, np.nan, dtype=np.float64)
    residual_standard_deviation = np.full(feature_count, np.nan, dtype=np.float64)
    residual_variance = np.full(feature_count, np.nan, dtype=np.float64)
    residual_degrees_of_freedom = np.zeros(feature_count, dtype=np.float64)
    residual_sum_of_squares = np.full(feature_count, np.nan, dtype=np.float64)
    fitted_values = np.full((sample_count, feature_count), np.nan, dtype=np.float64)
    residuals = np.full((sample_count, feature_count), np.nan, dtype=np.float64)
    average_expression = np.full(feature_count, np.nan, dtype=np.float64)
    observed_value_counts = np.zeros(feature_count, dtype=np.int64)
    feature_ranks = np.zeros(feature_count, dtype=np.int64)
    feature_covariances = np.full(
        (feature_count, coefficient_count, coefficient_count),
        np.nan,
        dtype=np.float64,
    )
    contrast_coefficients = (
        None
        if contrast_values is None
        else np.full((contrast_count, feature_count), np.nan, dtype=np.float64)
    )
    contrast_stdev_unscaled = (
        None
        if contrast_values is None
        else np.full((contrast_count, feature_count), np.nan, dtype=np.float64)
    )
    feature_contrast_covariances = (
        None
        if contrast_values is None
        else np.full(
            (feature_count, contrast_count, contrast_count),
            np.nan,
            dtype=np.float64,
        )
    )

    factorization_cache: dict[tuple[bool, ...], _ObservationFactorization] = {
        full_factorization.mask: full_factorization
    }
    cache_hits = 0
    feature_statuses: list[str] = []

    for feature_position in range(feature_count):
        values = matrix_values[feature_position, :]
        feature_observation_mask = tuple(
            bool(value) for value in observation_masks[feature_position, :].tolist()
        )
        factorization = factorization_cache.get(feature_observation_mask)
        if factorization is None:
            factorization = _build_observation_factorization(
                feature_observation_mask,
                design=design_values,
                block_ids=block_id_values,
                consensus_correlation=correlation,
            )
            factorization_cache[feature_observation_mask] = factorization
        else:
            cache_hits += 1

        status = factorization.status
        feature_statuses.append(status)
        observed_value_counts[feature_position] = factorization.observed_count
        feature_ranks[feature_position] = factorization.rank
        residual_degrees_of_freedom[feature_position] = (
            factorization.residual_degrees_of_freedom
        )
        feature_covariances[feature_position, :, :] = (
            factorization.coefficient_covariance
        )
        stdev_unscaled[:, feature_position] = factorization.stdev_unscaled

        if factorization.observed_count > 0:
            observed_values = values[factorization.positions]
            average_expression[feature_position] = float(np.mean(observed_values))
        if (
            factorization.status == COMPOUND_SYMMETRY_GLS_STATUS_NO_OBSERVATIONS
            or factorization.status
            == COMPOUND_SYMMETRY_GLS_STATUS_NO_ESTIMABLE_COEFFICIENTS
        ):
            continue
        if factorization.covariance_cholesky is None:
            raise PhosPyInputError(
                "compound_symmetry_gls internal factorization is missing the "
                "covariance Cholesky factor for an observed feature"
            )

        observed_values = values[factorization.positions]
        selected_columns = factorization.selected_column_positions
        whitened_values = _whiten_values(
            factorization.covariance_cholesky,
            observed_values,
        )
        selected_design = factorization.whitened_design[:, list(selected_columns)]
        selected_covariance = factorization.coefficient_covariance[
            np.ix_(selected_columns, selected_columns)
        ]
        beta_selected = cast(
            _FloatArray,
            np.asarray(
                selected_covariance @ (selected_design.T @ whitened_values),
                dtype=np.float64,
            ),
        )
        selected_array = np.asarray(selected_columns, dtype=np.int64)
        coefficients[selected_array, feature_position] = beta_selected
        coefficient_estimability[selected_array, feature_position] = True

        observed_design = design_values[factorization.positions, :]
        raw_fitted = cast(
            _FloatArray,
            np.asarray(
                observed_design[:, list(selected_columns)] @ beta_selected,
                dtype=np.float64,
            ),
        )
        raw_residuals = cast(
            _FloatArray,
            np.asarray(observed_values - raw_fitted, dtype=np.float64),
        )
        raw_residuals = _zero_residual_roundoff(
            residuals=raw_residuals,
            response_values=observed_values,
            fitted_values=raw_fitted,
            sample_count=factorization.observed_count,
            coefficient_count=coefficient_count,
        )
        whitened_residuals = _whiten_values(
            factorization.covariance_cholesky,
            raw_residuals,
        )
        rss = float(whitened_residuals @ whitened_residuals)
        if not math.isfinite(rss) or rss < 0.0:
            raise PhosPyInputError(
                "compound_symmetry_gls produced an invalid whitened residual sum "
                f"of squares for feature {feature_id_values[feature_position]!r}"
            )
        residual_sum_of_squares[feature_position] = rss
        if factorization.residual_degrees_of_freedom > 0.0:
            variance = rss / factorization.residual_degrees_of_freedom
            residual_variance[feature_position] = variance
            residual_standard_deviation[feature_position] = math.sqrt(
                max(variance, 0.0)
            )

        fitted_values[factorization.positions, feature_position] = raw_fitted
        residuals[factorization.positions, feature_position] = raw_residuals

        if contrast_values is not None:
            assert contrast_coefficients is not None
            assert contrast_stdev_unscaled is not None
            assert feature_contrast_covariances is not None
            _fill_feature_contrast_fit(
                contrast_values=contrast_values,
                beta=coefficients[:, feature_position],
                coefficient_estimability=coefficient_estimability[:, feature_position],
                coefficient_covariance=factorization.coefficient_covariance,
                common_coefficient_correlation=common_coefficient_correlation,
                output_coefficients=contrast_coefficients[:, feature_position],
                output_stdev_unscaled=contrast_stdev_unscaled[:, feature_position],
                output_covariance=feature_contrast_covariances[
                    feature_position,
                    :,
                    :,
                ],
            )

    common_contrast_covariance = None
    if contrast_values is not None:
        common_contrast_covariance = _contrast_covariance(
            contrasts=contrast_values,
            coefficient_covariance=full_factorization.coefficient_covariance,
        )

    return CompoundSymmetryGLSFit(
        feature_ids=feature_id_values,
        coefficient_names=coefficient_name_values,
        sample_count=sample_count,
        coefficient_count=coefficient_count,
        consensus_correlation=correlation,
        block_ids=block_id_values,
        solver=COMPOUND_SYMMETRY_GLS_SOLVER,
        covariance_model=COMPOUND_SYMMETRY_GLS_COVARIANCE_MODEL,
        rank_tolerance_policy=COMPOUND_SYMMETRY_GLS_RANK_TOLERANCE_POLICY,
        coefficients=_immutable_float_array(coefficients),
        coefficient_estimability=_immutable_bool_array(coefficient_estimability),
        stdev_unscaled=_immutable_float_array(stdev_unscaled),
        residual_standard_deviation=_immutable_float_array(residual_standard_deviation),
        residual_variance=_immutable_float_array(residual_variance),
        residual_degrees_of_freedom=_immutable_float_array(residual_degrees_of_freedom),
        residual_sum_of_squares=_immutable_float_array(residual_sum_of_squares),
        fitted_values=_immutable_float_array(fitted_values),
        residuals=_immutable_float_array(residuals),
        average_expression=_immutable_float_array(average_expression),
        observed_value_counts=_immutable_int_array(observed_value_counts),
        feature_ranks=_immutable_int_array(feature_ranks),
        feature_fit_statuses=tuple(feature_statuses),
        coefficient_covariance=_immutable_float_array(
            full_factorization.coefficient_covariance
        ),
        feature_coefficient_covariances=_immutable_float_array(feature_covariances),
        factorization_cache_size=len(factorization_cache),
        factorization_cache_hit_count=cache_hits,
        contrast_names=contrast_name_values,
        contrast_coefficients=(
            None
            if contrast_coefficients is None
            else _immutable_float_array(contrast_coefficients)
        ),
        contrast_stdev_unscaled=(
            None
            if contrast_stdev_unscaled is None
            else _immutable_float_array(contrast_stdev_unscaled)
        ),
        contrast_covariance=(
            None
            if common_contrast_covariance is None
            else _immutable_float_array(common_contrast_covariance)
        ),
        feature_contrast_covariances=(
            None
            if feature_contrast_covariances is None
            else _immutable_float_array(feature_contrast_covariances)
        ),
    )


def fit_duplicate_correlation_gls(
    matrix: npt.ArrayLike,
    design: npt.ArrayLike,
    block_ids: Sequence[object],
    consensus_correlation: object,
    *,
    feature_ids: Sequence[object] | None = None,
    coefficient_names: Sequence[object] | None = None,
    observation_mask: npt.ArrayLike | None = None,
    contrasts: npt.ArrayLike | None = None,
    contrast_names: Sequence[object] | None = None,
) -> CompoundSymmetryGLSFit:
    """Alias for duplicate-correlation GLS fitting."""

    return fit_compound_symmetry_gls(
        matrix,
        design,
        block_ids,
        consensus_correlation,
        feature_ids=feature_ids,
        coefficient_names=coefficient_names,
        observation_mask=observation_mask,
        contrasts=contrasts,
        contrast_names=contrast_names,
    )


def _build_observation_factorization(
    observation_mask: tuple[bool, ...],
    *,
    design: _FloatArray,
    block_ids: tuple[str, ...],
    consensus_correlation: float,
) -> _ObservationFactorization:
    positions = np.asarray(
        [position for position, observed in enumerate(observation_mask) if observed],
        dtype=np.int64,
    )
    observed_count = int(positions.size)
    coefficient_count = int(design.shape[1])
    empty_covariance = np.full(
        (coefficient_count, coefficient_count),
        np.nan,
        dtype=np.float64,
    )
    empty_stdev = np.full(coefficient_count, np.nan, dtype=np.float64)
    if observed_count == 0:
        return _ObservationFactorization(
            mask=observation_mask,
            positions=_immutable_int_array(positions),
            observed_count=0,
            rank=0,
            residual_degrees_of_freedom=0.0,
            selected_column_positions=(),
            whitened_design=_immutable_float_array(
                np.empty((0, coefficient_count), dtype=np.float64)
            ),
            covariance_cholesky=None,
            coefficient_covariance=_immutable_float_array(empty_covariance),
            stdev_unscaled=_immutable_float_array(empty_stdev),
            status=COMPOUND_SYMMETRY_GLS_STATUS_NO_OBSERVATIONS,
        )

    observed_design = cast(
        _FloatArray,
        np.asarray(design[positions, :], dtype=np.float64),
    )
    observed_blocks = tuple(block_ids[int(position)] for position in positions.tolist())
    covariance = _compound_symmetry_correlation_matrix(
        consensus_correlation,
        block_ids=observed_blocks,
    )
    try:
        covariance_cholesky = cast(
            _FloatArray,
            np.asarray(np.linalg.cholesky(covariance), dtype=np.float64),
        )
    except np.linalg.LinAlgError as error:
        raise PhosPyInputError(
            "compound_symmetry_gls covariance Cholesky decomposition failed for "
            "an observation subset; the block correlation covariance is not "
            "positive definite"
        ) from error

    whitened_design = _whiten_matrix(covariance_cholesky, observed_design)
    selected_column_positions, rank = _independent_column_positions(whitened_design)
    residual_degrees_of_freedom = float(observed_count - rank)
    if rank == 0:
        return _ObservationFactorization(
            mask=observation_mask,
            positions=_immutable_int_array(positions),
            observed_count=observed_count,
            rank=0,
            residual_degrees_of_freedom=float(observed_count),
            selected_column_positions=(),
            whitened_design=_immutable_float_array(whitened_design),
            covariance_cholesky=_immutable_float_array(covariance_cholesky),
            coefficient_covariance=_immutable_float_array(empty_covariance),
            stdev_unscaled=_immutable_float_array(empty_stdev),
            status=COMPOUND_SYMMETRY_GLS_STATUS_NO_ESTIMABLE_COEFFICIENTS,
        )

    selected_design = whitened_design[:, list(selected_column_positions)]
    coefficient_covariance = np.full(
        (coefficient_count, coefficient_count),
        np.nan,
        dtype=np.float64,
    )
    selected_covariance = _coefficient_covariance_from_full_rank_design(selected_design)
    coefficient_covariance[
        np.ix_(selected_column_positions, selected_column_positions)
    ] = selected_covariance
    stdev_unscaled = _stdev_from_covariance_diagonal(coefficient_covariance)
    return _ObservationFactorization(
        mask=observation_mask,
        positions=_immutable_int_array(positions),
        observed_count=observed_count,
        rank=rank,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
        selected_column_positions=selected_column_positions,
        whitened_design=_immutable_float_array(whitened_design),
        covariance_cholesky=_immutable_float_array(covariance_cholesky),
        coefficient_covariance=_immutable_float_array(coefficient_covariance),
        stdev_unscaled=_immutable_float_array(stdev_unscaled),
        status=(
            COMPOUND_SYMMETRY_GLS_STATUS_FIT
            if rank == coefficient_count
            else COMPOUND_SYMMETRY_GLS_STATUS_PARTIAL_RANK_LOSS
        ),
    )


def _coefficient_covariance_from_full_rank_design(
    design: _FloatArray,
) -> _FloatArray:
    gram = cast(_FloatArray, np.asarray(design.T @ design, dtype=np.float64))
    try:
        gram_cholesky = cast(
            _FloatArray,
            np.asarray(np.linalg.cholesky(gram), dtype=np.float64),
        )
    except np.linalg.LinAlgError as error:
        raise PhosPyInputError(
            "compound_symmetry_gls coefficient covariance Cholesky decomposition "
            "failed for the whitened fixed-effects design"
        ) from error
    identity = np.eye(int(gram.shape[0]), dtype=np.float64)
    return _cholesky_solve(gram_cholesky, identity)


def _fill_feature_contrast_fit(
    *,
    contrast_values: _FloatArray,
    beta: _FloatArray,
    coefficient_estimability: _BoolArray,
    coefficient_covariance: _FloatArray,
    common_coefficient_correlation: _FloatArray,
    output_coefficients: _FloatArray,
    output_stdev_unscaled: _FloatArray,
    output_covariance: _FloatArray,
) -> None:
    estimable_contrasts = _estimable_contrast_mask(
        contrast_values=contrast_values,
        coefficient_estimability=coefficient_estimability,
    )
    finite_beta = np.where(np.isfinite(beta), beta, 0.0)
    contrast_estimates = cast(
        _FloatArray,
        np.asarray(contrast_values.T @ finite_beta, dtype=np.float64),
    )
    contrast_covariance = _contrast_covariance(
        contrasts=contrast_values,
        coefficient_covariance=_limma_style_feature_coefficient_covariance(
            coefficient_covariance=coefficient_covariance,
            common_coefficient_correlation=common_coefficient_correlation,
        ),
    )
    stdev = _stdev_from_covariance_diagonal(contrast_covariance)
    output_covariance[:, :] = contrast_covariance
    non_estimable_contrasts = ~estimable_contrasts
    output_covariance[non_estimable_contrasts, :] = np.nan
    output_covariance[:, non_estimable_contrasts] = np.nan
    output_coefficients[estimable_contrasts] = contrast_estimates[estimable_contrasts]
    output_stdev_unscaled[estimable_contrasts] = stdev[estimable_contrasts]


def _estimable_contrast_mask(
    *,
    contrast_values: _FloatArray,
    coefficient_estimability: _BoolArray,
) -> _BoolArray:
    nonzero = np.abs(contrast_values) > 0.0
    unavailable = nonzero & ~coefficient_estimability[:, np.newaxis]
    return cast(_BoolArray, np.asarray(~np.any(unavailable, axis=0), dtype=np.bool_))


def _contrast_covariance(
    *,
    contrasts: _FloatArray,
    coefficient_covariance: _FloatArray,
) -> _FloatArray:
    finite_covariance = np.where(
        np.isfinite(coefficient_covariance),
        coefficient_covariance,
        0.0,
    )
    covariance = cast(
        _FloatArray,
        np.asarray(contrasts.T @ finite_covariance @ contrasts, dtype=np.float64),
    )
    covariance[~np.isfinite(covariance)] = np.nan
    return covariance


def _limma_style_feature_coefficient_covariance(
    *,
    coefficient_covariance: _FloatArray,
    common_coefficient_correlation: _FloatArray,
) -> _FloatArray:
    feature_stdev = _stdev_from_covariance_diagonal(coefficient_covariance)
    finite_stdev = np.where(np.isfinite(feature_stdev), feature_stdev, 0.0)
    covariance = (
        finite_stdev[:, np.newaxis]
        * common_coefficient_correlation
        * finite_stdev[np.newaxis, :]
    )
    covariance[~np.isfinite(covariance)] = np.nan
    return cast(_FloatArray, covariance)


def _correlation_from_covariance(covariance: _FloatArray) -> _FloatArray:
    diagonal = cast(
        _FloatArray,
        np.asarray(np.diagonal(covariance), dtype=np.float64),
    )
    if not np.isfinite(diagonal).all() or np.any(diagonal <= 0.0):
        raise PhosPyInputError(
            "compound_symmetry_gls common coefficient covariance must have "
            "finite positive diagonal entries"
        )
    scale = np.sqrt(diagonal)
    correlation = covariance / (scale[:, np.newaxis] * scale[np.newaxis, :])
    if not np.isfinite(correlation).all():
        raise PhosPyInputError(
            "compound_symmetry_gls common coefficient correlation is invalid"
        )
    return cast(_FloatArray, np.asarray(correlation, dtype=np.float64))


def _stdev_from_covariance_diagonal(covariance: _FloatArray) -> _FloatArray:
    diagonal = cast(
        _FloatArray,
        np.asarray(np.diagonal(covariance), dtype=np.float64),
    )
    stdev = np.full(diagonal.shape, np.nan, dtype=np.float64)
    finite = np.isfinite(diagonal)
    invalid_negative = finite & (diagonal < -_NEGATIVE_VARIANCE_ROUNDOFF_TOLERANCE)
    valid = finite & ~invalid_negative
    stdev[valid] = np.sqrt(np.clip(diagonal[valid], 0.0, None))
    return cast(_FloatArray, stdev)


def _compound_symmetry_correlation_matrix(
    correlation: float,
    *,
    block_ids: tuple[str, ...],
) -> _FloatArray:
    sample_count = len(block_ids)
    covariance = np.eye(sample_count, dtype=np.float64)
    for row in range(sample_count):
        for column in range(row + 1, sample_count):
            if block_ids[row] == block_ids[column]:
                covariance[row, column] = correlation
                covariance[column, row] = correlation
    return covariance


def _require_valid_consensus_correlation(
    value: object,
    *,
    block_ids: tuple[str, ...],
) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(
            "compound_symmetry_gls.consensus_correlation must be numeric"
        )
    correlation = float(value)
    if not math.isfinite(correlation) or not -1.0 < correlation < 1.0:
        raise PhosPyInputError(
            "compound_symmetry_gls.consensus_correlation must be finite and in "
            "(-1.0, 1.0)"
        )
    maximum_repeated_block_size = _maximum_repeated_block_size(block_ids)
    if maximum_repeated_block_size > 1:
        lower_bound = -1.0 / float(maximum_repeated_block_size - 1)
        if correlation <= lower_bound:
            raise PhosPyInputError(
                "compound_symmetry_gls.consensus_correlation is outside the "
                "positive-definite compound-symmetry interval for the supplied "
                "block structure; "
                f"correlation={correlation:.12g}, "
                f"lower_bound={lower_bound:.12g}, upper_bound=1"
            )
    covariance = _compound_symmetry_correlation_matrix(
        correlation,
        block_ids=block_ids,
    )
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as error:
        raise PhosPyInputError(
            "compound_symmetry_gls.consensus_correlation does not produce a "
            "positive-definite block covariance matrix"
        ) from error
    return correlation


def _maximum_repeated_block_size(block_ids: tuple[str, ...]) -> int:
    repeated_sizes = [count for count in Counter(block_ids).values() if count > 1]
    if not repeated_sizes:
        return 1
    return max(repeated_sizes)


def _require_full_design_rank(design: _FloatArray) -> None:
    rank = _scaled_design_rank(design)
    columns = int(design.shape[1])
    if rank < columns:
        raise PhosPyInputError(
            "compound_symmetry_gls.design is rank deficient under scaled SVD "
            "rank assessment; "
            f"rank={rank}, columns={columns}"
        )


def _independent_column_positions(design: _FloatArray) -> tuple[tuple[int, ...], int]:
    selected: list[int] = []
    current_rank = 0
    for column_position in range(int(design.shape[1])):
        candidate = (*selected, column_position)
        candidate_rank = _scaled_design_rank(design[:, list(candidate)])
        if candidate_rank > current_rank:
            selected.append(column_position)
            current_rank = candidate_rank
    return tuple(selected), current_rank


def _scaled_design_rank(design: _FloatArray) -> int:
    if design.ndim != 2:
        raise PhosPyInputError("compound_symmetry_gls.design must be two-dimensional")
    sample_count = int(design.shape[0])
    coefficient_count = int(design.shape[1])
    if sample_count < 1 or coefficient_count < 1:
        return 0
    if not np.isfinite(design).all():
        raise PhosPyInputError(
            "compound_symmetry_gls.design must contain only finite numeric values"
        )
    column_scales = cast(
        _FloatArray,
        np.asarray(np.linalg.norm(design, axis=0), dtype=np.float64),
    )
    usable_columns = np.isfinite(column_scales) & (column_scales > 0.0)
    if not np.any(usable_columns):
        return 0
    scaled_design = cast(
        _FloatArray,
        np.asarray(
            design[:, usable_columns] / column_scales[usable_columns][np.newaxis, :],
            dtype=np.float64,
        ),
    )
    try:
        singular_values = cast(
            _FloatArray,
            np.asarray(
                np.linalg.svd(scaled_design, compute_uv=False),
                dtype=np.float64,
            ),
        )
    except np.linalg.LinAlgError as error:
        raise PhosPyInputError(
            "compound_symmetry_gls.design scaled SVD rank assessment failed"
        ) from error
    if singular_values.size == 0 or not np.isfinite(singular_values).all():
        return 0
    largest = float(singular_values[0])
    tolerance = (
        np.finfo(np.float64).eps * float(max(sample_count, coefficient_count)) * largest
    )
    return int(np.count_nonzero(singular_values > tolerance))


def _as_float_matrix(values: npt.ArrayLike, *, field_name: str) -> _FloatArray:
    try:
        array = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise PhosPyInputError(f"{field_name} must be numeric") from error
    if array.ndim != 2:
        raise PhosPyInputError(f"{field_name} must be two-dimensional")
    return cast(_FloatArray, array)


def _coerce_unique_labels(
    labels: Sequence[object] | None,
    *,
    expected_count: int,
    field_name: str,
    default_prefix: str,
) -> tuple[str, ...]:
    if labels is None:
        return tuple(f"{default_prefix}_{index + 1}" for index in range(expected_count))
    if isinstance(labels, (str, bytes, bytearray)):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    coerced = tuple(
        _require_non_empty_string(value, field_name=f"{field_name}[]")
        for value in labels
    )
    if len(coerced) != expected_count:
        raise PhosPyInputError(
            f"{field_name} length must match expected count; "
            f"labels={len(coerced)}, expected={expected_count}"
        )
    if len(set(coerced)) != len(coerced):
        raise PhosPyInputError(f"{field_name} must be unique")
    return coerced


def _coerce_block_ids(
    block_ids: Sequence[object],
    *,
    expected_count: int,
) -> tuple[str, ...]:
    if isinstance(block_ids, (str, bytes, bytearray)):
        raise PhosPyInputError("compound_symmetry_gls.block_ids must be a sequence")
    coerced = tuple(
        _require_non_empty_string(
            value,
            field_name="compound_symmetry_gls.block_ids[]",
        )
        for value in block_ids
    )
    if len(coerced) != expected_count:
        raise PhosPyInputError(
            "compound_symmetry_gls.block_ids length must match design rows; "
            f"block_ids={len(coerced)}, design_rows={expected_count}"
        )
    return coerced


def _reject_block_design_columns(coefficient_names: tuple[str, ...]) -> None:
    block_columns = tuple(
        column_name
        for column_name in coefficient_names
        if column_name.lower().startswith("block[")
        or column_name.lower().startswith("block_")
    )
    if block_columns:
        raise PhosPyInputError(
            "compound_symmetry_gls.design must exclude fixed block dummy columns "
            "when using a block-correlation structure; "
            f"block_design_columns={list(block_columns)}"
        )


def _coerce_contrasts(
    contrasts: npt.ArrayLike,
    *,
    coefficient_count: int,
) -> _FloatArray:
    values = _as_float_matrix(
        contrasts,
        field_name="compound_symmetry_gls.contrasts",
    )
    if int(values.shape[0]) != coefficient_count:
        raise PhosPyInputError(
            "compound_symmetry_gls.contrasts row count must match design "
            "coefficients; "
            f"contrast_rows={int(values.shape[0])}, "
            f"design_coefficients={coefficient_count}"
        )
    if int(values.shape[1]) < 1:
        raise PhosPyInputError(
            "compound_symmetry_gls.contrasts must contain at least one contrast"
        )
    if not np.isfinite(values).all():
        raise PhosPyInputError(
            "compound_symmetry_gls.contrasts must contain only finite numeric values"
        )
    return values


def _coerce_observation_mask(
    observation_mask: npt.ArrayLike | None,
    *,
    matrix: _FloatArray,
) -> _BoolArray:
    if observation_mask is None:
        return cast(_BoolArray, np.asarray(np.isfinite(matrix), dtype=np.bool_))
    mask = np.asarray(observation_mask, dtype=np.bool_)
    if mask.ndim != 2:
        raise PhosPyInputError(
            "compound_symmetry_gls.observation_mask must be two-dimensional"
        )
    expected_shape = (int(matrix.shape[0]), int(matrix.shape[1]))
    actual_shape = (int(mask.shape[0]), int(mask.shape[1]))
    if actual_shape != expected_shape:
        raise PhosPyInputError(
            "compound_symmetry_gls.observation_mask shape must match matrix shape; "
            f"observation_mask_shape={actual_shape}, matrix_shape={expected_shape}"
        )
    observed_non_finite = mask & ~np.isfinite(matrix)
    if np.any(observed_non_finite):
        positions = np.argwhere(observed_non_finite)
        preview = ", ".join(
            f"({int(row)}, {int(column)})" for row, column in positions[:3]
        )
        suffix = "" if int(positions.shape[0]) <= 3 else ", ..."
        raise PhosPyInputError(
            "compound_symmetry_gls.observation_mask marks non-finite matrix "
            "values as observed; positions="
            f"{preview}{suffix}"
        )
    return cast(_BoolArray, mask)


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if value is None:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    if isinstance(value, float) and math.isnan(value):
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    text = str(value).strip()
    if not text or text.lower() == "nan":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return text


def _require_position(value: object, *, upper_bound: int, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an integer position")
    position = int(value)
    if position < 0 or position >= upper_bound:
        raise PhosPyInputError(f"{field_name} must be within [0, {upper_bound - 1}]")
    return position


def _whiten_matrix(cholesky_factor: _FloatArray, values: _FloatArray) -> _FloatArray:
    return cast(
        _FloatArray,
        solve_triangular(
            cholesky_factor,
            values,
            lower=True,
            check_finite=False,
        ),
    )


def _whiten_values(cholesky_factor: _FloatArray, values: _FloatArray) -> _FloatArray:
    return cast(
        _FloatArray,
        solve_triangular(
            cholesky_factor,
            values,
            lower=True,
            check_finite=False,
        ),
    )


def _cholesky_solve(cholesky_factor: _FloatArray, rhs: _FloatArray) -> _FloatArray:
    forward = cast(
        _FloatArray,
        solve_triangular(
            cholesky_factor,
            rhs,
            lower=True,
            check_finite=False,
        ),
    )
    return cast(
        _FloatArray,
        solve_triangular(
            cholesky_factor.T,
            forward,
            lower=False,
            check_finite=False,
        ),
    )


def _zero_residual_roundoff(
    *,
    residuals: _FloatArray,
    response_values: _FloatArray,
    fitted_values: _FloatArray,
    sample_count: int,
    coefficient_count: int,
) -> _FloatArray:
    response_scale = float(np.nanmax(np.abs(response_values)))
    fitted_scale = float(np.nanmax(np.abs(fitted_values)))
    scale = max(response_scale, fitted_scale, 1.0)
    tolerance = (
        np.finfo(np.float64).eps * float(max(sample_count, coefficient_count)) * scale
    )
    zero_mask = np.isfinite(residuals) & (np.abs(residuals) <= tolerance)
    return cast(_FloatArray, np.where(zero_mask, 0.0, residuals))


def _optional_array_equal(
    left: _FloatArray | None,
    right: _FloatArray | None,
) -> bool:
    if left is None or right is None:
        return left is right
    return bool(np.array_equal(left, right, equal_nan=True))


def _immutable_float_array(values: npt.ArrayLike) -> _FloatArray:
    array = np.array(values, dtype=np.float64, copy=True, order="C")
    immutable_buffer = array.tobytes(order="C")
    immutable_view = np.frombuffer(
        immutable_buffer,
        dtype=np.float64,
        count=int(array.size),
    ).reshape(array.shape, order="C")
    return cast(_FloatArray, immutable_view)


def _immutable_bool_array(values: npt.ArrayLike) -> _BoolArray:
    array = np.array(values, dtype=np.bool_, copy=True, order="C")
    immutable_buffer = array.tobytes(order="C")
    immutable_view = np.frombuffer(
        immutable_buffer,
        dtype=np.bool_,
        count=int(array.size),
    ).reshape(array.shape, order="C")
    return cast(_BoolArray, immutable_view)


def _immutable_int_array(values: npt.ArrayLike) -> _IntArray:
    array = np.array(values, dtype=np.int64, copy=True, order="C")
    immutable_buffer = array.tobytes(order="C")
    immutable_view = np.frombuffer(
        immutable_buffer,
        dtype=np.int64,
        count=int(array.size),
    ).reshape(array.shape, order="C")
    return cast(_IntArray, immutable_view)


__all__ = [
    "COMPOUND_SYMMETRY_GLS_COVARIANCE_MODEL",
    "COMPOUND_SYMMETRY_GLS_RANK_TOLERANCE_POLICY",
    "COMPOUND_SYMMETRY_GLS_SOLVER",
    "COMPOUND_SYMMETRY_GLS_STATUS_FIT",
    "COMPOUND_SYMMETRY_GLS_STATUS_NO_ESTIMABLE_COEFFICIENTS",
    "COMPOUND_SYMMETRY_GLS_STATUS_NO_OBSERVATIONS",
    "COMPOUND_SYMMETRY_GLS_STATUS_PARTIAL_RANK_LOSS",
    "CompoundSymmetryGLSFit",
    "fit_compound_symmetry_gls",
    "fit_duplicate_correlation_gls",
]
