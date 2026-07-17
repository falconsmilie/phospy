"""Stable linear-model decomposition for differential analysis."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import cast

import numpy as np
import numpy.typing as npt

DIFFERENTIAL_LINEAR_MODEL_DECOMPOSITION_METHOD = "scaled_svd"
DIFFERENTIAL_LINEAR_MODEL_SOLVER = "scaled_svd_least_squares"
DIFFERENTIAL_LINEAR_MODEL_COLUMN_SCALE_METHOD = "l2_norm"
DIFFERENTIAL_LINEAR_MODEL_RANK_TOLERANCE_POLICY = (
    "rank = count(singular_value > eps * max(n_samples, n_coefficients) * "
    "largest_singular_value) after L2 column scaling"
)
DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER = 1.0e10

_NEGATIVE_VARIANCE_ROUNDOFF_TOLERANCE = 1.0e-12


class DifferentialDesignDecompositionError(ValueError):
    """Raised when a differential design cannot support stable OLS fitting."""


@dataclass(frozen=True, slots=True)
class DifferentialLinearFit:
    """Feature-wise OLS fit returned by the shared design decomposition."""

    coefficients: npt.NDArray[np.float64]
    fitted_values: npt.NDArray[np.float64]
    residuals: npt.NDArray[np.float64]
    residual_sum_of_squares: npt.NDArray[np.float64]
    residual_variance: npt.NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class DifferentialDesignDecomposition:
    """Scaled-SVD decomposition owning differential linear-model numerics."""

    sample_count: int
    coefficient_count: int
    rank: int
    residual_degrees_of_freedom: float
    singular_values: tuple[float, ...]
    rank_tolerance: float
    condition_number: float
    max_condition_number: float
    decomposition_method: str
    solver: str
    column_scale_method: str
    rank_tolerance_policy: str
    _design_values: npt.NDArray[np.float64]
    _column_scales: npt.NDArray[np.float64]
    _u: npt.NDArray[np.float64]
    _vt: npt.NDArray[np.float64]
    _coefficient_covariance: npt.NDArray[np.float64]

    @property
    def coefficient_covariance(self) -> npt.NDArray[np.float64]:
        """Return the OLS coefficient covariance factor in original column units."""

        return np.array(self._coefficient_covariance, dtype=np.float64, copy=True)

    def assert_matches_design(
        self,
        design: npt.ArrayLike,
        *,
        field_name: str = "differential.design",
    ) -> None:
        """Raise if this decomposition was not built from the supplied design."""

        design_values = np.asarray(design, dtype=np.float64)
        if design_values.ndim != 2:
            raise DifferentialDesignDecompositionError(
                f"{field_name} must be two-dimensional"
            )
        expected_shape = (self.sample_count, self.coefficient_count)
        actual_shape = (int(design_values.shape[0]), int(design_values.shape[1]))
        if actual_shape != expected_shape:
            raise DifferentialDesignDecompositionError(
                f"{field_name} shape does not match differential design "
                "decomposition; "
                f"design_shape={actual_shape}, decomposition_shape={expected_shape}"
            )
        if not np.isfinite(design_values).all():
            raise DifferentialDesignDecompositionError(
                f"{field_name} must contain only finite numeric values"
            )
        if not np.array_equal(design_values, self._design_values):
            raise DifferentialDesignDecompositionError(
                f"{field_name} values do not match differential design decomposition"
            )

    def fit(self, response: npt.ArrayLike) -> DifferentialLinearFit:
        """Fit responses with samples on rows using the stored scaled SVD."""

        response_values = np.asarray(response, dtype=np.float64)
        if response_values.ndim != 2:
            raise DifferentialDesignDecompositionError(
                "differential response matrix must be two-dimensional"
            )
        if int(response_values.shape[0]) != self.sample_count:
            raise DifferentialDesignDecompositionError(
                "differential response matrix row count must match design samples; "
                f"response_rows={int(response_values.shape[0])}, "
                f"design_samples={self.sample_count}"
            )
        if not np.isfinite(response_values).all():
            raise DifferentialDesignDecompositionError(
                "differential response matrix must contain only finite numeric values"
            )

        singular_values = np.asarray(self.singular_values, dtype=np.float64)
        u_t_response = cast(
            npt.NDArray[np.float64],
            np.asarray(self._u.T @ response_values, dtype=np.float64),
        )
        scaled_coefficients = cast(
            npt.NDArray[np.float64],
            np.asarray(
                self._vt.T @ (u_t_response / singular_values[:, np.newaxis]),
                dtype=np.float64,
            ),
        )
        coefficients = scaled_coefficients / self._column_scales[:, np.newaxis]
        fitted_values = cast(
            npt.NDArray[np.float64],
            np.asarray(self._design_values @ coefficients, dtype=np.float64),
        )
        residuals = cast(
            npt.NDArray[np.float64],
            np.asarray(response_values - fitted_values, dtype=np.float64),
        )
        residuals = _zero_residual_roundoff(
            residuals=residuals,
            response_values=response_values,
            fitted_values=fitted_values,
            sample_count=self.sample_count,
            coefficient_count=self.coefficient_count,
        )
        rss = cast(
            npt.NDArray[np.float64],
            np.asarray(np.sum(np.square(residuals), axis=0), dtype=np.float64),
        )
        residual_variance = cast(
            npt.NDArray[np.float64],
            np.asarray(rss / self.residual_degrees_of_freedom, dtype=np.float64),
        )
        return DifferentialLinearFit(
            coefficients=coefficients,
            fitted_values=fitted_values,
            residuals=residuals,
            residual_sum_of_squares=rss,
            residual_variance=residual_variance,
        )

    def contrast_covariance(
        self,
        contrasts: npt.ArrayLike,
    ) -> npt.NDArray[np.float64]:
        """Return contrast covariance factors aligned to contrast columns."""

        contrast_values = _contrast_values(
            contrasts,
            coefficient_count=self.coefficient_count,
        )
        return cast(
            npt.NDArray[np.float64],
            np.asarray(
                contrast_values.T @ self._coefficient_covariance @ contrast_values,
                dtype=np.float64,
            ),
        )

    def contrast_scales(self, contrasts: npt.ArrayLike) -> npt.NDArray[np.float64]:
        """Return standard-error scale factors for each contrast column."""

        covariance = self.contrast_covariance(contrasts)
        diagonal = cast(
            npt.NDArray[np.float64],
            np.asarray(np.diagonal(covariance), dtype=np.float64),
        )
        invalid_negative = diagonal < -_NEGATIVE_VARIANCE_ROUNDOFF_TOLERANCE
        if np.any(invalid_negative):
            scale = np.full(diagonal.shape, np.nan, dtype=np.float64)
            scale[~invalid_negative] = np.sqrt(
                np.clip(diagonal[~invalid_negative], 0.0, None)
            )
            return scale
        return cast(
            npt.NDArray[np.float64],
            np.asarray(np.sqrt(np.clip(diagonal, 0.0, None)), dtype=np.float64),
        )

    def invalid_contrast_positions(self, contrasts: npt.ArrayLike) -> tuple[int, ...]:
        """Return contrast column positions that are zero-length or unstable."""

        scales = self.contrast_scales(contrasts)
        invalid_positions = np.flatnonzero(~np.isfinite(scales) | (scales <= 0.0))
        return tuple(int(position) for position in invalid_positions.tolist())

    def diagnostics_payload(self) -> dict[str, object]:
        """Return JSON-compatible numerical design diagnostics."""

        return {
            "decomposition_method": self.decomposition_method,
            "solver": self.solver,
            "column_scale_method": self.column_scale_method,
            "rank_tolerance_policy": self.rank_tolerance_policy,
            "rank_tolerance": float(self.rank_tolerance),
            "condition_number": float(self.condition_number),
            "max_condition_number": float(self.max_condition_number),
            "singular_values": [float(value) for value in self.singular_values],
        }


def decompose_differential_design(
    design: npt.ArrayLike,
    *,
    max_condition_number: float = DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER,
) -> DifferentialDesignDecomposition:
    """Build an admissible scaled-SVD decomposition for differential fitting."""

    design_values = np.asarray(design, dtype=np.float64)
    if design_values.ndim != 2:
        raise DifferentialDesignDecompositionError(
            "differential design matrix must be two-dimensional"
        )
    if not np.isfinite(design_values).all():
        raise DifferentialDesignDecompositionError(
            "differential design matrix must contain only finite numeric values"
        )

    sample_count = int(design_values.shape[0])
    coefficient_count = int(design_values.shape[1])
    if sample_count < 1 or coefficient_count < 1:
        raise DifferentialDesignDecompositionError(
            "differential design matrix must contain at least one sample and one "
            "coefficient"
        )
    max_condition = float(max_condition_number)
    if not math.isfinite(max_condition) or max_condition <= 1.0:
        raise DifferentialDesignDecompositionError(
            "differential design max_condition_number must be finite and > 1.0"
        )

    column_scales = cast(
        npt.NDArray[np.float64],
        np.asarray(np.linalg.norm(design_values, axis=0), dtype=np.float64),
    )
    invalid_scale_positions = np.flatnonzero(
        ~np.isfinite(column_scales) | (column_scales <= 0.0)
    )
    if invalid_scale_positions.size:
        raise DifferentialDesignDecompositionError(
            "differential design matrix contains zero-length or invalid coefficient "
            "columns; invalid_column_positions="
            + ", ".join(str(int(position)) for position in invalid_scale_positions)
        )

    scaled_design = design_values / column_scales[np.newaxis, :]
    try:
        u, singular_values, vt = np.linalg.svd(scaled_design, full_matrices=False)
    except np.linalg.LinAlgError as exc:
        raise DifferentialDesignDecompositionError(
            "scaled SVD failed for differential design matrix"
        ) from exc

    singular_values = cast(
        npt.NDArray[np.float64],
        np.asarray(singular_values, dtype=np.float64),
    )
    if singular_values.size == 0 or not np.isfinite(singular_values).all():
        raise DifferentialDesignDecompositionError(
            "scaled SVD produced invalid singular values for differential design matrix"
        )
    largest_singular_value = float(singular_values[0])
    rank_tolerance = (
        np.finfo(np.float64).eps
        * float(max(sample_count, coefficient_count))
        * largest_singular_value
    )
    rank = int(np.count_nonzero(singular_values > rank_tolerance))
    if rank < coefficient_count:
        raise DifferentialDesignDecompositionError(
            "differential design matrix is rank deficient under scaled SVD "
            "decomposition; "
            f"rank={rank}, columns={coefficient_count}, "
            f"rank_tolerance={rank_tolerance:.6g}, "
            f"singular_values={_format_singular_values(singular_values)}"
        )

    smallest_singular_value = float(singular_values[-1])
    condition_number = largest_singular_value / smallest_singular_value
    if not math.isfinite(condition_number) or condition_number > max_condition:
        raise DifferentialDesignDecompositionError(
            "differential design matrix is too ill-conditioned for stable "
            "moderated-contrast analysis under scaled SVD decomposition; "
            f"condition_number={condition_number:.6g}, "
            f"max_condition_number={max_condition:.6g}, "
            f"singular_values={_format_singular_values(singular_values)}"
        )

    residual_dof = float(sample_count - rank)
    if residual_dof <= 0.0:
        raise DifferentialDesignDecompositionError(
            "differential analysis residual degrees of freedom must be positive; "
            f"samples={sample_count}, rank={rank}, residual_dof={residual_dof}"
        )

    inverse_singular_squared = 1.0 / np.square(singular_values)
    scaled_covariance = cast(
        npt.NDArray[np.float64],
        np.asarray((vt.T * inverse_singular_squared[np.newaxis, :]) @ vt),
    )
    coefficient_covariance = scaled_covariance / (
        column_scales[:, np.newaxis] * column_scales[np.newaxis, :]
    )
    return DifferentialDesignDecomposition(
        sample_count=sample_count,
        coefficient_count=coefficient_count,
        rank=rank,
        residual_degrees_of_freedom=residual_dof,
        singular_values=tuple(float(value) for value in singular_values.tolist()),
        rank_tolerance=float(rank_tolerance),
        condition_number=float(condition_number),
        max_condition_number=max_condition,
        decomposition_method=DIFFERENTIAL_LINEAR_MODEL_DECOMPOSITION_METHOD,
        solver=DIFFERENTIAL_LINEAR_MODEL_SOLVER,
        column_scale_method=DIFFERENTIAL_LINEAR_MODEL_COLUMN_SCALE_METHOD,
        rank_tolerance_policy=DIFFERENTIAL_LINEAR_MODEL_RANK_TOLERANCE_POLICY,
        _design_values=_readonly_float_array(design_values),
        _column_scales=_readonly_float_array(column_scales),
        _u=_readonly_float_array(u),
        _vt=_readonly_float_array(vt),
        _coefficient_covariance=_readonly_float_array(coefficient_covariance),
    )


def _contrast_values(
    contrasts: npt.ArrayLike,
    *,
    coefficient_count: int,
) -> npt.NDArray[np.float64]:
    contrast_values = np.asarray(contrasts, dtype=np.float64)
    if contrast_values.ndim != 2:
        raise DifferentialDesignDecompositionError(
            "differential contrast matrix must be two-dimensional"
        )
    if int(contrast_values.shape[0]) != coefficient_count:
        raise DifferentialDesignDecompositionError(
            "differential contrast matrix row count must match design coefficients; "
            f"contrast_rows={int(contrast_values.shape[0])}, "
            f"design_coefficients={coefficient_count}"
        )
    if not np.isfinite(contrast_values).all():
        raise DifferentialDesignDecompositionError(
            "differential contrast matrix must contain only finite numeric values"
        )
    return cast(npt.NDArray[np.float64], contrast_values)


def _zero_residual_roundoff(
    *,
    residuals: npt.NDArray[np.float64],
    response_values: npt.NDArray[np.float64],
    fitted_values: npt.NDArray[np.float64],
    sample_count: int,
    coefficient_count: int,
) -> npt.NDArray[np.float64]:
    response_scale = np.nanmax(np.abs(response_values), axis=0)
    fitted_scale = np.nanmax(np.abs(fitted_values), axis=0)
    scale = cast(
        npt.NDArray[np.float64],
        np.maximum(np.maximum(response_scale, fitted_scale), 1.0),
    )
    tolerance = (
        np.finfo(np.float64).eps * float(max(sample_count, coefficient_count)) * scale
    )
    zero_mask = np.isfinite(residuals) & np.isfinite(tolerance[np.newaxis, :])
    zero_mask &= np.abs(residuals) <= tolerance[np.newaxis, :]
    return cast(
        npt.NDArray[np.float64],
        np.where(zero_mask, 0.0, residuals),
    )


def _format_singular_values(values: npt.NDArray[np.float64]) -> str:
    preview_values = [f"{float(value):.6g}" for value in values[:5]]
    suffix = "" if int(values.size) <= 5 else ", ..."
    return "[" + ", ".join(preview_values) + suffix + "]"


def _readonly_float_array(values: npt.ArrayLike) -> npt.NDArray[np.float64]:
    array = np.array(values, dtype=np.float64, copy=True)
    array.setflags(write=False)
    return cast(npt.NDArray[np.float64], array)


__all__ = [
    "DIFFERENTIAL_LINEAR_MODEL_COLUMN_SCALE_METHOD",
    "DIFFERENTIAL_LINEAR_MODEL_DECOMPOSITION_METHOD",
    "DIFFERENTIAL_LINEAR_MODEL_MAX_CONDITION_NUMBER",
    "DIFFERENTIAL_LINEAR_MODEL_RANK_TOLERANCE_POLICY",
    "DIFFERENTIAL_LINEAR_MODEL_SOLVER",
    "DifferentialDesignDecomposition",
    "DifferentialDesignDecompositionError",
    "DifferentialLinearFit",
    "decompose_differential_design",
]
