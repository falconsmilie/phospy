"""Diagnostic payload models for differential analysis."""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.comparison import series_equals
from phospy.frames.ownership import export_series, own_series
from phospy.science.differential.models.provenance import (
    DifferentialContrastDefinition,
)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class EmpiricalBayesPriorDiagnostics:
    """Diagnostics for prior-variance and prior-df estimation.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit diagnostics-content comparison.
    """

    __hash__ = object.__hash__

    method: str
    robust: bool
    trend: bool
    winsor_tail_p: tuple[float, float]
    base_prior_variance: float
    base_prior_degrees_of_freedom: float
    robust_outlier_count: int
    robust_outlier_fraction: float
    winsorized_low_count: int
    winsorized_high_count: int
    prior_variance: pd.Series
    prior_degrees_of_freedom: pd.Series

    def __init__(
        self,
        *,
        method: str,
        robust: bool,
        trend: bool,
        winsor_tail_p: tuple[float, float],
        base_prior_variance: float,
        base_prior_degrees_of_freedom: float,
        robust_outlier_count: int,
        robust_outlier_fraction: float,
        winsorized_low_count: int,
        winsorized_high_count: int,
        prior_variance: pd.Series,
        prior_degrees_of_freedom: pd.Series,
        _assume_owned: bool = False,
    ) -> None:
        prior_variance = own_series(
            prior_variance,
            field_name="differential_result.prior_diagnostics.prior_variance",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        prior_degrees_of_freedom = own_series(
            prior_degrees_of_freedom,
            field_name=(
                "differential_result.prior_diagnostics.prior_degrees_of_freedom"
            ),
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        if not prior_variance.index.equals(prior_degrees_of_freedom.index):
            raise PhosPyInputError(
                "differential_result prior diagnostics index mismatch between "
                "prior_variance and prior_degrees_of_freedom"
            )
        object.__setattr__(self, "method", str(method))
        object.__setattr__(self, "robust", bool(robust))
        object.__setattr__(self, "trend", bool(trend))
        object.__setattr__(
            self,
            "winsor_tail_p",
            (float(winsor_tail_p[0]), float(winsor_tail_p[1])),
        )
        object.__setattr__(self, "base_prior_variance", float(base_prior_variance))
        object.__setattr__(
            self,
            "base_prior_degrees_of_freedom",
            float(base_prior_degrees_of_freedom),
        )
        object.__setattr__(self, "robust_outlier_count", int(robust_outlier_count))
        object.__setattr__(
            self,
            "robust_outlier_fraction",
            float(robust_outlier_fraction),
        )
        object.__setattr__(self, "winsorized_low_count", int(winsorized_low_count))
        object.__setattr__(self, "winsorized_high_count", int(winsorized_high_count))
        object.__setattr__(self, "prior_variance", prior_variance)
        object.__setattr__(self, "prior_degrees_of_freedom", prior_degrees_of_freedom)

    def prior_variance_series(self) -> pd.Series:
        return export_series(self.prior_variance)

    def prior_degrees_of_freedom_series(self) -> pd.Series:
        return export_series(self.prior_degrees_of_freedom)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another prior diagnostics object has same content."""

        if not isinstance(other, EmpiricalBayesPriorDiagnostics):
            return False
        return (
            self.method == other.method
            and self.robust == other.robust
            and self.trend == other.trend
            and self.winsor_tail_p == other.winsor_tail_p
            and self.base_prior_variance == other.base_prior_variance
            and self.base_prior_degrees_of_freedom
            == other.base_prior_degrees_of_freedom
            and self.robust_outlier_count == other.robust_outlier_count
            and self.robust_outlier_fraction == other.robust_outlier_fraction
            and self.winsorized_low_count == other.winsorized_low_count
            and self.winsorized_high_count == other.winsorized_high_count
            and series_equals(self.prior_variance, other.prior_variance)
            and series_equals(
                self.prior_degrees_of_freedom,
                other.prior_degrees_of_freedom,
            )
        )


@dataclass(frozen=True, slots=True, init=False, eq=False)
class MeanVarianceTrendDiagnostics:
    """Diagnostics payload for mean-intensity vs variance trend fitting.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit diagnostics-content comparison.
    """

    __hash__ = object.__hash__

    mean_intensity: pd.Series
    log_residual_variance: pd.Series
    fitted_log_prior_variance: pd.Series

    def __init__(
        self,
        *,
        mean_intensity: pd.Series,
        log_residual_variance: pd.Series,
        fitted_log_prior_variance: pd.Series,
        _assume_owned: bool = False,
    ) -> None:
        mean_intensity = own_series(
            mean_intensity,
            field_name="differential_result.mean_variance_trend.mean_intensity",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        log_residual_variance = own_series(
            log_residual_variance,
            field_name=(
                "differential_result.mean_variance_trend.log_residual_variance"
            ),
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        fitted_log_prior_variance = own_series(
            fitted_log_prior_variance,
            field_name=(
                "differential_result.mean_variance_trend.fitted_log_prior_variance"
            ),
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        if not mean_intensity.index.equals(log_residual_variance.index):
            raise PhosPyInputError(
                "mean-variance trend diagnostics index mismatch for mean_intensity and "
                "log_residual_variance"
            )
        if not mean_intensity.index.equals(fitted_log_prior_variance.index):
            raise PhosPyInputError(
                "mean-variance trend diagnostics index mismatch for mean_intensity and "
                "fitted_log_prior_variance"
            )
        object.__setattr__(self, "mean_intensity", mean_intensity)
        object.__setattr__(self, "log_residual_variance", log_residual_variance)
        object.__setattr__(
            self,
            "fitted_log_prior_variance",
            fitted_log_prior_variance,
        )

    def mean_intensity_series(self) -> pd.Series:
        return export_series(self.mean_intensity)

    def log_residual_variance_series(self) -> pd.Series:
        return export_series(self.log_residual_variance)

    def fitted_log_prior_variance_series(self) -> pd.Series:
        return export_series(self.fitted_log_prior_variance)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another trend diagnostics object has same content."""

        if not isinstance(other, MeanVarianceTrendDiagnostics):
            return False
        return (
            series_equals(self.mean_intensity, other.mean_intensity)
            and series_equals(
                self.log_residual_variance,
                other.log_residual_variance,
            )
            and series_equals(
                self.fitted_log_prior_variance,
                other.fitted_log_prior_variance,
            )
        )


@dataclass(frozen=True, slots=True)
class DifferentialModelDiagnostics:
    """User-visible scope and model diagnostics for differential results."""

    model_type: str
    design_columns: tuple[str, ...]
    contrast_definitions: tuple[DifferentialContrastDefinition, ...]
    rank: int
    n_samples: int
    n_sites: int
    residual_degrees_of_freedom: float
    variance_method: str
    moderation_method: str
    multiple_testing_method: str
    imputation_policy: str
    missing_value_policy: str
    intensity_scale: str
    normalisation_state: str
    batch_or_covariate_terms: tuple[str, ...]
    unsupported_assumptions: tuple[str, ...]
    warnings: tuple[str, ...]
    decomposition_method: str = "not_recorded"
    solver: str = "not_recorded"
    column_scale_method: str = "not_recorded"
    rank_tolerance_policy: str = "not_recorded"
    rank_tolerance: float = 0.0
    condition_number: float = 0.0
    max_condition_number: float = 0.0
    singular_values: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        rank = _require_non_negative_int(
            self.rank,
            field_name="differential_result.diagnostics.rank",
        )
        n_samples = _require_non_negative_int(
            self.n_samples,
            field_name="differential_result.diagnostics.n_samples",
        )
        n_sites = _require_non_negative_int(
            self.n_sites,
            field_name="differential_result.diagnostics.n_sites",
        )
        residual_dof = _require_finite_float(
            self.residual_degrees_of_freedom,
            field_name=("differential_result.diagnostics.residual_degrees_of_freedom"),
        )
        if residual_dof < 0.0:
            raise PhosPyInputError(
                "differential_result.diagnostics.residual_degrees_of_freedom "
                "must be >= 0.0"
            )
        rank_tolerance = _require_non_negative_finite_float(
            self.rank_tolerance,
            field_name="differential_result.diagnostics.rank_tolerance",
        )
        condition_number = _require_non_negative_finite_float(
            self.condition_number,
            field_name="differential_result.diagnostics.condition_number",
        )
        max_condition_number = _require_non_negative_finite_float(
            self.max_condition_number,
            field_name="differential_result.diagnostics.max_condition_number",
        )
        singular_values = tuple(float(value) for value in self.singular_values)
        if any(value < 0.0 or not math.isfinite(value) for value in singular_values):
            raise PhosPyInputError(
                "differential_result.diagnostics.singular_values must contain "
                "finite values >= 0.0"
            )
        contrast_definitions = tuple(self.contrast_definitions)
        for definition in contrast_definitions:
            if not isinstance(
                cast(object, definition),
                DifferentialContrastDefinition,
            ):
                raise PhosPyInputError(
                    "differential_result.diagnostics.contrast_definitions must "
                    "contain DifferentialContrastDefinition values"
                )
        object.__setattr__(
            self,
            "model_type",
            _require_non_empty_text(
                self.model_type,
                field_name="differential_result.diagnostics.model_type",
            ),
        )
        object.__setattr__(
            self,
            "design_columns",
            _text_tuple(
                self.design_columns,
                field_name="differential_result.diagnostics.design_columns",
            ),
        )
        object.__setattr__(self, "contrast_definitions", contrast_definitions)
        object.__setattr__(self, "rank", rank)
        object.__setattr__(self, "n_samples", n_samples)
        object.__setattr__(self, "n_sites", n_sites)
        object.__setattr__(self, "residual_degrees_of_freedom", residual_dof)
        object.__setattr__(
            self,
            "decomposition_method",
            _require_non_empty_text(
                self.decomposition_method,
                field_name="differential_result.diagnostics.decomposition_method",
            ),
        )
        object.__setattr__(
            self,
            "solver",
            _require_non_empty_text(
                self.solver,
                field_name="differential_result.diagnostics.solver",
            ),
        )
        object.__setattr__(
            self,
            "column_scale_method",
            _require_non_empty_text(
                self.column_scale_method,
                field_name="differential_result.diagnostics.column_scale_method",
            ),
        )
        object.__setattr__(
            self,
            "rank_tolerance_policy",
            _require_non_empty_text(
                self.rank_tolerance_policy,
                field_name="differential_result.diagnostics.rank_tolerance_policy",
            ),
        )
        object.__setattr__(self, "rank_tolerance", rank_tolerance)
        object.__setattr__(self, "condition_number", condition_number)
        object.__setattr__(self, "max_condition_number", max_condition_number)
        object.__setattr__(self, "singular_values", singular_values)
        object.__setattr__(
            self,
            "variance_method",
            _require_non_empty_text(
                self.variance_method,
                field_name="differential_result.diagnostics.variance_method",
            ),
        )
        object.__setattr__(
            self,
            "moderation_method",
            _require_non_empty_text(
                self.moderation_method,
                field_name="differential_result.diagnostics.moderation_method",
            ),
        )
        object.__setattr__(
            self,
            "multiple_testing_method",
            _require_non_empty_text(
                self.multiple_testing_method,
                field_name="differential_result.diagnostics.multiple_testing_method",
            ),
        )
        object.__setattr__(
            self,
            "imputation_policy",
            _require_non_empty_text(
                self.imputation_policy,
                field_name="differential_result.diagnostics.imputation_policy",
            ),
        )
        object.__setattr__(
            self,
            "missing_value_policy",
            _require_non_empty_text(
                self.missing_value_policy,
                field_name="differential_result.diagnostics.missing_value_policy",
            ),
        )
        object.__setattr__(
            self,
            "intensity_scale",
            _require_non_empty_text(
                self.intensity_scale,
                field_name="differential_result.diagnostics.intensity_scale",
            ),
        )
        object.__setattr__(
            self,
            "normalisation_state",
            _require_non_empty_text(
                self.normalisation_state,
                field_name="differential_result.diagnostics.normalisation_state",
            ),
        )
        object.__setattr__(
            self,
            "batch_or_covariate_terms",
            _text_tuple(
                self.batch_or_covariate_terms,
                field_name=("differential_result.diagnostics.batch_or_covariate_terms"),
            ),
        )
        object.__setattr__(
            self,
            "unsupported_assumptions",
            _text_tuple(
                self.unsupported_assumptions,
                field_name=("differential_result.diagnostics.unsupported_assumptions"),
            ),
        )
        object.__setattr__(
            self,
            "warnings",
            _text_tuple(
                self.warnings,
                field_name="differential_result.diagnostics.warnings",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible diagnostics payload."""

        return {
            "model_type": self.model_type,
            "design_columns": list(self.design_columns),
            "contrast_definitions": [
                _contrast_definition_payload(definition)
                for definition in self.contrast_definitions
            ],
            "rank": self.rank,
            "n_samples": self.n_samples,
            "n_sites": self.n_sites,
            "residual_degrees_of_freedom": self.residual_degrees_of_freedom,
            "decomposition_method": self.decomposition_method,
            "solver": self.solver,
            "column_scale_method": self.column_scale_method,
            "rank_tolerance_policy": self.rank_tolerance_policy,
            "rank_tolerance": self.rank_tolerance,
            "condition_number": self.condition_number,
            "max_condition_number": self.max_condition_number,
            "singular_values": list(self.singular_values),
            "variance_method": self.variance_method,
            "moderation_method": self.moderation_method,
            "multiple_testing_method": self.multiple_testing_method,
            "imputation_policy": self.imputation_policy,
            "missing_value_policy": self.missing_value_policy,
            "intensity_scale": self.intensity_scale,
            "normalisation_state": self.normalisation_state,
            "batch_or_covariate_terms": list(self.batch_or_covariate_terms),
            "unsupported_assumptions": list(self.unsupported_assumptions),
            "warnings": list(self.warnings),
        }


def _contrast_definition_payload(
    definition: DifferentialContrastDefinition,
) -> dict[str, object]:
    return {
        "name": definition.name,
        "numerator_condition": definition.numerator_condition,
        "denominator_condition": definition.denominator_condition,
        "coefficients": [
            {"coefficient": coefficient, "weight": weight}
            for coefficient, weight in definition.coefficients
        ],
        "description": definition.description,
    }


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return value.strip()


def _text_tuple(value: object, *, field_name: str) -> tuple[str, ...]:
    if isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    if not isinstance(value, Iterable):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    values = tuple(cast(Iterable[object], value))
    normalized: list[str] = []
    for item in values:
        normalized.append(
            _require_non_empty_text(
                item,
                field_name=f"{field_name}[]",
            )
        )
    return tuple(normalized)


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be a non-negative integer")
    if value < 0:
        raise PhosPyInputError(f"{field_name} must be >= 0")
    return int(value)


def _require_finite_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be a finite numeric value")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise PhosPyInputError(f"{field_name} must be finite")
    return numeric


def _require_non_negative_finite_float(value: object, *, field_name: str) -> float:
    numeric = _require_finite_float(value, field_name=field_name)
    if numeric < 0.0:
        raise PhosPyInputError(f"{field_name} must be >= 0.0")
    return numeric


__all__ = [
    "DifferentialModelDiagnostics",
    "EmpiricalBayesPriorDiagnostics",
    "MeanVarianceTrendDiagnostics",
]
