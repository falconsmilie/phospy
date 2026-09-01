"""Diagnostic payload models for differential analysis."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.comparison import dataframe_equals, series_equals
from phospy.frames.ownership import (
    export_dataframe,
    export_series,
    own_dataframe,
    own_series,
)
from phospy.frames.validation import (
    require_columns,
    require_dataframe,
    require_non_empty_dataframe,
    require_string_index,
    require_unique_columns,
    require_unique_index,
)
from phospy.science.configs.differential import (
    DIFFERENTIAL_PROTEIN_AWARE_METHOD_PROTEIN_COVARIATE_ADJUSTED_MODERATED_LINEAR_MODEL_V1 as _PROTEIN_AWARE_SUPPORTED_METHOD_ID,
)
from phospy.science.differential.models.provenance import (
    DifferentialContrastDefinition,
)
from phospy.science.differential.models.tables import (
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_TESTED,
    DIFFERENTIAL_RESULT_WITHHELD_STATUSES,
    validate_result_status_reason_contract,
)

PROTEIN_AWARE_DIFFERENTIAL_CLAIM_STATUS_EXPERIMENTAL = "experimental"
PROTEIN_AWARE_DIFFERENTIAL_MODEL_TYPE = "protein_covariate_adjusted_moderated_ols"
PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_SCOPE = (
    "successfully_fitted_augmented_designs_by_total_protein_row"
)
PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_STATISTIC = (
    "min_median_max_condition_number_across_fitted_augmented_designs"
)
PROTEIN_AWARE_DIFFERENTIAL_CENTERING_POLICY = "mean_centered_no_standardization"
PROTEIN_AWARE_DIFFERENTIAL_IMPUTATION_POLICY = "none"
PROTEIN_AWARE_DIFFERENTIAL_FALLBACK_POLICY = "no_fallback_to_ordinary_differential_lane"
PROTEIN_AWARE_DIFFERENTIAL_PER_SITE_DIAGNOSTIC_COLUMNS: tuple[str, ...] = (
    "site_key",
    "protein_identifier",
    "total_protein_row_key",
    "protein_aware_preparation_eligibility",
    "protein_aware_preparation_reasons",
    "protein_covariate_raw_mean",
    "protein_covariate_raw_standard_deviation",
    "protein_covariate_centered_standard_deviation",
    "protein_augmented_design_rank",
    "protein_augmented_design_residual_degrees_of_freedom",
    "protein_augmented_design_condition_number",
    "protein_contrast_estimability_status",
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
    "protein_covariate_coefficient",
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


@dataclass(frozen=True, slots=True, init=False, eq=False)
class ProteinAwareDifferentialDiagnostics:
    """Protein-aware differential diagnostics with per-site model facts.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit diagnostics-content comparison.
    """

    __hash__ = object.__hash__

    method_id: str
    claim_status: str
    preparation_policy: str
    protein_mapping_policy: str
    protein_covariate_centered: bool
    protein_covariate_standardized: bool
    protein_covariate_centering_policy: str
    protein_covariate_imputation_policy: str
    fallback_policy: str
    total_site_count: int
    ordinary_eligible_site_count: int
    protein_preparation_eligible_site_count: int
    tested_site_count: int
    withheld_site_count: int
    status_counts: tuple[tuple[str, int], ...]
    reason_counts: tuple[tuple[str, int], ...]
    execution_sample_order: tuple[str, ...]
    base_design_rank: int
    base_residual_degrees_of_freedom: float
    expected_augmented_rank: int
    common_augmented_rank: int
    common_augmented_residual_degrees_of_freedom: float
    distinct_matched_protein_row_count: int
    fitted_protein_row_count: int
    condition_number_summary_scope: str
    condition_number_summary_statistic: str
    min_augmented_condition_number: float | None
    median_augmented_condition_number: float | None
    max_augmented_condition_number: float | None
    _per_site_diagnostics: pd.DataFrame

    def __init__(
        self,
        *,
        method_id: str,
        claim_status: str,
        preparation_policy: str,
        protein_mapping_policy: str,
        protein_covariate_centered: bool,
        protein_covariate_standardized: bool,
        protein_covariate_centering_policy: str,
        protein_covariate_imputation_policy: str,
        fallback_policy: str,
        total_site_count: int,
        ordinary_eligible_site_count: int,
        protein_preparation_eligible_site_count: int,
        tested_site_count: int,
        withheld_site_count: int,
        status_counts: Iterable[tuple[str, int]],
        reason_counts: Iterable[tuple[str, int]],
        execution_sample_order: Iterable[str],
        base_design_rank: int,
        base_residual_degrees_of_freedom: float,
        expected_augmented_rank: int,
        common_augmented_rank: int,
        common_augmented_residual_degrees_of_freedom: float,
        distinct_matched_protein_row_count: int,
        fitted_protein_row_count: int,
        condition_number_summary_scope: str,
        condition_number_summary_statistic: str,
        min_augmented_condition_number: float | None,
        median_augmented_condition_number: float | None,
        max_augmented_condition_number: float | None,
        per_site_diagnostics: pd.DataFrame,
        _assume_owned: bool = False,
    ) -> None:
        total_site_count = _require_non_negative_int(
            total_site_count,
            field_name="protein_aware_diagnostics.total_site_count",
        )
        ordinary_eligible_site_count = _require_non_negative_int(
            ordinary_eligible_site_count,
            field_name="protein_aware_diagnostics.ordinary_eligible_site_count",
        )
        protein_preparation_eligible_site_count = _require_non_negative_int(
            protein_preparation_eligible_site_count,
            field_name=(
                "protein_aware_diagnostics.protein_preparation_eligible_site_count"
            ),
        )
        tested_site_count = _require_non_negative_int(
            tested_site_count,
            field_name="protein_aware_diagnostics.tested_site_count",
        )
        withheld_site_count = _require_non_negative_int(
            withheld_site_count,
            field_name="protein_aware_diagnostics.withheld_site_count",
        )
        if tested_site_count + withheld_site_count != total_site_count:
            raise PhosPyInputError(
                "protein_aware_diagnostics total_site_count must equal "
                "tested_site_count + withheld_site_count"
            )
        if ordinary_eligible_site_count > total_site_count:
            raise PhosPyInputError(
                "protein_aware_diagnostics.ordinary_eligible_site_count cannot exceed "
                "total_site_count"
            )
        if protein_preparation_eligible_site_count > ordinary_eligible_site_count:
            raise PhosPyInputError(
                "protein_aware_diagnostics.protein_preparation_eligible_site_count "
                "cannot exceed ordinary_eligible_site_count"
            )
        if tested_site_count > protein_preparation_eligible_site_count:
            raise PhosPyInputError(
                "protein_aware_diagnostics.tested_site_count cannot exceed "
                "protein_preparation_eligible_site_count"
            )

        base_residual_dof = _require_finite_float(
            base_residual_degrees_of_freedom,
            field_name=("protein_aware_diagnostics.base_residual_degrees_of_freedom"),
        )
        common_augmented_dof = _require_finite_float(
            common_augmented_residual_degrees_of_freedom,
            field_name=(
                "protein_aware_diagnostics.common_augmented_residual_degrees_of_freedom"
            ),
        )
        if base_residual_dof < 0.0 or common_augmented_dof < 0.0:
            raise PhosPyInputError(
                "protein_aware_diagnostics residual degrees of freedom values must "
                "be >= 0.0"
            )

        matched_protein_rows = _require_non_negative_int(
            distinct_matched_protein_row_count,
            field_name=("protein_aware_diagnostics.distinct_matched_protein_row_count"),
        )
        fitted_protein_rows = _require_non_negative_int(
            fitted_protein_row_count,
            field_name="protein_aware_diagnostics.fitted_protein_row_count",
        )
        if fitted_protein_rows > matched_protein_rows:
            raise PhosPyInputError(
                "protein_aware_diagnostics.fitted_protein_row_count cannot exceed "
                "distinct_matched_protein_row_count"
            )

        owned_per_site = own_dataframe(
            per_site_diagnostics,
            field_name="protein_aware_diagnostics.per_site_diagnostics",
            error_type=PhosPyInputError,
            assume_owned=_assume_owned,
        )
        _validate_protein_aware_per_site_diagnostics(
            owned_per_site,
            expected_total_site_count=total_site_count,
            expected_tested_site_count=tested_site_count,
            expected_withheld_site_count=withheld_site_count,
        )

        normalized_status_counts = _normalize_count_pairs(
            status_counts,
            field_name="protein_aware_diagnostics.status_counts",
        )
        normalized_reason_counts = _normalize_count_pairs(
            reason_counts,
            field_name="protein_aware_diagnostics.reason_counts",
        )
        _require_counts_match_per_site_diagnostics(
            per_site_diagnostics=owned_per_site,
            status_counts=normalized_status_counts,
            reason_counts=normalized_reason_counts,
        )

        object.__setattr__(
            self,
            "method_id",
            _require_exact_text(
                method_id,
                field_name="protein_aware_diagnostics.method_id",
                expected=_PROTEIN_AWARE_SUPPORTED_METHOD_ID,
            ),
        )
        object.__setattr__(
            self,
            "claim_status",
            _require_exact_text(
                claim_status,
                field_name="protein_aware_diagnostics.claim_status",
                expected=PROTEIN_AWARE_DIFFERENTIAL_CLAIM_STATUS_EXPERIMENTAL,
            ),
        )
        object.__setattr__(
            self,
            "preparation_policy",
            _require_non_empty_text(
                preparation_policy,
                field_name="protein_aware_diagnostics.preparation_policy",
            ),
        )
        object.__setattr__(
            self,
            "protein_mapping_policy",
            _require_non_empty_text(
                protein_mapping_policy,
                field_name="protein_aware_diagnostics.protein_mapping_policy",
            ),
        )
        object.__setattr__(
            self,
            "protein_covariate_centered",
            _require_exact_bool(
                protein_covariate_centered,
                field_name="protein_aware_diagnostics.protein_covariate_centered",
                expected=True,
            ),
        )
        object.__setattr__(
            self,
            "protein_covariate_standardized",
            _require_exact_bool(
                protein_covariate_standardized,
                field_name="protein_aware_diagnostics.protein_covariate_standardized",
                expected=False,
            ),
        )
        object.__setattr__(
            self,
            "protein_covariate_centering_policy",
            _require_exact_text(
                protein_covariate_centering_policy,
                field_name=(
                    "protein_aware_diagnostics.protein_covariate_centering_policy"
                ),
                expected=PROTEIN_AWARE_DIFFERENTIAL_CENTERING_POLICY,
            ),
        )
        object.__setattr__(
            self,
            "protein_covariate_imputation_policy",
            _require_exact_text(
                protein_covariate_imputation_policy,
                field_name=(
                    "protein_aware_diagnostics.protein_covariate_imputation_policy"
                ),
                expected=PROTEIN_AWARE_DIFFERENTIAL_IMPUTATION_POLICY,
            ),
        )
        object.__setattr__(
            self,
            "fallback_policy",
            _require_exact_text(
                fallback_policy,
                field_name="protein_aware_diagnostics.fallback_policy",
                expected=PROTEIN_AWARE_DIFFERENTIAL_FALLBACK_POLICY,
            ),
        )
        object.__setattr__(self, "total_site_count", total_site_count)
        object.__setattr__(
            self,
            "ordinary_eligible_site_count",
            ordinary_eligible_site_count,
        )
        object.__setattr__(
            self,
            "protein_preparation_eligible_site_count",
            protein_preparation_eligible_site_count,
        )
        object.__setattr__(self, "tested_site_count", tested_site_count)
        object.__setattr__(self, "withheld_site_count", withheld_site_count)
        object.__setattr__(self, "status_counts", normalized_status_counts)
        object.__setattr__(self, "reason_counts", normalized_reason_counts)
        object.__setattr__(
            self,
            "execution_sample_order",
            _text_tuple(
                execution_sample_order,
                field_name="protein_aware_diagnostics.execution_sample_order",
            ),
        )
        object.__setattr__(
            self,
            "base_design_rank",
            _require_non_negative_int(
                base_design_rank,
                field_name="protein_aware_diagnostics.base_design_rank",
            ),
        )
        object.__setattr__(
            self,
            "base_residual_degrees_of_freedom",
            base_residual_dof,
        )
        object.__setattr__(
            self,
            "expected_augmented_rank",
            _require_non_negative_int(
                expected_augmented_rank,
                field_name="protein_aware_diagnostics.expected_augmented_rank",
            ),
        )
        object.__setattr__(
            self,
            "common_augmented_rank",
            _require_non_negative_int(
                common_augmented_rank,
                field_name="protein_aware_diagnostics.common_augmented_rank",
            ),
        )
        object.__setattr__(
            self,
            "common_augmented_residual_degrees_of_freedom",
            common_augmented_dof,
        )
        object.__setattr__(
            self,
            "distinct_matched_protein_row_count",
            matched_protein_rows,
        )
        object.__setattr__(self, "fitted_protein_row_count", fitted_protein_rows)
        object.__setattr__(
            self,
            "condition_number_summary_scope",
            _require_exact_text(
                condition_number_summary_scope,
                field_name=("protein_aware_diagnostics.condition_number_summary_scope"),
                expected=PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_SCOPE,
            ),
        )
        object.__setattr__(
            self,
            "condition_number_summary_statistic",
            _require_exact_text(
                condition_number_summary_statistic,
                field_name=(
                    "protein_aware_diagnostics.condition_number_summary_statistic"
                ),
                expected=(
                    PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_STATISTIC
                ),
            ),
        )
        object.__setattr__(
            self,
            "min_augmented_condition_number",
            _optional_non_negative_finite_float(
                min_augmented_condition_number,
                field_name=("protein_aware_diagnostics.min_augmented_condition_number"),
            ),
        )
        object.__setattr__(
            self,
            "median_augmented_condition_number",
            _optional_non_negative_finite_float(
                median_augmented_condition_number,
                field_name=(
                    "protein_aware_diagnostics.median_augmented_condition_number"
                ),
            ),
        )
        object.__setattr__(
            self,
            "max_augmented_condition_number",
            _optional_non_negative_finite_float(
                max_augmented_condition_number,
                field_name=("protein_aware_diagnostics.max_augmented_condition_number"),
            ),
        )
        object.__setattr__(self, "_per_site_diagnostics", owned_per_site)

    @property
    def per_site_diagnostics(self) -> pd.DataFrame:
        return export_dataframe(self._per_site_diagnostics)

    def per_site_diagnostics_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self._per_site_diagnostics)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another diagnostics object has same content."""

        if not isinstance(other, ProteinAwareDifferentialDiagnostics):
            return False
        return (
            self.method_id == other.method_id
            and self.claim_status == other.claim_status
            and self.preparation_policy == other.preparation_policy
            and self.protein_mapping_policy == other.protein_mapping_policy
            and self.protein_covariate_centered == other.protein_covariate_centered
            and self.protein_covariate_standardized
            == other.protein_covariate_standardized
            and self.protein_covariate_centering_policy
            == other.protein_covariate_centering_policy
            and self.protein_covariate_imputation_policy
            == other.protein_covariate_imputation_policy
            and self.fallback_policy == other.fallback_policy
            and self.total_site_count == other.total_site_count
            and self.ordinary_eligible_site_count == other.ordinary_eligible_site_count
            and self.protein_preparation_eligible_site_count
            == other.protein_preparation_eligible_site_count
            and self.tested_site_count == other.tested_site_count
            and self.withheld_site_count == other.withheld_site_count
            and self.status_counts == other.status_counts
            and self.reason_counts == other.reason_counts
            and self.execution_sample_order == other.execution_sample_order
            and self.base_design_rank == other.base_design_rank
            and self.base_residual_degrees_of_freedom
            == other.base_residual_degrees_of_freedom
            and self.expected_augmented_rank == other.expected_augmented_rank
            and self.common_augmented_rank == other.common_augmented_rank
            and self.common_augmented_residual_degrees_of_freedom
            == other.common_augmented_residual_degrees_of_freedom
            and self.distinct_matched_protein_row_count
            == other.distinct_matched_protein_row_count
            and self.fitted_protein_row_count == other.fitted_protein_row_count
            and self.condition_number_summary_scope
            == other.condition_number_summary_scope
            and self.condition_number_summary_statistic
            == other.condition_number_summary_statistic
            and self.min_augmented_condition_number
            == other.min_augmented_condition_number
            and self.median_augmented_condition_number
            == other.median_augmented_condition_number
            and self.max_augmented_condition_number
            == other.max_augmented_condition_number
            and dataframe_equals(
                self._per_site_diagnostics,
                other._per_site_diagnostics,
            )
        )

    def to_payload(self) -> dict[str, object]:
        """Return a deterministic JSON-compatible diagnostics payload."""

        return {
            "method_id": self.method_id,
            "claim_status": self.claim_status,
            "preparation_policy": self.preparation_policy,
            "protein_mapping_policy": self.protein_mapping_policy,
            "protein_covariate": {
                "centered": self.protein_covariate_centered,
                "standardized": self.protein_covariate_standardized,
                "centering_policy": self.protein_covariate_centering_policy,
                "imputation_policy": self.protein_covariate_imputation_policy,
                "fallback_policy": self.fallback_policy,
            },
            "counts": {
                "total_site_count": self.total_site_count,
                "ordinary_eligible_site_count": self.ordinary_eligible_site_count,
                "protein_preparation_eligible_site_count": (
                    self.protein_preparation_eligible_site_count
                ),
                "tested_site_count": self.tested_site_count,
                "withheld_site_count": self.withheld_site_count,
                "distinct_matched_protein_row_count": (
                    self.distinct_matched_protein_row_count
                ),
                "fitted_protein_row_count": self.fitted_protein_row_count,
            },
            "status_counts": [
                {"status": status, "count": count}
                for status, count in self.status_counts
            ],
            "reason_counts": [
                {"reason": reason, "count": count}
                for reason, count in self.reason_counts
            ],
            "execution_sample_order": list(self.execution_sample_order),
            "design_diagnostics": {
                "base_design_rank": self.base_design_rank,
                "base_residual_degrees_of_freedom": (
                    self.base_residual_degrees_of_freedom
                ),
                "expected_augmented_rank": self.expected_augmented_rank,
                "common_augmented_rank": self.common_augmented_rank,
                "common_augmented_residual_degrees_of_freedom": (
                    self.common_augmented_residual_degrees_of_freedom
                ),
            },
            "condition_number_summary": {
                "scope": self.condition_number_summary_scope,
                "statistic": self.condition_number_summary_statistic,
                "min": self.min_augmented_condition_number,
                "median": self.median_augmented_condition_number,
                "max": self.max_augmented_condition_number,
            },
            "per_site_diagnostics": _dataframe_records_payload(
                self._per_site_diagnostics
            ),
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


def _require_exact_text(value: object, *, field_name: str, expected: str) -> str:
    text = _require_non_empty_text(value, field_name=field_name)
    if text != expected:
        raise PhosPyInputError(f"{field_name} must be {expected!r}")
    return text


def _require_exact_bool(value: object, *, field_name: str, expected: bool) -> bool:
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a boolean")
    if value is not expected:
        raise PhosPyInputError(f"{field_name} must be {expected!r}")
    return value


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


def _optional_non_negative_finite_float(
    value: object | None,
    *,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    return _require_non_negative_finite_float(value, field_name=field_name)


def _normalize_count_pairs(
    value: object,
    *,
    field_name: str,
) -> tuple[tuple[str, int], ...]:
    if isinstance(value, str) or not isinstance(value, Iterable):
        raise PhosPyInputError(f"{field_name} must be a sequence of count pairs")
    normalized: list[tuple[str, int]] = []
    seen: set[str] = set()
    for pair in cast(Iterable[object], value):
        if not isinstance(pair, tuple):
            raise PhosPyInputError(
                f"{field_name} must contain (non_empty_string, count) pairs"
            )
        raw_pair = cast(tuple[object, ...], pair)
        if len(raw_pair) != 2:
            raise PhosPyInputError(
                f"{field_name} must contain (non_empty_string, count) pairs"
            )
        name_raw: object = raw_pair[0]
        count_raw: object = raw_pair[1]
        if (
            isinstance(name_raw, bool)
            or not isinstance(name_raw, str)
            or not name_raw.strip()
        ):
            raise PhosPyInputError(
                f"{field_name} must contain (non_empty_string, count) pairs"
            )
        name = name_raw.strip()
        if name in seen:
            raise PhosPyInputError(f"{field_name} contains duplicate key {name!r}")
        count = _require_non_negative_int(
            count_raw,
            field_name=f"{field_name}[{name!r}]",
        )
        normalized.append((name, count))
        seen.add(name)
    return tuple(normalized)


def _validate_protein_aware_per_site_diagnostics(
    table: pd.DataFrame,
    *,
    expected_total_site_count: int,
    expected_tested_site_count: int,
    expected_withheld_site_count: int,
) -> None:
    field_name = "protein_aware_diagnostics.per_site_diagnostics"
    require_dataframe(
        table,
        field_name=field_name,
        allow_empty=False,
        error_type=PhosPyInputError,
    )
    require_non_empty_dataframe(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_string_index(
        table.index,
        field_name=f"{field_name}.index",
        error_type=PhosPyInputError,
    )
    require_unique_index(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_unique_columns(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_columns(
        table,
        field_name=field_name,
        required_columns=PROTEIN_AWARE_DIFFERENTIAL_PER_SITE_DIAGNOSTIC_COLUMNS,
        error_type=PhosPyInputError,
    )
    _reject_nuisance_inferential_columns(table=table, field_name=field_name)

    if int(table.shape[0]) != expected_total_site_count:
        raise PhosPyInputError(
            "protein_aware_diagnostics.per_site_diagnostics row count must match "
            "total_site_count"
        )
    site_key_values = [str(value) for value in table["site_key"].tolist()]
    index_values = [str(value) for value in table.index.tolist()]
    if site_key_values != index_values:
        raise PhosPyInputError(
            "protein_aware_diagnostics.per_site_diagnostics.site_key must exactly "
            "match the diagnostics index"
        )
    validate_result_status_reason_contract(table, field_name=field_name)

    status_values = np.asarray(
        table[DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str).to_numpy(),
        dtype=str,
    )
    tested_mask = np.asarray(
        status_values == DIFFERENTIAL_RESULT_STATUS_TESTED,
        dtype=bool,
    )
    withheld_mask = np.asarray(
        np.isin(status_values, DIFFERENTIAL_RESULT_WITHHELD_STATUSES),
        dtype=bool,
    )
    if int(tested_mask.sum()) != expected_tested_site_count:
        raise PhosPyInputError(
            "protein_aware_diagnostics.per_site_diagnostics tested row count must "
            "match tested_site_count"
        )
    if int(withheld_mask.sum()) != expected_withheld_site_count:
        raise PhosPyInputError(
            "protein_aware_diagnostics.per_site_diagnostics withheld row count must "
            "match withheld_site_count"
        )

    required_tested_numeric_columns = (
        "protein_covariate_raw_mean",
        "protein_covariate_raw_standard_deviation",
        "protein_covariate_centered_standard_deviation",
        "protein_augmented_design_rank",
        "protein_augmented_design_residual_degrees_of_freedom",
        "protein_augmented_design_condition_number",
        "protein_covariate_coefficient",
    )
    for column_name in required_tested_numeric_columns:
        values = table[column_name]
        finite_mask = np.isfinite(values.to_numpy(dtype=float, na_value=np.nan))
        if not bool(finite_mask[tested_mask].all()):
            raise PhosPyInputError(
                "protein_aware_diagnostics.per_site_diagnostics tested rows must "
                f"contain finite {column_name} values"
            )
    protein_coefficients = table["protein_covariate_coefficient"].to_numpy(
        dtype=object,
    )
    withheld_coefficients = protein_coefficients[withheld_mask]
    if bool(np.asarray(~pd.isna(withheld_coefficients), dtype=bool).any()):
        raise PhosPyInputError(
            "protein_aware_diagnostics.per_site_diagnostics withheld rows must not "
            "contain protein_covariate_coefficient values"
        )


def _reject_nuisance_inferential_columns(
    *,
    table: pd.DataFrame,
    field_name: str,
) -> None:
    disallowed_names = {
        "protein_covariate_p_value",
        "protein_covariate_q_value",
        "protein_covariate_adj_p_value",
        "protein_covariate_adjusted_p_value",
        "protein_covariate_p.value",
        "protein_covariate_adj.p.val",
    }
    present = [
        str(column)
        for column in table.columns
        if str(column).lower() in disallowed_names
    ]
    if present:
        raise PhosPyInputError(
            f"{field_name} must not contain inferential nuisance columns: "
            + ", ".join(present)
        )


def _require_counts_match_per_site_diagnostics(
    *,
    per_site_diagnostics: pd.DataFrame,
    status_counts: tuple[tuple[str, int], ...],
    reason_counts: tuple[tuple[str, int], ...],
) -> None:
    status_values = [
        str(value)
        for value in per_site_diagnostics[DIFFERENTIAL_RESULT_STATUS_COLUMN].tolist()
    ]
    expected_status_counts = tuple(
        (status, status_values.count(status)) for status in dict.fromkeys(status_values)
    )
    if status_counts != expected_status_counts:
        raise PhosPyInputError(
            "protein_aware_diagnostics.status_counts must match "
            "per_site_diagnostics.result_status counts"
        )

    reason_values = [
        str(value).strip()
        for value in per_site_diagnostics[
            DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN
        ].tolist()
        if str(value).strip()
    ]
    expected_reason_counts = tuple(
        (reason, reason_values.count(reason)) for reason in dict.fromkeys(reason_values)
    )
    if reason_counts != expected_reason_counts:
        raise PhosPyInputError(
            "protein_aware_diagnostics.reason_counts must match "
            "per_site_diagnostics.result_status_reason counts"
        )


def _dataframe_records_payload(frame: pd.DataFrame) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    columns = [str(column) for column in frame.columns]
    values = cast(npt.NDArray[np.object_], frame.to_numpy(dtype=object))
    for row_position in range(int(values.shape[0])):
        record: dict[str, object] = {}
        for column_position, column in enumerate(columns):
            record[column] = _json_payload(
                cast(object, values[row_position, column_position])
            )
        records.append(record)
    return records


def _json_payload(value: object) -> object:
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, object], value)
        return {str(key): _json_payload(item) for key, item in mapping.items()}
    if isinstance(value, tuple):
        values = cast(tuple[object, ...], value)
        return [_json_payload(item) for item in values]
    if isinstance(value, list):
        values = cast(list[object], value)
        return [_json_payload(item) for item in values]
    if value is None:
        return None
    item: object = (
        cast(object, value.item()) if isinstance(value, np.generic) else value
    )
    if item is pd.NA or item is pd.NaT:
        return None
    if isinstance(item, float) and not math.isfinite(item):
        return None
    return item


__all__ = [
    "DifferentialModelDiagnostics",
    "EmpiricalBayesPriorDiagnostics",
    "MeanVarianceTrendDiagnostics",
    "PROTEIN_AWARE_DIFFERENTIAL_CLAIM_STATUS_EXPERIMENTAL",
    "PROTEIN_AWARE_DIFFERENTIAL_CONDITION_NUMBER_SUMMARY_SCOPE",
    "PROTEIN_AWARE_DIFFERENTIAL_MODEL_TYPE",
    "PROTEIN_AWARE_DIFFERENTIAL_PER_SITE_DIAGNOSTIC_COLUMNS",
    "ProteinAwareDifferentialDiagnostics",
]
