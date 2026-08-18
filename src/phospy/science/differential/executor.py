"""Differential-analysis execution engine."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd
from scipy import stats

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.compound_symmetry_gls import (
    COMPOUND_SYMMETRY_GLS_STATUS_FIT,
    CompoundSymmetryGLSFit,
    fit_duplicate_correlation_gls,
)
from phospy.science.differential.duplicate_correlation import (
    estimate_duplicate_correlation_reml_consensus,
)
from phospy.science.differential.empirical_bayes import fit_empirical_bayes
from phospy.science.differential.linear_model import (
    DifferentialDesignDecomposition,
    DifferentialDesignDecompositionError,
)
from phospy.science.differential.models import (
    DifferentialAnalysisRequest,
    DifferentialComputationResult,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.models.duplicate_correlation import (
    DuplicateCorrelationConsensusResult,
)
from phospy.science.statistics.multiple_testing import adjust_p_values


@dataclass(frozen=True, slots=True)
class DuplicateCorrelationDifferentialFit:
    """Duplicate-correlation fit artifact for workflow provenance assembly."""

    computation_result: DifferentialComputationResult
    consensus: DuplicateCorrelationConsensusResult
    gls_fit: CompoundSymmetryGLSFit


class DifferentialAnalysisExecutor:
    """Run OLS fitting and empirical-Bayes moderation."""

    def run(
        self, request: DifferentialAnalysisRequest
    ) -> DifferentialComputationResult:
        matrix_aligned, design_aligned, contrast_aligned = _align_execution_inputs(
            request
        )

        design_decomposition = cast(
            DifferentialDesignDecomposition,
            request.design_decomposition,
        )
        try:
            design_decomposition.assert_matches_design(
                design_aligned.to_numpy(dtype=float),
                field_name="differential.design",
            )
        except DifferentialDesignDecompositionError as error:
            raise PhosPyInputError(
                "differential.design_decomposition does not match the aligned "
                f"differential design: {error}"
            ) from error
        residual_dof = design_decomposition.residual_degrees_of_freedom

        response = matrix_aligned.to_numpy(dtype=float).T
        try:
            linear_fit = design_decomposition.fit(response)
        except DifferentialDesignDecompositionError as error:
            raise PhosPyInputError(
                "differential model fitting failed under the shared design "
                f"decomposition: {error}"
            ) from error
        coefficients = linear_fit.coefficients
        residual_variance = linear_fit.residual_variance
        mean_intensity = np.mean(response, axis=0)

        try:
            eb_fit = fit_empirical_bayes(
                variances=residual_variance,
                residual_dof=residual_dof,
                method=request.empirical_bayes.method,
                trend=request.empirical_bayes.trend,
                winsor_tail_p=request.empirical_bayes.winsor_tail_p,
                mean_intensity=mean_intensity,
            )
        except ValueError as error:
            raise PhosPyInputError(
                "empirical-Bayes prior estimation failed for differential analysis"
            ) from error

        prior_variance = eb_fit.prior_variance
        prior_dof = eb_fit.prior_degrees_of_freedom

        total_residual_dof = residual_dof * float(matrix_aligned.shape[0])
        posterior_variance = np.empty_like(residual_variance, dtype=float)
        finite_prior_dof = np.isfinite(prior_dof)

        if np.any(finite_prior_dof):
            posterior_variance[finite_prior_dof] = (
                prior_dof[finite_prior_dof] * prior_variance[finite_prior_dof]
                + residual_dof * residual_variance[finite_prior_dof]
            ) / (prior_dof[finite_prior_dof] + residual_dof)
        if np.any(~finite_prior_dof):
            posterior_variance[~finite_prior_dof] = prior_variance[~finite_prior_dof]

        moderated_dof = np.where(
            finite_prior_dof,
            np.minimum(residual_dof + prior_dof, total_residual_dof),
            total_residual_dof,
        )
        invalid_moderated_dof = ~np.isfinite(moderated_dof) | (moderated_dof <= 0.0)
        if np.any(invalid_moderated_dof):
            raise PhosPyInputError(
                "differential analysis produced invalid moderated degrees of freedom; "
                "degrees of freedom must be finite and > 0.0; "
                f"{_preview_invalid_entries(invalid_moderated_dof, moderated_dof, row_index=matrix_aligned.index)}"
            )

        contrasts = contrast_aligned.to_numpy(dtype=float)
        contrast_effects = coefficients.T @ contrasts
        contrast_scale = design_decomposition.contrast_scales(contrasts)
        if np.any(~np.isfinite(contrast_scale)) or np.any(contrast_scale <= 0.0):
            raise PhosPyInputError(
                "differential.contrasts contains non-estimable or zero-length contrast "
                "vectors under the provided design matrix"
            )

        row_index = matrix_aligned.index.copy()
        contrast_tables: dict[str, pd.DataFrame] = {}
        posterior_sd = np.sqrt(posterior_variance)

        for column_idx, contrast_name in enumerate(contrast_aligned.columns):
            log_fc = contrast_effects[:, column_idx]
            standard_error = posterior_sd * contrast_scale[column_idx]
            invalid_standard_error = ~np.isfinite(standard_error) | (
                standard_error <= 0.0
            )
            if np.any(invalid_standard_error):
                raise PhosPyInputError(
                    "differential analysis produced unstable standard errors for "
                    f"contrast {contrast_name!r}; standard errors must be finite "
                    "and > 0.0; "
                    f"{_preview_invalid_entries(invalid_standard_error, standard_error, row_index=row_index)}"
                )
            moderated_t = log_fc / standard_error
            p_values = np.asarray(
                2.0 * stats.t.sf(np.abs(moderated_t), df=moderated_dof),
                dtype=float,
            )
            invalid_p_values = (
                ~np.isfinite(p_values) | (p_values < 0.0) | (p_values > 1.0)
            )
            if np.any(invalid_p_values):
                raise PhosPyInputError(
                    "differential analysis produced invalid p-values for contrast "
                    f"{contrast_name!r}; P.Value must be finite and within [0, 1]; "
                    f"{_preview_invalid_entries(invalid_p_values, p_values, row_index=row_index)}"
                )
            adjusted = adjust_p_values(
                p_values,
                method=request.multiple_testing_method,
            )
            contrast_tables[str(contrast_name)] = pd.DataFrame(
                {
                    "logFC": log_fc.astype(float),
                    "t": moderated_t.astype(float),
                    "P.Value": p_values,
                    "adj.P.Val": adjusted.astype(float),
                },
                index=row_index.copy(),
            )

        residual_variance_series = pd.Series(
            residual_variance.astype(float),
            index=row_index.copy(),
            name="residual_variance",
        )
        posterior_variance_series = pd.Series(
            posterior_variance.astype(float),
            index=row_index.copy(),
            name="posterior_residual_variance",
        )
        prior_variance_series = pd.Series(
            prior_variance.astype(float),
            index=row_index.copy(),
            name="prior_residual_variance",
        )
        prior_dof_series = pd.Series(
            prior_dof.astype(float),
            index=row_index.copy(),
            name="prior_degrees_of_freedom",
        )
        prior_diagnostics = EmpiricalBayesPriorDiagnostics(
            method=request.empirical_bayes.method,
            robust=request.empirical_bayes.method == "robust",
            trend=request.empirical_bayes.trend,
            winsor_tail_p=request.empirical_bayes.winsor_tail_p,
            base_prior_variance=eb_fit.base_prior_variance,
            base_prior_degrees_of_freedom=eb_fit.base_prior_degrees_of_freedom,
            robust_outlier_count=eb_fit.robust_outlier_count,
            robust_outlier_fraction=eb_fit.robust_outlier_fraction,
            winsorized_low_count=eb_fit.winsorized_low_count,
            winsorized_high_count=eb_fit.winsorized_high_count,
            prior_variance=prior_variance_series,
            prior_degrees_of_freedom=prior_dof_series,
            _assume_owned=True,
        )
        trend_diagnostics: MeanVarianceTrendDiagnostics | None = None
        if request.empirical_bayes.trend:
            trend_diagnostics = MeanVarianceTrendDiagnostics(
                mean_intensity=pd.Series(
                    eb_fit.mean_intensity,
                    index=row_index.copy(),
                    name="mean_intensity",
                ),
                log_residual_variance=pd.Series(
                    eb_fit.log_residual_variance,
                    index=row_index.copy(),
                    name="log_residual_variance",
                ),
                fitted_log_prior_variance=pd.Series(
                    eb_fit.fitted_log_prior_variance,
                    index=row_index.copy(),
                    name="fitted_log_prior_variance",
                ),
                _assume_owned=True,
            )

        return DifferentialComputationResult._from_owned(
            design_decomposition=design_decomposition,
            residual_variance=residual_variance_series,
            posterior_residual_variance=posterior_variance_series,
            prior_residual_variance=prior_variance_series,
            prior_degrees_of_freedom_series_value=prior_dof_series,
            prior_variance=float(np.nanmedian(prior_variance)),
            prior_degrees_of_freedom=float(np.nanmedian(prior_dof)),
            residual_degrees_of_freedom=float(residual_dof),
            empirical_bayes_method=request.empirical_bayes.method,
            empirical_bayes_robust=request.empirical_bayes.method == "robust",
            empirical_bayes_trend=request.empirical_bayes.trend,
            prior_diagnostics=prior_diagnostics,
            mean_variance_trend_diagnostics=trend_diagnostics,
            contrast_tables=contrast_tables,
        )


class DuplicateCorrelationDifferentialAnalysisExecutor:
    """Run duplicate-correlation REML, GLS, and empirical-Bayes moderation."""

    def run(
        self,
        request: DifferentialAnalysisRequest,
        *,
        block_ids: Sequence[object],
    ) -> DuplicateCorrelationDifferentialFit:
        matrix_aligned, design_aligned, contrast_aligned = _align_execution_inputs(
            request
        )
        design_decomposition = cast(
            DifferentialDesignDecomposition,
            request.design_decomposition,
        )
        try:
            design_decomposition.assert_matches_design(
                design_aligned.to_numpy(dtype=float),
                field_name="differential.design",
            )
        except DifferentialDesignDecompositionError as error:
            raise PhosPyInputError(
                "differential.design_decomposition does not match the aligned "
                f"differential duplicate-correlation design: {error}"
            ) from error

        residual_dof = float(design_decomposition.residual_degrees_of_freedom)
        matrix_values = matrix_aligned.to_numpy(dtype=float)
        design_values = design_aligned.to_numpy(dtype=float)
        contrast_values = contrast_aligned.to_numpy(dtype=float)
        feature_ids = tuple(str(label) for label in matrix_aligned.index)
        coefficient_names = tuple(str(label) for label in design_aligned.columns)
        contrast_names = tuple(str(label) for label in contrast_aligned.columns)

        consensus = estimate_duplicate_correlation_reml_consensus(
            matrix=matrix_values,
            design=design_values,
            block_ids=block_ids,
            feature_ids=feature_ids,
            design_column_names=coefficient_names,
            block_id_field_name="block_id",
        )
        if not consensus.success or consensus.consensus_correlation is None:
            raise PhosPyInputError(
                "duplicate-correlation REML consensus estimation failed; "
                f"failure_reason={consensus.failure_reason}, "
                f"eligible_feature_count={consensus.eligible_feature_count}, "
                f"estimated_feature_count={consensus.estimated_feature_count}"
            )

        gls_fit = fit_duplicate_correlation_gls(
            matrix=matrix_values,
            design=design_values,
            block_ids=block_ids,
            consensus_correlation=consensus.consensus_correlation,
            feature_ids=feature_ids,
            coefficient_names=coefficient_names,
            contrasts=contrast_values,
            contrast_names=contrast_names,
        )
        _require_successful_duplicate_correlation_gls_fit(
            gls_fit=gls_fit,
            row_index=matrix_aligned.index,
        )
        computation_result = _duplicate_correlation_computation_result(
            request=request,
            design_decomposition=design_decomposition,
            residual_dof=residual_dof,
            residual_variance=gls_fit.residual_variance,
            mean_intensity=gls_fit.average_expression,
            contrast_effects=cast(np.ndarray, gls_fit.contrast_coefficients).T,
            contrast_stdev_unscaled=cast(np.ndarray, gls_fit.contrast_stdev_unscaled).T,
            contrast_names=contrast_names,
            row_index=matrix_aligned.index.copy(),
        )
        return DuplicateCorrelationDifferentialFit(
            computation_result=computation_result,
            consensus=consensus,
            gls_fit=gls_fit,
        )


def _align_execution_inputs(
    request: DifferentialAnalysisRequest,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    matrix = cast(pd.DataFrame, request.matrix)
    design_frame = cast(pd.DataFrame, request.design.frame)
    contrast_frame = cast(pd.DataFrame, request.contrasts.frame)

    if matrix.shape[1] != design_frame.shape[0]:
        raise PhosPyInputError(
            "differential.matrix.columns count must match differential.design rows; "
            f"matrix_columns={int(matrix.shape[1])}, "
            f"design_rows={int(design_frame.shape[0])}"
        )

    matrix_samples = pd.Index(matrix.columns)
    design_samples = pd.Index(design_frame.index)
    if not matrix_samples.equals(design_samples) and (
        not matrix_samples.isin(design_samples).all()
        or not design_samples.isin(matrix_samples).all()
    ):
        missing_matrix_samples = [
            str(sample) for sample in design_samples if sample not in matrix_samples
        ]
        extra_matrix_samples = [
            str(sample) for sample in matrix_samples if sample not in design_samples
        ]
        raise PhosPyInputError(
            "differential.matrix.columns must match differential.design.index "
            "exactly as execution sample labels; "
            f"missing_matrix_samples={missing_matrix_samples}, "
            f"extra_matrix_samples={extra_matrix_samples}"
        )

    design_terms = pd.Index(design_frame.columns)
    contrast_terms = pd.Index(contrast_frame.index)
    if not design_terms.equals(contrast_terms) and (
        not design_terms.isin(contrast_terms).all()
        or not contrast_terms.isin(design_terms).all()
    ):
        missing_contrast_terms = [
            str(term) for term in design_terms if term not in contrast_terms
        ]
        extra_contrast_terms = [
            str(term) for term in contrast_terms if term not in design_terms
        ]
        raise PhosPyInputError(
            "differential.contrasts.index must match differential.design.columns "
            "exactly as execution design-term labels; "
            f"missing_contrast_terms={missing_contrast_terms}, "
            f"extra_contrast_terms={extra_contrast_terms}"
        )

    if int(contrast_frame.shape[1]) < 1:
        raise PhosPyInputError(
            "differential.contrasts must contain at least one contrast column"
        )

    matrix_aligned = cast(pd.DataFrame, matrix.loc[:, list(design_frame.index)])
    contrast_aligned = cast(
        pd.DataFrame,
        contrast_frame.loc[list(design_frame.columns), :],
    )
    return matrix_aligned, design_frame, contrast_aligned


def _duplicate_correlation_computation_result(
    *,
    request: DifferentialAnalysisRequest,
    design_decomposition: DifferentialDesignDecomposition,
    residual_dof: float,
    residual_variance: np.ndarray,
    mean_intensity: np.ndarray,
    contrast_effects: np.ndarray,
    contrast_stdev_unscaled: np.ndarray,
    contrast_names: tuple[str, ...],
    row_index: pd.Index,
) -> DifferentialComputationResult:
    residual_variance = np.asarray(residual_variance, dtype=float)
    mean_intensity = np.asarray(mean_intensity, dtype=float)
    invalid_residual_variance = ~np.isfinite(residual_variance) | (
        residual_variance <= 0.0
    )
    if np.any(invalid_residual_variance):
        raise PhosPyInputError(
            "duplicate-correlation GLS produced invalid residual variances; "
            "residual variances must be finite and > 0.0; "
            f"{_preview_invalid_entries(invalid_residual_variance, residual_variance, row_index=row_index)}"
        )
    try:
        eb_fit = fit_empirical_bayes(
            variances=residual_variance,
            residual_dof=residual_dof,
            method=request.empirical_bayes.method,
            trend=request.empirical_bayes.trend,
            winsor_tail_p=request.empirical_bayes.winsor_tail_p,
            mean_intensity=mean_intensity,
        )
    except ValueError as error:
        raise PhosPyInputError(
            "empirical-Bayes prior estimation failed for duplicate-correlation "
            "differential analysis"
        ) from error

    prior_variance = eb_fit.prior_variance
    prior_dof = eb_fit.prior_degrees_of_freedom

    total_residual_dof = residual_dof * float(row_index.size)
    posterior_variance = np.empty_like(residual_variance, dtype=float)
    finite_prior_dof = np.isfinite(prior_dof)

    if np.any(finite_prior_dof):
        posterior_variance[finite_prior_dof] = (
            prior_dof[finite_prior_dof] * prior_variance[finite_prior_dof]
            + residual_dof * residual_variance[finite_prior_dof]
        ) / (prior_dof[finite_prior_dof] + residual_dof)
    if np.any(~finite_prior_dof):
        posterior_variance[~finite_prior_dof] = prior_variance[~finite_prior_dof]

    moderated_dof = np.where(
        finite_prior_dof,
        np.minimum(residual_dof + prior_dof, total_residual_dof),
        total_residual_dof,
    )
    invalid_moderated_dof = ~np.isfinite(moderated_dof) | (moderated_dof <= 0.0)
    if np.any(invalid_moderated_dof):
        raise PhosPyInputError(
            "duplicate-correlation differential analysis produced invalid moderated "
            "degrees of freedom; degrees of freedom must be finite and > 0.0; "
            f"{_preview_invalid_entries(invalid_moderated_dof, moderated_dof, row_index=row_index)}"
        )

    contrast_effects = np.asarray(contrast_effects, dtype=float)
    contrast_stdev_unscaled = np.asarray(contrast_stdev_unscaled, dtype=float)
    expected_shape = (int(row_index.size), len(contrast_names))
    if tuple(contrast_effects.shape) != expected_shape:
        raise PhosPyInputError(
            "duplicate-correlation GLS contrast coefficients are not aligned to "
            "feature rows and contrast columns"
        )
    if tuple(contrast_stdev_unscaled.shape) != expected_shape:
        raise PhosPyInputError(
            "duplicate-correlation GLS contrast standard errors are not aligned to "
            "feature rows and contrast columns"
        )

    contrast_tables: dict[str, pd.DataFrame] = {}
    posterior_sd = np.sqrt(posterior_variance)
    for column_idx, contrast_name in enumerate(contrast_names):
        log_fc = contrast_effects[:, column_idx]
        standard_error = posterior_sd * contrast_stdev_unscaled[:, column_idx]
        invalid_standard_error = ~np.isfinite(standard_error) | (standard_error <= 0.0)
        if np.any(invalid_standard_error):
            raise PhosPyInputError(
                "duplicate-correlation differential analysis produced unstable "
                f"standard errors for contrast {contrast_name!r}; standard errors "
                "must be finite and > 0.0; "
                f"{_preview_invalid_entries(invalid_standard_error, standard_error, row_index=row_index)}"
            )
        moderated_t = log_fc / standard_error
        p_values = np.asarray(
            2.0 * stats.t.sf(np.abs(moderated_t), df=moderated_dof),
            dtype=float,
        )
        invalid_p_values = ~np.isfinite(p_values) | (p_values < 0.0) | (p_values > 1.0)
        if np.any(invalid_p_values):
            raise PhosPyInputError(
                "duplicate-correlation differential analysis produced invalid "
                f"p-values for contrast {contrast_name!r}; P.Value must be finite "
                "and within [0, 1]; "
                f"{_preview_invalid_entries(invalid_p_values, p_values, row_index=row_index)}"
            )
        adjusted = adjust_p_values(
            p_values,
            method=request.multiple_testing_method,
        )
        contrast_tables[str(contrast_name)] = pd.DataFrame(
            {
                "logFC": log_fc.astype(float),
                "t": moderated_t.astype(float),
                "P.Value": p_values,
                "adj.P.Val": adjusted.astype(float),
            },
            index=row_index.copy(),
        )

    residual_variance_series = pd.Series(
        residual_variance.astype(float),
        index=row_index.copy(),
        name="residual_variance",
    )
    posterior_variance_series = pd.Series(
        posterior_variance.astype(float),
        index=row_index.copy(),
        name="posterior_residual_variance",
    )
    prior_variance_series = pd.Series(
        prior_variance.astype(float),
        index=row_index.copy(),
        name="prior_residual_variance",
    )
    prior_dof_series = pd.Series(
        prior_dof.astype(float),
        index=row_index.copy(),
        name="prior_degrees_of_freedom",
    )
    prior_diagnostics = EmpiricalBayesPriorDiagnostics(
        method=request.empirical_bayes.method,
        robust=request.empirical_bayes.method == "robust",
        trend=request.empirical_bayes.trend,
        winsor_tail_p=request.empirical_bayes.winsor_tail_p,
        base_prior_variance=eb_fit.base_prior_variance,
        base_prior_degrees_of_freedom=eb_fit.base_prior_degrees_of_freedom,
        robust_outlier_count=eb_fit.robust_outlier_count,
        robust_outlier_fraction=eb_fit.robust_outlier_fraction,
        winsorized_low_count=eb_fit.winsorized_low_count,
        winsorized_high_count=eb_fit.winsorized_high_count,
        prior_variance=prior_variance_series,
        prior_degrees_of_freedom=prior_dof_series,
        _assume_owned=True,
    )
    trend_diagnostics: MeanVarianceTrendDiagnostics | None = None
    if request.empirical_bayes.trend:
        trend_diagnostics = MeanVarianceTrendDiagnostics(
            mean_intensity=pd.Series(
                eb_fit.mean_intensity,
                index=row_index.copy(),
                name="mean_intensity",
            ),
            log_residual_variance=pd.Series(
                eb_fit.log_residual_variance,
                index=row_index.copy(),
                name="log_residual_variance",
            ),
            fitted_log_prior_variance=pd.Series(
                eb_fit.fitted_log_prior_variance,
                index=row_index.copy(),
                name="fitted_log_prior_variance",
            ),
            _assume_owned=True,
        )

    return DifferentialComputationResult._from_owned(
        design_decomposition=design_decomposition,
        residual_variance=residual_variance_series,
        posterior_residual_variance=posterior_variance_series,
        prior_residual_variance=prior_variance_series,
        prior_degrees_of_freedom_series_value=prior_dof_series,
        prior_variance=float(np.nanmedian(prior_variance)),
        prior_degrees_of_freedom=float(np.nanmedian(prior_dof)),
        residual_degrees_of_freedom=float(residual_dof),
        empirical_bayes_method=request.empirical_bayes.method,
        empirical_bayes_robust=request.empirical_bayes.method == "robust",
        empirical_bayes_trend=request.empirical_bayes.trend,
        prior_diagnostics=prior_diagnostics,
        mean_variance_trend_diagnostics=trend_diagnostics,
        contrast_tables=contrast_tables,
    )


def _require_successful_duplicate_correlation_gls_fit(
    *,
    gls_fit: CompoundSymmetryGLSFit,
    row_index: pd.Index,
) -> None:
    if gls_fit.contrast_coefficients is None or gls_fit.contrast_stdev_unscaled is None:
        raise PhosPyInputError(
            "duplicate-correlation GLS fit did not produce contrast-level outputs"
        )
    failed_statuses = sorted(
        {
            str(status)
            for status in gls_fit.feature_fit_statuses
            if str(status) != COMPOUND_SYMMETRY_GLS_STATUS_FIT
        }
    )
    if failed_statuses:
        examples = [
            str(feature_id)
            for feature_id, status in zip(
                row_index,
                gls_fit.feature_fit_statuses,
                strict=True,
            )
            if str(status) != COMPOUND_SYMMETRY_GLS_STATUS_FIT
        ][:5]
        raise PhosPyInputError(
            "duplicate-correlation GLS fitting failed for one or more features; "
            "there is no fallback to fixed_block, ordinary least squares, or "
            "correlation zero. Feature statuses: "
            + ", ".join(failed_statuses)
            + "; examples="
            + ", ".join(examples)
        )


def _preview_invalid_entries(
    invalid_mask: np.ndarray,
    values: np.ndarray,
    *,
    row_index: pd.Index,
) -> str:
    invalid_positions = np.flatnonzero(invalid_mask)
    preview = ", ".join(
        f"({row_index[position]!r}, {values[position]:.6g})"
        for position in invalid_positions[:3]
    )
    suffix = (
        ""
        if invalid_positions.size <= 3
        else f", +{int(invalid_positions.size - 3)} more"
    )
    return f"invalid values: {preview}{suffix}; invalid_entry_count={int(invalid_positions.size)}"


__all__ = [
    "DifferentialAnalysisExecutor",
    "DuplicateCorrelationDifferentialAnalysisExecutor",
    "DuplicateCorrelationDifferentialFit",
]
