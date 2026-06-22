"""Differential-analysis execution engine."""

from __future__ import annotations

from typing import cast

import numpy as np
import pandas as pd
from scipy import stats

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.empirical_bayes import fit_empirical_bayes
from phospy.science.differential.models import (
    DifferentialAnalysisRequest,
    DifferentialComputationResult,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.multiple_testing import benjamini_hochberg


class DifferentialAnalysisExecutor:
    """Run OLS fitting and empirical-Bayes moderation."""

    def run(
        self, request: DifferentialAnalysisRequest
    ) -> DifferentialComputationResult:
        matrix_aligned, design_aligned, contrast_aligned = _align_execution_inputs(
            request
        )

        design_values = design_aligned.to_numpy(dtype=float)
        if design_values.ndim != 2:
            raise PhosPyInputError("differential.design must be two-dimensional")
        design_shape = cast(tuple[int, int], design_values.shape)
        sample_count = int(design_shape[0])
        coefficient_count = int(design_shape[1])
        rank = int(np.linalg.matrix_rank(design_values))
        if rank < coefficient_count:
            raise PhosPyInputError(
                "differential.design must be full column rank for moderated-contrast "
                f"analysis; rank={rank}, columns={coefficient_count}"
            )
        residual_dof = float(sample_count - rank)
        if residual_dof <= 0.0:
            raise PhosPyInputError(
                "differential analysis residual degrees of freedom must be positive; "
                f"samples={sample_count}, rank={rank}, residual_dof={residual_dof}"
            )

        response = matrix_aligned.to_numpy(dtype=float).T
        xtx_inv = np.linalg.pinv(design_values.T @ design_values)
        coefficients = xtx_inv @ design_values.T @ response
        residuals = response - design_values @ coefficients
        rss = np.sum(residuals**2, axis=0)
        residual_variance = rss / residual_dof
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
        contrast_covariance = contrasts.T @ xtx_inv @ contrasts
        contrast_scale = np.sqrt(np.diag(contrast_covariance))
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
            adjusted = benjamini_hochberg(p_values)
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
