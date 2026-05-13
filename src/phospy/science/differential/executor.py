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
    DifferentialAnalysisResult,
    EmpiricalBayesPriorDiagnostics,
    MeanVarianceTrendDiagnostics,
)
from phospy.science.differential.multiple_testing import benjamini_hochberg


class DifferentialAnalysisExecutor:
    """Run OLS fitting and limma-style empirical-Bayes moderation."""

    def run(self, request: DifferentialAnalysisRequest) -> DifferentialAnalysisResult:
        matrix = request.matrix
        design_frame = request.design.frame
        contrast_frame = request.contrasts.frame

        matrix_aligned = matrix.loc[:, design_frame.index]
        design_aligned = design_frame
        contrast_aligned = contrast_frame.loc[design_aligned.columns]

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

        contrasts = contrast_aligned.to_numpy(dtype=float)
        contrast_effects = response.T @ design_values @ xtx_inv @ contrasts
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
            moderated_t = log_fc / standard_error
            p_values = np.asarray(
                2.0 * stats.t.sf(np.abs(moderated_t), df=moderated_dof),
                dtype=float,
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

        return DifferentialAnalysisResult._from_owned(
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
