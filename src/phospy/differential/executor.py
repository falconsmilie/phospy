"""Differential-analysis execution engine."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from phospy.differential.empirical_bayes import fit_f_dist
from phospy.differential.models import (
    EMPIRICAL_BAYES_METHOD_STANDARD,
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
)
from phospy.differential.multiple_testing import benjamini_hochberg
from phospy.errors.input import PhosPyInputError


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
        sample_count, coefficient_count = design_values.shape
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

        if request.empirical_bayes.method != EMPIRICAL_BAYES_METHOD_STANDARD:
            raise PhosPyInputError(
                f"unsupported empirical_bayes.method={request.empirical_bayes.method!r}"
            )
        prior_variance, prior_dof = fit_f_dist(
            residual_variance,
            residual_dof=residual_dof,
        )
        if not np.isfinite(prior_variance) or prior_variance <= 0.0:
            raise PhosPyInputError(
                "empirical-Bayes prior variance estimation failed for "
                "differential analysis"
            )
        if np.isnan(prior_dof) or prior_dof < 0.0:
            raise PhosPyInputError(
                "empirical-Bayes prior degrees of freedom estimation failed for "
                "differential analysis"
            )

        total_residual_dof = residual_dof * float(matrix_aligned.shape[0])
        if np.isfinite(prior_dof):
            posterior_variance = (
                prior_dof * prior_variance + residual_dof * residual_variance
            ) / (prior_dof + residual_dof)
            moderated_dof = min(residual_dof + prior_dof, total_residual_dof)
        else:
            posterior_variance = np.full_like(residual_variance, prior_variance)
            moderated_dof = total_residual_dof

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
            p_values = 2.0 * stats.t.sf(np.abs(moderated_t), df=moderated_dof)
            adjusted = benjamini_hochberg(p_values)
            contrast_tables[str(contrast_name)] = pd.DataFrame(
                {
                    "logFC": log_fc.astype(float),
                    "t": moderated_t.astype(float),
                    "P.Value": p_values.astype(float),
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
        return DifferentialAnalysisResult._from_owned(
            residual_variance=residual_variance_series,
            posterior_residual_variance=posterior_variance_series,
            prior_variance=float(prior_variance),
            prior_degrees_of_freedom=float(prior_dof),
            residual_degrees_of_freedom=float(residual_dof),
            contrast_tables=contrast_tables,
        )
