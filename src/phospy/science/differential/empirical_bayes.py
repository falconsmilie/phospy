"""Empirical-Bayes moderation helpers for differential OLS."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats
from scipy.special import digamma, polygamma


@dataclass(frozen=True, slots=True)
class EmpiricalBayesFit:
    """Estimated empirical-Bayes prior parameters and diagnostics."""

    prior_variance: np.ndarray
    prior_degrees_of_freedom: np.ndarray
    base_prior_variance: float
    base_prior_degrees_of_freedom: float
    robust_outlier_count: int
    robust_outlier_fraction: float
    winsorized_low_count: int
    winsorized_high_count: int
    mean_intensity: np.ndarray | None = None
    log_residual_variance: np.ndarray | None = None
    fitted_log_prior_variance: np.ndarray | None = None


def fit_f_dist(
    variances: np.ndarray,
    *,
    residual_dof: float,
) -> tuple[float, float]:
    """Estimate scaled-F hyperparameters following limma's `fitFDist` moments."""

    if variances.size == 0:
        return float("nan"), float("nan")
    if variances.size == 1:
        return float(variances[0]), 0.0

    variances = np.asarray(variances, dtype=float)
    df1 = np.asarray(residual_dof, dtype=float)
    if not np.isfinite(df1) or df1 <= 1e-15:
        return float("nan"), float("nan")

    ok = np.isfinite(variances) & (variances >= -1e-15)
    if not np.any(ok):
        return float("nan"), float("nan")
    variances = variances[ok]
    if variances.size == 1:
        return float(variances[0]), 0.0

    variances = np.maximum(variances, 0.0)
    median_variance = float(np.median(variances))
    if median_variance == 0.0:
        median_variance = 1.0
    variances = np.maximum(variances, 1e-5 * median_variance)

    z = np.log(variances)
    e = z - digamma(df1 / 2.0) + np.log(df1 / 2.0)
    emean = float(np.mean(e))
    evar = float(np.sum((e - emean) ** 2) / float(variances.size - 1))
    evar = evar - float(polygamma(1, df1 / 2.0))

    if evar > 0.0:
        df2 = 2.0 * trigamma_inverse(evar)
        s20 = float(np.exp(emean + digamma(df2 / 2.0) - np.log(df2 / 2.0)))
        return s20, float(df2)

    return float(np.mean(variances)), float("inf")


def fit_empirical_bayes(
    *,
    variances: np.ndarray,
    residual_dof: float,
    method: str,
    trend: bool,
    winsor_tail_p: tuple[float, float],
    mean_intensity: np.ndarray | None = None,
) -> EmpiricalBayesFit:
    """Estimate prior variance/df with optional robust and trend modes."""

    residual_dof = float(residual_dof)
    if not np.isfinite(residual_dof) or residual_dof <= 0.0:
        raise ValueError(
            "residual_dof must be finite and > 0.0 for empirical-Bayes moderation"
        )

    variances = np.asarray(variances, dtype=float)
    if variances.ndim != 1:
        raise ValueError("variances must be one-dimensional")
    if variances.size == 0:
        return EmpiricalBayesFit(
            prior_variance=np.array([], dtype=float),
            prior_degrees_of_freedom=np.array([], dtype=float),
            base_prior_variance=float("nan"),
            base_prior_degrees_of_freedom=float("nan"),
            robust_outlier_count=0,
            robust_outlier_fraction=0.0,
            winsorized_low_count=0,
            winsorized_high_count=0,
        )

    variances = _stabilize_variances(variances)
    log_variances = np.log(variances)

    if trend:
        if mean_intensity is None:
            raise ValueError("mean_intensity is required when trend=True")
        mean_intensity = np.asarray(mean_intensity, dtype=float)
        if mean_intensity.shape != variances.shape:
            raise ValueError("mean_intensity must match variances length")
        trend_component = _fit_mean_variance_trend(mean_intensity, log_variances)
        trend_component = trend_component - float(np.mean(trend_component))
    else:
        trend_component = np.zeros_like(log_variances)
        mean_intensity = None

    detrended_log = log_variances - trend_component
    robust_outlier_count = 0
    robust_outlier_fraction = 0.0
    winsorized_low_count = 0
    winsorized_high_count = 0

    if method == "robust":
        (
            detrended_for_fit,
            winsorized_low_count,
            winsorized_high_count,
            lower_bound,
            upper_bound,
        ) = _winsorize_log_values(
            detrended_log,
            left_tail_p=winsor_tail_p[0],
            right_tail_p=winsor_tail_p[1],
        )
    else:
        detrended_for_fit = detrended_log
        lower_bound = float("nan")
        upper_bound = float("nan")

    base_prior_variance, base_prior_dof = fit_f_dist(
        np.exp(detrended_for_fit),
        residual_dof=residual_dof,
    )
    if not np.isfinite(base_prior_variance) or base_prior_variance <= 0.0:
        raise ValueError("failed to estimate prior variance")
    if np.isnan(base_prior_dof) or base_prior_dof < 0.0:
        raise ValueError("failed to estimate prior degrees of freedom")

    prior_variance = np.exp(trend_component) * float(base_prior_variance)
    prior_dof = np.full_like(variances, float(base_prior_dof))

    if method == "robust" and np.isfinite(base_prior_dof) and base_prior_dof > 0.0:
        f_stat = variances / np.maximum(prior_variance, np.finfo(float).tiny)
        log_tail_p = stats.f.logsf(f_stat, dfn=float(residual_dof), dfd=base_prior_dof)
        ranks = stats.rankdata(f_stat, method="ordinal")
        log_empirical_tail = np.log(variances.size - ranks + 0.5) - np.log(
            variances.size
        )
        log_prob_not_outlier = np.minimum(log_tail_p - log_empirical_tail, 0.0)
        prob_not_outlier = np.exp(log_prob_not_outlier)
        prob_outlier = -np.expm1(log_prob_not_outlier)
        robust_outlier_count = int(np.sum(log_prob_not_outlier < 0.0))
        robust_outlier_fraction = robust_outlier_count / float(variances.size)

        if robust_outlier_count > 0:
            min_log_tail = float(np.min(log_tail_p))
            if np.isneginf(min_log_tail):
                df2_outlier = 0.0
                prior_dof = prob_not_outlier * base_prior_dof
            else:
                df2_outlier = float(np.log(0.5) / min_log_tail * base_prior_dof)
                max_f = float(np.max(f_stat))
                new_log_tail = float(
                    stats.f.logsf(
                        max_f,
                        dfn=float(residual_dof),
                        dfd=max(df2_outlier, np.finfo(float).tiny),
                    )
                )
                if np.isfinite(new_log_tail) and new_log_tail < 0.0:
                    df2_outlier = float(np.log(0.5) / new_log_tail * df2_outlier)
                prior_dof = prob_not_outlier * base_prior_dof + prob_outlier * max(
                    df2_outlier,
                    0.0,
                )

            # Keep shrunk df monotone in tail probability for stability.
            order = np.argsort(log_tail_p)
            ordered = prior_dof[order]
            averaged = np.cumsum(ordered) / np.arange(1, ordered.size + 1, dtype=float)
            min_idx = int(np.argmin(averaged))
            ordered[: min_idx + 1] = averaged[min_idx]
            prior_dof[order] = np.maximum.accumulate(ordered)

            if np.isfinite(lower_bound) and np.isfinite(upper_bound):
                outside = (detrended_log < lower_bound) | (detrended_log > upper_bound)
                robust_outlier_count = max(robust_outlier_count, int(np.sum(outside)))
                robust_outlier_fraction = robust_outlier_count / float(variances.size)

    invalid_prior_variance = ~np.isfinite(prior_variance) | (prior_variance <= 0.0)
    if np.any(invalid_prior_variance):
        raise ValueError(
            "empirical-Bayes prior variance became invalid; expected finite values > 0.0"
        )
    invalid_prior_dof = np.isnan(prior_dof) | (prior_dof < 0.0)
    if np.any(invalid_prior_dof):
        raise ValueError(
            "empirical-Bayes prior degrees of freedom became invalid; expected values >= 0.0 or +inf"
        )

    return EmpiricalBayesFit(
        prior_variance=prior_variance.astype(float),
        prior_degrees_of_freedom=prior_dof.astype(float),
        base_prior_variance=float(base_prior_variance),
        base_prior_degrees_of_freedom=float(base_prior_dof),
        robust_outlier_count=robust_outlier_count,
        robust_outlier_fraction=float(robust_outlier_fraction),
        winsorized_low_count=winsorized_low_count,
        winsorized_high_count=winsorized_high_count,
        mean_intensity=mean_intensity,
        log_residual_variance=log_variances,
        fitted_log_prior_variance=np.log(prior_variance) if trend else None,
    )


def trigamma_inverse(value: float) -> float:
    """Solve `trigamma(y) == value` by Newton iteration."""

    if not np.isfinite(value):
        return value
    if value < 0:
        return float("nan")
    if value > 1e7:
        return 1.0 / value
    if value < 1e-6:
        return 1.0 / value

    y = 0.5 + 1.0 / value
    for _ in range(50):
        tri = float(polygamma(1, y))
        dif = tri * (1.0 - tri / value) / float(polygamma(2, y))
        y = y + dif
        if abs(dif / y) < 1e-8:
            break
    return float(y)


def _stabilize_variances(variances: np.ndarray) -> np.ndarray:
    variances = np.maximum(variances, 0.0)
    finite = variances[np.isfinite(variances)]
    if finite.size == 0:
        return np.full_like(variances, 1e-8)
    median_variance = float(np.median(finite))
    if median_variance <= 0.0:
        median_variance = 1.0
    floor = 1e-12 * median_variance
    clipped = np.where(np.isfinite(variances), variances, median_variance)
    return np.maximum(clipped, floor)


def _winsorize_log_values(
    values: np.ndarray,
    *,
    left_tail_p: float,
    right_tail_p: float,
) -> tuple[np.ndarray, int, int, float, float]:
    lower = float(np.quantile(values, left_tail_p))
    upper = float(np.quantile(values, 1.0 - right_tail_p))
    winsorized = np.minimum(np.maximum(values, lower), upper)
    low_count = int(np.sum(values < lower))
    high_count = int(np.sum(values > upper))
    return winsorized, low_count, high_count, lower, upper


def _fit_mean_variance_trend(
    mean_intensity: np.ndarray,
    log_variances: np.ndarray,
    *,
    span: float = 0.4,
) -> np.ndarray:
    """Fit a smooth mean-variance trend using local linear regression."""

    if mean_intensity.size < 3:
        return np.full_like(log_variances, float(np.mean(log_variances)))

    order = np.argsort(mean_intensity)
    x = mean_intensity[order]
    y = log_variances[order]
    n = x.size
    window = min(n, max(5, int(np.ceil(span * n))))
    fitted = np.empty(n, dtype=float)

    for idx in range(n):
        distance = np.abs(x - x[idx])
        bandwidth = float(np.partition(distance, window - 1)[window - 1])
        if bandwidth <= 0.0:
            same = distance == 0.0
            fitted[idx] = float(np.mean(y[same]))
            continue
        u = distance / bandwidth
        weights = np.where(u < 1.0, (1.0 - u**3) ** 3, 0.0)
        x_centered = x - x[idx]
        xw = x_centered * weights
        s0 = float(np.sum(weights))
        s1 = float(np.sum(xw))
        s2 = float(np.sum(x_centered * xw))
        t0 = float(np.sum(weights * y))
        t1 = float(np.sum(xw * y))
        determinant = s0 * s2 - s1 * s1
        if determinant <= 1e-12 or not np.isfinite(determinant):
            fitted[idx] = t0 / max(s0, 1e-12)
            continue
        intercept = (t0 * s2 - t1 * s1) / determinant
        fitted[idx] = intercept

    unsorted = np.empty_like(fitted)
    unsorted[order] = fitted
    finite = np.isfinite(unsorted)
    if not finite.all():
        fill = float(np.mean(log_variances))
        unsorted[~finite] = fill
    return unsorted
