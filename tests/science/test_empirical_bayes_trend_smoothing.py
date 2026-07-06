from __future__ import annotations

import numpy as np

from phospy.science.differential.empirical_bayes import (
    fit_empirical_bayes,
    fit_f_dist,
)
from tests.support.performance_contracts import DEFAULT_PERFORMANCE_SEED


def _build_trend_inputs(*, n_features: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(DEFAULT_PERFORMANCE_SEED + int(n_features))
    mean_intensity = np.sort(rng.beta(2.0, 5.0, int(n_features))) * 6.0 + 7.5
    feature_phase = np.linspace(0.0, 2.0 * np.pi, int(n_features), dtype=float)
    trend = (
        0.18 * np.sin(mean_intensity * 1.4)
        + 0.08 * np.cos(feature_phase * 3.0)
        - 0.05 * (mean_intensity - float(np.mean(mean_intensity)))
    )
    log_variances = (
        -1.4
        + trend
        + rng.normal(
            loc=0.0,
            scale=0.25,
            size=int(n_features),
        )
    )
    permutation = rng.permutation(int(n_features))
    return np.exp(log_variances[permutation]).astype(float), mean_intensity[permutation]


def _previous_exact_mean_variance_trend(
    mean_intensity: np.ndarray,
    log_variances: np.ndarray,
    *,
    span: float = 0.4,
) -> np.ndarray:
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
            fitted[idx] = float(np.mean(y[distance == 0.0]))
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
        fitted[idx] = (t0 * s2 - t1 * s1) / determinant

    unsorted = np.empty_like(fitted)
    unsorted[order] = fitted
    return unsorted


def _fit_with_trend(
    *,
    variances: np.ndarray,
    mean_intensity: np.ndarray,
):
    return fit_empirical_bayes(
        variances=variances,
        residual_dof=8.0,
        method="standard",
        trend=True,
        winsor_tail_p=(0.05, 0.10),
        mean_intensity=mean_intensity,
    )


def test_empirical_bayes_trend_smoothing_is_deterministic() -> None:
    variances, mean_intensity = _build_trend_inputs(n_features=3_000)

    first = _fit_with_trend(variances=variances, mean_intensity=mean_intensity)
    second = _fit_with_trend(variances=variances, mean_intensity=mean_intensity)

    np.testing.assert_array_equal(first.prior_variance, second.prior_variance)
    np.testing.assert_array_equal(
        first.prior_degrees_of_freedom,
        second.prior_degrees_of_freedom,
    )
    assert first.fitted_log_prior_variance is not None
    assert second.fitted_log_prior_variance is not None
    np.testing.assert_array_equal(
        first.fitted_log_prior_variance,
        second.fitted_log_prior_variance,
    )


def test_empirical_bayes_trend_smoothing_matches_previous_behavior_within_tolerance() -> (
    None
):
    variances, mean_intensity = _build_trend_inputs(n_features=2_500)
    log_variances = np.log(variances)
    previous_trend = _previous_exact_mean_variance_trend(
        mean_intensity,
        log_variances,
    )
    previous_trend = previous_trend - float(np.mean(previous_trend))
    previous_base_prior_variance, previous_base_prior_dof = fit_f_dist(
        np.exp(log_variances - previous_trend),
        residual_dof=8.0,
    )
    previous_fitted_log_prior_variance = previous_trend + np.log(
        previous_base_prior_variance
    )

    result = _fit_with_trend(variances=variances, mean_intensity=mean_intensity)

    assert result.fitted_log_prior_variance is not None
    # Large feature counts use deterministic anchor interpolation. Keep the
    # tolerance in fitted-log-prior units below the fixture noise scale.
    np.testing.assert_allclose(
        result.fitted_log_prior_variance,
        previous_fitted_log_prior_variance,
        rtol=0.0,
        atol=2.0e-4,
    )
    np.testing.assert_allclose(
        result.base_prior_variance,
        previous_base_prior_variance,
        rtol=1.0e-5,
        atol=0.0,
    )
    np.testing.assert_allclose(
        result.base_prior_degrees_of_freedom,
        previous_base_prior_dof,
        rtol=1.0e-4,
        atol=0.0,
    )
