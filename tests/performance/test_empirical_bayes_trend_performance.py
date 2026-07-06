from __future__ import annotations

from typing import cast

import numpy as np
import pytest

from phospy.science.differential.empirical_bayes import (
    EmpiricalBayesFit,
    fit_empirical_bayes,
)
from tests.support.performance_contracts import (
    DEFAULT_PERFORMANCE_SEED,
    EMPIRICAL_BAYES_TREND_LARGE_N_FEATURES,
    EMPIRICAL_BAYES_TREND_LARGE_PEAK_MIB_MAX,
    EMPIRICAL_BAYES_TREND_LARGE_RUNTIME_SECONDS_MAX,
    EMPIRICAL_BAYES_TREND_MEDIUM_N_FEATURES,
    EMPIRICAL_BAYES_TREND_MEDIUM_PEAK_MIB_MAX,
    EMPIRICAL_BAYES_TREND_MEDIUM_RUNTIME_SECONDS_MAX,
    EMPIRICAL_BAYES_TREND_SMALL_N_FEATURES,
    EMPIRICAL_BAYES_TREND_SMALL_PEAK_MIB_MAX,
    EMPIRICAL_BAYES_TREND_SMALL_RUNTIME_SECONDS_MAX,
    measure_runtime_and_peak_mib,
)

pytestmark = [pytest.mark.performance, pytest.mark.release_gate]

_FEATURE_COUNT_CASES = (
    pytest.param(
        "small",
        EMPIRICAL_BAYES_TREND_SMALL_N_FEATURES,
        EMPIRICAL_BAYES_TREND_SMALL_RUNTIME_SECONDS_MAX,
        EMPIRICAL_BAYES_TREND_SMALL_PEAK_MIB_MAX,
        id="small-1000-features",
    ),
    pytest.param(
        "medium",
        EMPIRICAL_BAYES_TREND_MEDIUM_N_FEATURES,
        EMPIRICAL_BAYES_TREND_MEDIUM_RUNTIME_SECONDS_MAX,
        EMPIRICAL_BAYES_TREND_MEDIUM_PEAK_MIB_MAX,
        id="medium-10000-features",
    ),
)


def _build_trend_inputs(*, n_features: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(DEFAULT_PERFORMANCE_SEED + int(n_features))
    mean_intensity = np.linspace(7.5, 13.5, int(n_features), dtype=float)
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
    return np.exp(log_variances).astype(float), mean_intensity


def _assert_trend_fit_performance(
    *,
    n_features: int,
    runtime_seconds_max: float,
    peak_mib_max: float,
) -> None:
    variances, mean_intensity = _build_trend_inputs(n_features=n_features)

    def run_fit() -> EmpiricalBayesFit:
        return fit_empirical_bayes(
            variances=variances,
            residual_dof=8.0,
            method="standard",
            trend=True,
            winsor_tail_p=(0.05, 0.10),
            mean_intensity=mean_intensity,
        )

    measured_result, runtime_seconds, peak_mib = measure_runtime_and_peak_mib(
        run_fit,
        warmup=True,
    )
    result = cast(EmpiricalBayesFit, measured_result)

    assert result.prior_variance.shape == (n_features,)
    assert result.prior_degrees_of_freedom.shape == (n_features,)
    assert result.mean_intensity is not None
    assert result.mean_intensity.shape == (n_features,)
    assert result.log_residual_variance is not None
    assert result.log_residual_variance.shape == (n_features,)
    assert result.fitted_log_prior_variance is not None
    assert result.fitted_log_prior_variance.shape == (n_features,)
    assert np.isfinite(result.prior_variance).all()
    assert np.isfinite(result.fitted_log_prior_variance).all()
    assert runtime_seconds < runtime_seconds_max
    assert peak_mib < peak_mib_max


@pytest.mark.parametrize(
    (
        "_scale_name",
        "n_features",
        "runtime_seconds_max",
        "peak_mib_max",
    ),
    _FEATURE_COUNT_CASES,
)
def test_empirical_bayes_trend_smoothing_performance_realistic_feature_counts(
    _scale_name: str,
    n_features: int,
    runtime_seconds_max: float,
    peak_mib_max: float,
) -> None:
    _assert_trend_fit_performance(
        n_features=n_features,
        runtime_seconds_max=runtime_seconds_max,
        peak_mib_max=peak_mib_max,
    )


def test_empirical_bayes_trend_smoothing_performance_improves_large_matrix() -> None:
    _assert_trend_fit_performance(
        n_features=EMPIRICAL_BAYES_TREND_LARGE_N_FEATURES,
        runtime_seconds_max=EMPIRICAL_BAYES_TREND_LARGE_RUNTIME_SECONDS_MAX,
        peak_mib_max=EMPIRICAL_BAYES_TREND_LARGE_PEAK_MIB_MAX,
    )
