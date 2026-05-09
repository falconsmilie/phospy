"""Limma-style empirical-Bayes helpers."""

from __future__ import annotations

import numpy as np
from scipy.special import digamma, polygamma


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
