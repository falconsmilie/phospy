#!/usr/bin/env python3
"""Compare fixed-order profile dot products with the former BLAS reduction."""

from __future__ import annotations

import json
import platform
import statistics
import sys
import time

import numpy as np

from phospy.kinase.scoring.profile import _fixed_order_dot_products


def _measure(operation, *, repeats: int) -> tuple[float, np.ndarray]:
    timings: list[float] = []
    result: np.ndarray | None = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = operation()
        timings.append(time.perf_counter() - start)
    if result is None:
        raise AssertionError("benchmark requires at least one repeat")
    return statistics.median(timings), result


def main() -> None:
    site_count = 5_000
    kinase_count = 100
    sample_count = 12
    repeats = 21
    generator = np.random.default_rng(20261018)
    left = generator.normal(size=(site_count, sample_count))
    right = generator.normal(size=(kinase_count, sample_count))
    left -= left.mean(axis=1, keepdims=True)
    right -= right.mean(axis=1, keepdims=True)

    _ = left @ right.T
    _fixed_order_dot_products(left, right)
    blas_seconds, blas_result = _measure(lambda: left @ right.T, repeats=repeats)
    fixed_seconds, fixed_result = _measure(
        lambda: _fixed_order_dot_products(left, right),
        repeats=repeats,
    )
    np.testing.assert_allclose(fixed_result, blas_result, rtol=1e-14, atol=1e-14)
    record = {
        "sites": site_count,
        "kinases": kinase_count,
        "samples": sample_count,
        "repeats": repeats,
        "blas_median_seconds": blas_seconds,
        "fixed_order_median_seconds": fixed_seconds,
        "fixed_over_blas_ratio": fixed_seconds / blas_seconds,
        "maximum_absolute_numerator_difference": float(
            np.max(np.abs(fixed_result - blas_result))
        ),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "executable": sys.executable,
    }
    print("profile_score_portability=" + json.dumps(record, sort_keys=True))


if __name__ == "__main__":
    main()
