#!/usr/bin/env python3
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _time_call(repeats: int, func, *args, **kwargs) -> float:
    durations: list[float] = []
    for _ in range(repeats):
        start = time.perf_counter()
        func(*args, **kwargs)
        durations.append(time.perf_counter() - start)
    return sum(durations) / len(durations)


def _build_benchmark_inputs(
    *,
    n_sites: int = 240,
    n_features: int = 32,
    random_state: int = 11,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    rng = np.random.default_rng(random_state)
    positive = rng.normal(loc=0.8, scale=0.1, size=(n_sites // 2, n_features))
    negative = rng.normal(loc=0.2, scale=0.1, size=(n_sites // 2, n_features))
    train_values = np.vstack([positive, negative]).clip(0.0, 1.0)
    train_mat = pd.DataFrame(
        train_values,
        index=[f"SITE_{index + 1}" for index in range(n_sites)],
        columns=[f"feature_{index + 1}" for index in range(n_features)],
    )
    test_mat = train_mat.copy(deep=True)
    labels = np.asarray([1] * (n_sites // 2) + [2] * (n_sites // 2), dtype=int)
    return train_mat, test_mat, labels


def main() -> None:
    from phospy.prediction.sampling import multi_ada_sampling

    train_mat, test_mat, labels = _build_benchmark_inputs()
    train_values = train_mat.to_numpy(dtype=float)
    test_values = test_mat.to_numpy(dtype=float)
    repeats = 3
    common_kwargs = dict(
        labels=labels,
        kernel="rbf",
        n_iterations=5,
        capture_trace=False,
        trace_level="none",
        trace_sink=None,
        kinase="KINASE_A",
        ensemble_index=1,
        initial_negative_sites=train_mat.index[len(train_mat) // 2 :].tolist(),
        debug_top_n=10,
        svm_mode="default",
        sampling_override=None,
    )

    dataframe_runtime = _time_call(
        repeats,
        multi_ada_sampling,
        train_mat=train_mat,
        test_mat=test_mat,
        resampling_rng=np.random.default_rng(17),
        **common_kwargs,
    )
    array_runtime = _time_call(
        repeats,
        multi_ada_sampling,
        train_mat=None,
        test_mat=None,
        resampling_rng=np.random.default_rng(17),
        train_values=train_values,
        train_index=train_mat.index,
        test_values=test_values,
        test_index=test_mat.index,
        **common_kwargs,
    )

    speedup = dataframe_runtime / array_runtime if array_runtime > 0 else float("inf")
    print(f"dataframe_path_mean_seconds={dataframe_runtime:.6f}")
    print(f"array_path_mean_seconds={array_runtime:.6f}")
    print(f"array_path_speedup={speedup:.3f}x")


if __name__ == "__main__":
    main()
