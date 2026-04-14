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


def _build_inputs(
    *,
    n_sites: int = 3000,
    n_kinases: int = 200,
    random_state: int = 31,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    rng = np.random.default_rng(random_state)
    site_index = pd.Index([f"SITE_{i + 1}" for i in range(n_sites)], dtype=object)
    kinase_names = [f"KINASE_{i + 1}" for i in range(n_kinases)]
    feature_mat = pd.DataFrame(
        rng.uniform(0.0, 1.0, size=(n_sites, 4)),
        index=site_index,
        columns=["f1", "f2", "f3", "f4"],
    )
    score_map = {
        kinase: rng.uniform(0.0, 20.0, size=n_sites).astype(float)
        for kinase in kinase_names
    }
    return feature_mat, score_map


def _legacy_accumulate(
    *,
    feature_mat: pd.DataFrame,
    score_map: dict[str, np.ndarray],
) -> pd.DataFrame:
    pred_matrix = pd.DataFrame(
        0.0,
        index=feature_mat.index.copy(),
        columns=list(score_map),
        dtype=float,
    )
    for kinase, score_values in score_map.items():
        batch_scores = pd.Series(
            score_values,
            index=feature_mat.index.copy(),
            dtype=float,
        )
        pred_matrix.loc[:, kinase] += batch_scores.to_numpy(dtype=float, copy=False)
    return pred_matrix


def _array_accumulate(
    *,
    feature_mat: pd.DataFrame,
    score_map: dict[str, np.ndarray],
) -> pd.DataFrame:
    from phospy.prediction.aggregation import PredictionAggregator
    from phospy.prediction.execution import KinasePredictionBatch

    aggregator = PredictionAggregator()
    pred_matrix = aggregator.initialize_prediction_matrix(
        feature_mat=feature_mat,
        substrate_list={kinase: [] for kinase in score_map},
    )
    for kinase, score_values in score_map.items():
        aggregator.add_kinase_scores(
            pred_matrix=pred_matrix,
            batch=KinasePredictionBatch(
                kinase=kinase,
                score_values=score_values,
                score_index=feature_mat.index,
            ),
        )
    return pd.DataFrame(
        pred_matrix.values,
        index=pred_matrix.index,
        columns=pred_matrix.columns,
        copy=False,
    )


def main() -> None:
    feature_mat, score_map = _build_inputs()
    repeats = 5

    baseline = _legacy_accumulate(feature_mat=feature_mat, score_map=score_map)
    optimized = _array_accumulate(feature_mat=feature_mat, score_map=score_map)
    pd.testing.assert_frame_equal(baseline, optimized, check_dtype=False)

    legacy_runtime = _time_call(
        repeats,
        _legacy_accumulate,
        feature_mat=feature_mat,
        score_map=score_map,
    )
    array_runtime = _time_call(
        repeats,
        _array_accumulate,
        feature_mat=feature_mat,
        score_map=score_map,
    )
    speedup = legacy_runtime / array_runtime if array_runtime > 0 else float("inf")
    print(f"legacy_series_dataframe_mean_seconds={legacy_runtime:.6f}")
    print(f"array_buffer_mean_seconds={array_runtime:.6f}")
    print(f"array_buffer_speedup={speedup:.3f}x")


if __name__ == "__main__":
    main()
