"""Kinase adaptive ensemble prediction execution helpers."""

from __future__ import annotations

from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.kinase.prediction.config import KinaseAdaptivePolicy
from phospy.kinase.prediction.policies import (
    PredictionSamplingPolicy,
    resolve_prediction_sampling_policy,
)
from phospy.kinase.prediction.sampling_core import run_adaptive_sampling_ensemble
from phospy.kinase.prediction.sampling_runtime import PredictionSamplingRandomSource

_STABLE_PARALLEL_MIN_KINASES = 2
_STABLE_MAX_WORKERS = 2


@dataclass(frozen=True, slots=True)
class AdaptiveEnsemblePredictionExecutionConfig:
    """Resolved adaptive prediction policy consumed by numerical execution."""

    adaptive_policy: KinaseAdaptivePolicy
    adaptive_ensemble_runs: int
    n_iterations: int


@dataclass(frozen=True, slots=True)
class _AdaptiveKinaseTask:
    kinase: str
    feature_values: np.ndarray
    positive_values: np.ndarray
    negative_positions: np.ndarray
    negative_sampling_rng: np.random.Generator
    resampling_rng: np.random.Generator
    prediction_config: AdaptiveEnsemblePredictionExecutionConfig
    sampling_policy: PredictionSamplingPolicy
    kernel: str


def _score_candidate_kinase(task: _AdaptiveKinaseTask) -> tuple[str, np.ndarray]:
    positive_count = int(task.positive_values.shape[0])
    labels = np.concatenate(
        [
            np.repeat(1, positive_count),
            np.repeat(2, positive_count),
        ]
    )
    aggregated_scores = np.zeros(int(task.feature_values.shape[0]), dtype=float)
    for _ in range(task.prediction_config.adaptive_ensemble_runs):
        sampled_negative_positions = task.negative_sampling_rng.choice(
            task.negative_positions,
            size=positive_count,
            replace=int(task.negative_positions.size) < positive_count,
        )
        train_values = np.concatenate(
            [
                task.positive_values,
                task.feature_values[sampled_negative_positions, :],
            ],
            axis=0,
        )
        ensemble_scores = run_adaptive_sampling_ensemble(
            train_values=train_values,
            train_labels=labels,
            test_values=task.feature_values,
            kernel=task.kernel,
            n_iterations=task.prediction_config.n_iterations,
            resampling_rng=task.resampling_rng,
            sampling_policy=task.sampling_policy,
        )
        aggregated_scores += ensemble_scores

    mean_scores = aggregated_scores / float(
        task.prediction_config.adaptive_ensemble_runs
    )
    finite = np.isfinite(mean_scores)
    mean_scores[finite] = np.clip(mean_scores[finite], 0.0, 1.0)
    return task.kinase, mean_scores


def _score_candidate_kinases(
    tasks: list[_AdaptiveKinaseTask],
    *,
    parallel: bool,
) -> dict[str, np.ndarray]:
    if not parallel or len(tasks) < _STABLE_PARALLEL_MIN_KINASES:
        return dict(_score_candidate_kinase(task) for task in tasks)

    max_workers = min(_STABLE_MAX_WORKERS, len(tasks))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        scored = executor.map(_score_candidate_kinase, tasks)
        return dict(scored)


def run_adaptive_ensemble_prediction(
    *,
    prediction_score_matrix: pd.DataFrame,
    candidate_substrates: Mapping[str, list[str]],
    prediction_config: AdaptiveEnsemblePredictionExecutionConfig,
    random_state: int,
    kernel: str = "rbf",
) -> pd.DataFrame:
    """Run adaptive ensemble execution for all candidate kinases."""

    if prediction_score_matrix.empty:
        return pd.DataFrame(index=prediction_score_matrix.index.copy(), dtype=float)

    feature_values = prediction_score_matrix.to_numpy(dtype=float, copy=False)
    feature_index = prediction_score_matrix.index
    all_positions = np.arange(len(feature_index), dtype=int)
    site_positions = {site: i for i, site in enumerate(feature_index)}
    sampling_policy = resolve_prediction_sampling_policy(
        prediction_config.adaptive_policy
    )
    random_source = PredictionSamplingRandomSource(
        policy=sampling_policy,
        random_state=random_state,
    )
    tasks: list[_AdaptiveKinaseTask] = []

    for kinase, substrates in candidate_substrates.items():
        if kinase not in prediction_score_matrix.columns:
            continue
        unique_substrates = list(dict.fromkeys(substrates))
        positive_position_list = [
            site_positions[site] for site in unique_substrates if site in site_positions
        ]
        positive_positions = np.asarray(
            positive_position_list,
            dtype=int,
        )
        if positive_positions.size == 0:
            continue

        negative_mask = np.ones(len(all_positions), dtype=bool)
        negative_mask[positive_positions] = False
        negative_positions = all_positions[negative_mask]
        if negative_positions.size == 0:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "prediction.adaptive_negative_pool; "
                f"kinase={kinase} has no negative pool rows"
            )

        negative_sampling_rng, resampling_rng = random_source.generators_for_kinase(
            kinase=kinase
        )
        tasks.append(
            _AdaptiveKinaseTask(
                kinase=str(kinase),
                feature_values=feature_values,
                positive_values=feature_values[positive_positions, :],
                negative_positions=negative_positions,
                negative_sampling_rng=negative_sampling_rng,
                resampling_rng=resampling_rng,
                prediction_config=prediction_config,
                sampling_policy=sampling_policy,
                kernel=kernel,
            )
        )

    kinase_scores = _score_candidate_kinases(
        tasks,
        parallel=sampling_policy.seed_strategy == "stable_by_kinase",
    )

    return pd.DataFrame(
        kinase_scores,
        index=feature_index.copy(),
        copy=False,
    )


__all__ = [
    "AdaptiveEnsemblePredictionExecutionConfig",
    "run_adaptive_ensemble_prediction",
]
