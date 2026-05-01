"""Adaptive ensemble prediction execution helpers."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd

from phospy.api.configs import KinasePredictionConfig
from phospy.errors.workflows import WorkflowStageError
from phospy.prediction.policies import (
    PredictionSamplingRandomSource,
    resolve_prediction_sampling_policy,
)
from phospy.prediction.sampling_core import run_adaptive_sampling_ensemble


def run_adaptive_ensemble_prediction(
    *,
    prediction_score_matrix: pd.DataFrame,
    candidate_substrates: Mapping[str, list[str]],
    prediction_config: KinasePredictionConfig,
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
    kinase_scores: dict[str, np.ndarray] = {}

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

        positive_values = feature_values[positive_positions, :]
        labels = np.concatenate(
            [
                np.repeat(1, len(positive_positions)),
                np.repeat(2, len(positive_positions)),
            ]
        )
        negative_sampling_rng, resampling_rng = random_source.generators_for_kinase(
            kinase=kinase
        )
        aggregated_scores = np.zeros(len(all_positions), dtype=float)
        for _ in range(prediction_config.adaptive_ensemble_runs):
            sampled_negative_positions = negative_sampling_rng.choice(
                negative_positions,
                size=len(positive_positions),
                replace=len(negative_positions) < len(positive_positions),
            )
            train_values = np.concatenate(
                [
                    positive_values,
                    feature_values[sampled_negative_positions, :],
                ],
                axis=0,
            )
            ensemble_scores = run_adaptive_sampling_ensemble(
                train_values=train_values,
                train_labels=labels,
                test_values=feature_values,
                kernel=kernel,
                n_iterations=prediction_config.n_iterations,
                resampling_rng=resampling_rng,
                sampling_policy=sampling_policy,
            )
            aggregated_scores += ensemble_scores

        mean_scores = aggregated_scores / float(
            prediction_config.adaptive_ensemble_runs
        )
        finite = np.isfinite(mean_scores)
        mean_scores[finite] = np.clip(mean_scores[finite], 0.0, 1.0)
        kinase_scores[str(kinase)] = mean_scores

    return pd.DataFrame(
        kinase_scores,
        index=feature_index.copy(),
        copy=False,
    )


__all__ = ["run_adaptive_ensemble_prediction"]
