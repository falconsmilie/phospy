from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..types import PredictionSvmMode

PredictionSamplingSeedStrategy = Literal["stable_by_kinase", "global_parity"]
PredictionResamplingWeightMode = Literal["default", "r_parity"]
PredictionFinalScoreMode = Literal["mean_probability", "decision_sigmoid"]


@dataclass(frozen=True, slots=True)
class PredictionSamplingPolicy:
    """Central prediction-sampling contract for a public prediction mode."""

    name: PredictionSvmMode
    seed_strategy: PredictionSamplingSeedStrategy
    resampling_weight_mode: PredictionResamplingWeightMode
    final_score_mode: PredictionFinalScoreMode


DEFAULT_PREDICTION_SAMPLING_POLICY = PredictionSamplingPolicy(
    name="default",
    seed_strategy="stable_by_kinase",
    resampling_weight_mode="default",
    final_score_mode="mean_probability",
)

R_PARITY_PREDICTION_SAMPLING_POLICY = PredictionSamplingPolicy(
    name="r_parity",
    seed_strategy="global_parity",
    resampling_weight_mode="r_parity",
    final_score_mode="decision_sigmoid",
)


def resolve_prediction_sampling_policy(
    svm_mode: PredictionSvmMode,
) -> PredictionSamplingPolicy:
    if svm_mode == "default":
        return DEFAULT_PREDICTION_SAMPLING_POLICY
    return R_PARITY_PREDICTION_SAMPLING_POLICY


class PredictionSamplingRandomSource:
    """Resolve per-kinase RNG streams using a central sampling policy."""

    def __init__(
        self,
        *,
        policy: PredictionSamplingPolicy,
        random_state: int | None,
    ) -> None:
        self.policy = policy
        self.random_state = random_state
        self._global_rng: np.random.Generator | None = None
        if self.policy.seed_strategy == "global_parity":
            self._global_rng = np.random.default_rng(random_state)

    def generators_for_kinase(
        self,
        *,
        kinase: str,
    ) -> tuple[np.random.Generator, np.random.Generator]:
        if self.policy.seed_strategy == "stable_by_kinase":
            from .sampling_runtime import make_kinase_prediction_random_generators

            return make_kinase_prediction_random_generators(
                random_state=self.random_state,
                kinase=kinase,
            )

        if self._global_rng is None:
            self._global_rng = np.random.default_rng(self.random_state)

        from .sampling_runtime import make_prediction_random_generators

        return make_prediction_random_generators(self._global_rng)


__all__ = [
    "DEFAULT_PREDICTION_SAMPLING_POLICY",
    "PredictionFinalScoreMode",
    "PredictionResamplingWeightMode",
    "PredictionSamplingPolicy",
    "PredictionSamplingRandomSource",
    "PredictionSamplingSeedStrategy",
    "R_PARITY_PREDICTION_SAMPLING_POLICY",
    "resolve_prediction_sampling_policy",
]
