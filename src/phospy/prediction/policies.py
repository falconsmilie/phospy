"""Adaptive prediction sampling policies and RNG source resolution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

import numpy as np

from phospy.api.configs import (
    KINASE_ADAPTIVE_POLICY_R_PARITY,
    KINASE_ADAPTIVE_POLICY_STABLE,
    KinaseAdaptivePolicy,
)
from phospy.scientific_policies import ScientificPolicyId, ScientificPolicyRecord

PredictionSamplingSeedStrategy = Literal["stable_by_kinase", "global_parity"]
PredictionResamplingWeightMode = Literal["default", "r_parity"]
PredictionFinalScoreMode = Literal["mean_probability", "decision_sigmoid"]


@dataclass(frozen=True, slots=True)
class PredictionSamplingPolicy:
    """Resolved adaptive-sampling policy contract for one public mode."""

    name: str
    version: str
    parameters: Mapping[str, object]
    description: str
    adaptive_policy: KinaseAdaptivePolicy
    seed_strategy: PredictionSamplingSeedStrategy
    resampling_weight_mode: PredictionResamplingWeightMode
    final_score_mode: PredictionFinalScoreMode

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "parameters",
            MappingProxyType(
                {str(key): value for key, value in self.parameters.items()}
            ),
        )

    @property
    def record(self) -> ScientificPolicyRecord:
        return ScientificPolicyRecord(
            id=ScientificPolicyId.ADAPTIVE_PREDICTION_SAMPLING,
            name=self.name,
            version=self.version,
            description=self.description,
            parameters=self.parameters,
            assumptions=(
                "Adaptive sampling policy changes stochastic training-set draws and "
                "can alter kinase ranking outputs.",
            ),
            output_scale="Adaptive prediction score matrix in [0, 1].",
            quantitative_meaning="relative_prediction_support",
        )


DEFAULT_PREDICTION_SAMPLING_POLICY = PredictionSamplingPolicy(
    name="adaptive_prediction_sampling_stable_v1",
    version="1",
    parameters={
        "adaptive_policy": KINASE_ADAPTIVE_POLICY_STABLE,
        "seed_strategy": "stable_by_kinase",
        "resampling_weight_mode": "default",
        "final_score_mode": "mean_probability",
    },
    description=(
        "Deterministic adaptive sampling with per-kinase seeded RNG streams and "
        "default resampling-weight flattening."
    ),
    adaptive_policy=KINASE_ADAPTIVE_POLICY_STABLE,
    seed_strategy="stable_by_kinase",
    resampling_weight_mode="default",
    final_score_mode="mean_probability",
)

R_PARITY_PREDICTION_SAMPLING_POLICY = PredictionSamplingPolicy(
    name="adaptive_prediction_sampling_r_parity_v1",
    version="1",
    parameters={
        "adaptive_policy": KINASE_ADAPTIVE_POLICY_R_PARITY,
        "seed_strategy": "global_parity",
        "resampling_weight_mode": "r_parity",
        "final_score_mode": "decision_sigmoid",
    },
    description=(
        "Adaptive sampling mode aligned to R-parity behavior with global RNG stream "
        "ordering and parity-compatible score handling."
    ),
    adaptive_policy=KINASE_ADAPTIVE_POLICY_R_PARITY,
    seed_strategy="global_parity",
    resampling_weight_mode="r_parity",
    final_score_mode="decision_sigmoid",
)


def resolve_prediction_sampling_policy(
    adaptive_policy: KinaseAdaptivePolicy,
) -> PredictionSamplingPolicy:
    """Resolve adaptive-sampling policy for the configured public mode."""

    if adaptive_policy == KINASE_ADAPTIVE_POLICY_R_PARITY:
        return R_PARITY_PREDICTION_SAMPLING_POLICY
    return DEFAULT_PREDICTION_SAMPLING_POLICY


class PredictionSamplingRandomSource:
    """Resolve per-kinase RNG streams from the selected sampling policy."""

    def __init__(
        self,
        *,
        policy: PredictionSamplingPolicy,
        random_state: int,
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
            from phospy.prediction.sampling_runtime import (
                make_kinase_prediction_random_generators,
            )

            return make_kinase_prediction_random_generators(
                random_state=self.random_state,
                kinase=kinase,
            )

        if self._global_rng is None:
            self._global_rng = np.random.default_rng(self.random_state)

        from phospy.prediction.sampling_runtime import make_prediction_random_generators

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
