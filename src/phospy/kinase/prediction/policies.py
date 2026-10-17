"""Kinase prediction candidate-selection and adaptive-sampling policies."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

from phospy.kinase.policies import ThresholdMode
from phospy.kinase.prediction.config import (
    KINASE_ADAPTIVE_POLICY_R_PARITY,
    KINASE_ADAPTIVE_POLICY_STABLE,
    KinaseAdaptivePolicy,
)
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)

CANDIDATE_SCORE_THRESHOLD = 0.0
CANDIDATE_MIN_INCLUSION = 1

PredictionSamplingSeedStrategy = Literal["stable_by_kinase", "global_parity"]
PredictionResamplingWeightMode = Literal["default", "r_parity"]
PredictionFinalScoreMode = Literal["mean_probability", "decision_sigmoid"]


@dataclass(frozen=True, slots=True)
class CandidateSubstrateSelectionPolicy:
    """Executable metadata policy for candidate substrate selection behavior."""

    top_k: int
    score_threshold: float
    inclusion: int
    threshold_operator: ThresholdMode = ThresholdMode.GREATER_THAN
    ranking_rule: str = "top_n_scores_per_kinase_then_threshold"
    site_restriction: str = "none"

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "threshold_operator",
            ThresholdMode.parse(
                self.threshold_operator,
                field_name=("candidate substrate selection policy threshold_operator"),
            ),
        )

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_candidate_substrate_selection_policy(
            top_k=self.top_k,
            score_threshold=self.score_threshold,
            inclusion=self.inclusion,
            threshold_operator=self.threshold_operator,
            ranking_rule=self.ranking_rule,
            site_restriction=self.site_restriction,
        )


def build_candidate_substrate_selection_policy(
    *,
    top_k: int,
    score_threshold: float,
    inclusion: int,
    threshold_operator: ThresholdMode | str = ThresholdMode.GREATER_THAN,
    ranking_rule: str = "top_n_scores_per_kinase_then_threshold",
    site_restriction: str = "none",
) -> ScientificPolicyRecord:
    resolved_threshold_mode = ThresholdMode.parse(
        threshold_operator,
        field_name="candidate substrate selection policy threshold_operator",
    )
    return ScientificPolicyRecord(
        id=ScientificPolicyId.CANDIDATE_SUBSTRATE_SELECTION,
        name="Candidate Substrate Selection Policy",
        version="1",
        description=(
            "Selects per-kinase candidate substrate sites from downstream support "
            "scores using top-k ranking, threshold filtering, and minimum inclusion."
        ),
        parameters={
            "top_k": int(top_k),
            "score_threshold": float(score_threshold),
            "inclusion": int(inclusion),
            "threshold_operator": resolved_threshold_mode.value,
            "ranking_rule": str(ranking_rule),
            "site_restriction": str(site_restriction),
        },
        assumptions=(
            "Only finite scores are eligible for candidate selection.",
            "Threshold and inclusion rules jointly determine which kinases are "
            "considered to have usable candidate substrate support.",
            "Changing selection thresholds changes downstream kinase ranking and "
            "prediction outputs.",
        ),
        output_scale=(
            "Per-kinase candidate substrate sets for downstream ranking/prediction."
        ),
        quantitative_meaning="candidate_support_set",
    )


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


__all__ = [
    "CANDIDATE_MIN_INCLUSION",
    "CANDIDATE_SCORE_THRESHOLD",
    "CandidateSubstrateSelectionPolicy",
    "DEFAULT_PREDICTION_SAMPLING_POLICY",
    "PredictionFinalScoreMode",
    "PredictionResamplingWeightMode",
    "PredictionSamplingPolicy",
    "PredictionSamplingSeedStrategy",
    "R_PARITY_PREDICTION_SAMPLING_POLICY",
    "build_candidate_substrate_selection_policy",
    "resolve_prediction_sampling_policy",
]
