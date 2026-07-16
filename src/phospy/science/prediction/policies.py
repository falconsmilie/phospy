"""Adaptive prediction sampling policies and RNG source resolution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

from phospy.contracts.configs.prediction import (
    KINASE_ADAPTIVE_POLICY_R_PARITY,
    KINASE_ADAPTIVE_POLICY_STABLE,
    KinaseAdaptivePolicy,
)
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)

PredictionSamplingSeedStrategy = Literal["stable_by_kinase", "global_parity"]
PredictionResamplingWeightMode = Literal["default", "r_parity"]
PredictionFinalScoreMode = Literal["mean_probability", "decision_sigmoid"]


@dataclass(frozen=True, slots=True)
class KinaseLibraryMotifScoringPolicy:
    """Metadata policy for pure Kinase Library-style motif scoring."""

    score_scale: str
    residue_classes: tuple[str, ...]
    upstream_residues: int
    downstream_residues: int
    sequence_semantics: str
    reference_distributions_supplied: bool = False
    higher_is_better: bool = True

    @property
    def record(self) -> ScientificPolicyRecord:
        return ScientificPolicyRecord(
            id=ScientificPolicyId.KINASE_LIBRARY_MOTIF_SCORING,
            name="kinase_library_motif_scoring_v1",
            version="1",
            description=(
                "Scores phosphosite sequence windows against Kinase "
                "Library-style position-specific matrices without using known "
                "kinase-substrate edges."
            ),
            parameters={
                "score_scale": str(self.score_scale),
                "residue_classes": "|".join(str(item) for item in self.residue_classes),
                "upstream_residues": int(self.upstream_residues),
                "downstream_residues": int(self.downstream_residues),
                "sequence_semantics": str(self.sequence_semantics),
                "reference_distributions_supplied": bool(
                    self.reference_distributions_supplied
                ),
                "higher_is_better": bool(self.higher_is_better),
                "raw_score_formula": (
                    "sum(matrix[amino_acid_at_relative_position, relative_position])"
                ),
                "percentile_method": (
                    (
                        "100 * count(reference_score <= site_score) / n"
                        if self.higher_is_better
                        else "100 * count(reference_score >= site_score) / n"
                    )
                    if self.reference_distributions_supplied
                    else None
                ),
                "rank_method": (
                    (
                        "1 + count(reference_score > site_score)"
                        if self.higher_is_better
                        else "1 + count(reference_score < site_score)"
                    )
                    if self.reference_distributions_supplied
                    else None
                ),
            },
            assumptions=(
                "Ser/Thr and Tyr residue-class lanes are not interchangeable.",
                "Missing or invalid sequence windows remain unscored rather than "
                "receiving neutral scores.",
                "Percentiles and ranks are empirical summaries of caller-supplied "
                "reference distributions when those distributions are available.",
            ),
            output_scale=(
                "Raw provider-scale motif score sums, with optional empirical "
                "percentile and rank matrices."
            ),
            quantitative_meaning="motif_sequence_match_score",
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
    "DEFAULT_PREDICTION_SAMPLING_POLICY",
    "KinaseLibraryMotifScoringPolicy",
    "PredictionFinalScoreMode",
    "PredictionResamplingWeightMode",
    "PredictionSamplingPolicy",
    "PredictionSamplingSeedStrategy",
    "R_PARITY_PREDICTION_SAMPLING_POLICY",
    "resolve_prediction_sampling_policy",
]
