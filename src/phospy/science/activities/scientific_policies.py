"""Scientific policy records for activity inference summaries."""

from __future__ import annotations

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.science.scoring.policy_models import ThresholdMode


def build_simplified_weighted_substrate_activity_policy(
    *,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
        name="Simplified Weighted Substrate Activity",
        version="1",
        description=(
            "Computes weighted activity and thresholded substrate-mean activity "
            "from predicted substrate support."
        ),
        parameters={
            "threshold": float(threshold),
            "min_substrates": int(min_substrates),
            "top_n_substrates": int(top_n_substrates),
            "weighted_activity_rule": (
                "prediction-weighted mean over top-N predicted substrates"
            ),
            "thresholded_activity_rule": (
                "mean phospho over predicted substrates with "
                f"{ThresholdMode.GREATER_THAN_OR_EQUAL.value}"
            ),
        },
        assumptions=(
            "Predicted substrate support approximates kinase-substrate relevance.",
            "Higher weighted/thresholded values indicate stronger relative activity "
            "support in-run.",
            "This is not full KSEA-style statistical enrichment.",
        ),
        output_scale=(
            "Sample-by-kinase relative activity summaries (weighted mean and "
            "thresholded mean)."
        ),
        quantitative_meaning="relative_activity_support",
    )


def build_ksea_zscore_activity_policy(
    *,
    evidence_threshold: float,
    min_substrates: int,
    p_value_method: str,
    adjust_p_values: bool,
    q_value_method: str | None,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.KSEA_ZSCORE_ACTIVITY,
        name="ksea_zscore_activity_v1",
        version="1",
        description=(
            "Computes KSEA-style substrate-set enrichment activity z-scores using "
            "unweighted substrate membership after evidence thresholding."
        ),
        parameters={
            "evidence_threshold": float(evidence_threshold),
            "min_substrates": int(min_substrates),
            "membership_rule": (
                "finite_evidence "
                f"{_threshold_operator_token(ThresholdMode.GREATER_THAN_OR_EQUAL)} "
                "evidence_threshold"
            ),
            "weighting_rule": "unweighted_membership",
            "z_score_formula": "(mean_S - mean_U) * sqrt(n) / sd_U",
            "background_sd_ddof": 1,
            "p_value_method": str(p_value_method),
            "adjust_p_values": bool(adjust_p_values),
            "q_value_method": None if q_value_method is None else str(q_value_method),
        },
        assumptions=(
            "Substrate evidence contributes as binary membership after thresholding.",
            "Background phosphosite values define per-condition mean and sample variance.",
            "Scores with insufficient substrates or invalid background variance are not computable.",
            "KSEA z-scores are statistical enrichment summaries and are not PhosR-equivalent activity inference.",
        ),
        output_scale=(
            "Condition-by-kinase z-score substrate-set enrichment activity matrix "
            "with normal-approximation p-values."
        ),
        quantitative_meaning="substrate_set_enrichment_z_score",
    )


def build_ssgsea_substrate_enrichment_activity_policy(
    *,
    min_substrates: int,
    ranking_direction: str,
    permutation_count: int,
    random_seed: int | None,
    adjust_p_values: bool,
    q_value_method: str | None,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY,
        name="ssgsea_substrate_enrichment_activity_v1",
        version="1",
        description=(
            "Computes a PhosPy ssGSEA-style kinase substrate-set enrichment "
            "score over ranked phosphosite effect values."
        ),
        parameters={
            "min_substrates": int(min_substrates),
            "ranking_direction": str(ranking_direction),
            "rank_walk_rule": (
                "stable rank order; hits increment by 1/n_substrates and misses "
                "decrement by 1/n_non_substrates"
            ),
            "score_formula": "sum(cumulative_hit - cumulative_miss) / n_background",
            "membership_rule": "explicit kinase-substrate membership table",
            "permutation_count": int(permutation_count),
            "p_value_method": (
                "seeded site-label permutation two-sided empirical p-value"
                if int(permutation_count) > 0
                else None
            ),
            "random_seed": None if random_seed is None else int(random_seed),
            "adjust_p_values": bool(adjust_p_values),
            "q_value_method": None if q_value_method is None else str(q_value_method),
        },
        assumptions=(
            "Kinase substrate membership defines the tested phosphosite set.",
            "Rank concentration of substrate effects summarizes relative kinase activity support.",
            "Permutation p-values, when requested, use seeded random substrate-set label permutations.",
            "This is a validated PhosPy implementation and is not a PTM-SEA parity claim.",
        ),
        output_scale=(
            "Condition-by-kinase rank-walk substrate-set enrichment activity "
            "score matrix with optional empirical p-values."
        ),
        quantitative_meaning="rank_based_substrate_set_enrichment_score",
    )


def _threshold_operator_token(mode: ThresholdMode) -> str:
    if mode is ThresholdMode.GREATER_THAN:
        return ">"
    if mode is ThresholdMode.GREATER_THAN_OR_EQUAL:
        return ">="
    return ">"


__all__ = [
    "build_ksea_zscore_activity_policy",
    "build_simplified_weighted_substrate_activity_policy",
    "build_ssgsea_substrate_enrichment_activity_policy",
]
