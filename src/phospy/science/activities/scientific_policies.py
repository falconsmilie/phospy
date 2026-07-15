"""Scientific policy records for exploratory activity-like score summaries."""

from __future__ import annotations

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.science.scoring.policy_models import ThresholdMode

SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION = "1"
SSGSEA_PERMUTATION_RNG_SEED_POLICY = "stable_by_method_condition_kinase"
SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION = "1"
SSGSEA_PERMUTATION_RNG_SEED_MATERIAL = (
    "blake2b-128-json(method_id, method_version, seed_policy, "
    "seed_policy_version, random_seed, condition, kinase, stream)"
)


def build_simplified_weighted_substrate_activity_policy(
    *,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
        name="Simplified Weighted Substrate Activity Score",
        version="1",
        description=(
            "Computes a substrate-supported kinase activity score and "
            "thresholded substrate-mean activity-like summary from predicted "
            "substrate support."
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
            "Predicted substrate support approximates kinase-substrate relevance "
            "for exploratory scoring.",
            "Higher weighted/thresholded values indicate stronger relative "
            "candidate kinase support in-run.",
            "Sparse or missing substrate support weakens interpretation.",
            "The score does not prove kinase activation or causal regulation; "
            "causal kinase activity claims require external validation.",
            "This is not full KSEA-style statistical enrichment.",
        ),
        output_scale=(
            "Sample-by-kinase exploratory activity-like summaries (weighted mean "
            "and thresholded mean)."
        ),
        quantitative_meaning="relative_substrate_supported_kinase_score",
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
            "Computes KSEA-style inferred kinase activity z-scores using "
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
            "Background phosphosite values define per-condition mean and sample "
            "variance.",
            "Scores with insufficient substrates or invalid background variance "
            "are not computable.",
            "Sparse or missing substrate support weakens interpretation.",
            "KSEA z-scores are statistical enrichment summaries, not validated "
            "causal kinase activation.",
            "Causal kinase activity claims require external validation and a "
            "study design that supports them.",
            "KSEA z-scores are not PhosR-equivalent activity inference.",
        ),
        output_scale=(
            "Condition-by-kinase inferred kinase activity score matrix "
            "(z-score substrate-set enrichment) with normal-approximation p-values."
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
    has_permutations = int(permutation_count) > 0
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY,
        name="ssgsea_substrate_enrichment_activity_v1",
        version=SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION,
        description=(
            "Computes a PhosPy ssGSEA-style kinase substrate-set enrichment "
            "score over ranked phosphosite effect values."
        ),
        parameters={
            "method_id": ScientificPolicyId.SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY.value,
            "method_version": SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION,
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
            "permutation_rng_seed_policy": (
                SSGSEA_PERMUTATION_RNG_SEED_POLICY if has_permutations else None
            ),
            "permutation_rng_seed_policy_version": (
                SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION if has_permutations else None
            ),
            "permutation_rng_seed_material": (
                SSGSEA_PERMUTATION_RNG_SEED_MATERIAL if has_permutations else None
            ),
            "adjust_p_values": bool(adjust_p_values),
            "q_value_method": None if q_value_method is None else str(q_value_method),
        },
        assumptions=(
            "Kinase substrate membership defines the tested phosphosite set.",
            "Rank concentration of substrate effects summarizes candidate kinase "
            "support.",
            "Sparse or missing substrate support weakens interpretation.",
            "Permutation p-values, when requested, use seeded random substrate-set "
            "label permutations with deterministic child RNG streams keyed by "
            "method, condition, kinase, and user seed.",
            "The enrichment score does not prove kinase activation or causal "
            "regulation; causal kinase activity claims require external validation.",
            "This is a validated PhosPy implementation and is not a PTM-SEA "
            "parity claim.",
        ),
        output_scale=(
            "Condition-by-kinase rank-walk substrate-supported kinase score matrix "
            "with optional empirical p-values."
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
    "SSGSEA_PERMUTATION_RNG_SEED_MATERIAL",
    "SSGSEA_PERMUTATION_RNG_SEED_POLICY",
    "SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION",
    "SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION",
]
