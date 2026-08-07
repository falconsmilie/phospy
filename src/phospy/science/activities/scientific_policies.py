"""Scientific policy records for exploratory activity-like score summaries."""

from __future__ import annotations

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.science.activities.membership import (
    ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION,
    ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
    KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION,
)
from phospy.science.scoring.policy_models import ThresholdMode

SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION = "2"
KSEA_ZSCORE_ACTIVITY_POLICY_VERSION = "5"
SSGSEA_PERMUTATION_RNG_SEED_POLICY = "stable_by_method_profile_kinase"
SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION = "1"
SSGSEA_PERMUTATION_RNG_SEED_MATERIAL = (
    "blake2b-128-json(method_id, method_version, profile_id, kinase, "
    "stream, random_seed; v1 compatibility salt retained internally)"
)
SSGSEA_TIE_POLICY = "midrank_block_expectation"
SSGSEA_TIE_BLOCK_CONTRIBUTION_RULE = (
    "Equal finite values form one tie block. A block with h substrates and m "
    "non-substrates contributes the expected rank-walk area over all "
    "within-block orders: b * running_before + ((b + 1) / 2) * "
    "(h / n_substrates - m / n_non_substrates), where b = h + m. The walk then "
    "advances by that same block delta."
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
            "Profile-by-kinase exploratory activity-like summaries (weighted mean "
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
    membership_inferential_eligible: bool | None = None,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.KSEA_ZSCORE_ACTIVITY,
        name="ksea_zscore_activity_v1",
        version=KSEA_ZSCORE_ACTIVITY_POLICY_VERSION,
        description=(
            "Computes KSEA-style inferred kinase activity z-scores using "
            "unweighted substrate membership after evidence thresholding. "
            "Ordinary normal-approximation p-values and BH q-values are emitted "
            "only when typed membership-selection provenance declares the "
            "membership independent of the tested quantitative matrix."
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
            "membership_selection_policy_version": (
                "activity_membership_selection_v"
                f"{ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION}"
            ),
            "membership_selection_payload_schema_version": (
                ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
            ),
            "ksea_membership_inferential_policy_version": (
                KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION
            ),
            "ordinary_p_q_requires_inferentially_eligible_membership": True,
            "ordinary_p_q_availability_derivation": (
                "science_domain_membership_policy_from_typed_provenance"
            ),
            "adaptive_membership_p_q_policy": ("unavailable_descriptive_z_scores_only"),
            "membership_inferential_eligible": membership_inferential_eligible,
        },
        assumptions=(
            "Substrate evidence contributes as binary membership after thresholding.",
            "Normal-approximation p-values assume substrate membership was fixed "
            "independently of the tested quantitative values.",
            "When membership was selected using the tested quantitative matrix, "
            "the method reports descriptive z-scores and substrate counts only; "
            "ordinary p-values and q-values are unavailable unless a valid nested "
            "resampling or sample-splitting procedure is implemented.",
            "Background phosphosite values define per-profile mean and sample "
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
            "Profile-by-kinase inferred kinase activity score matrix "
            "(z-score substrate-set enrichment); normal-approximation p-values "
            "are present only for inferentially eligible membership."
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
            "score over ranked phosphosite effect values with explicit "
            "equal-value tie-block handling."
        ),
        parameters={
            "method_id": ScientificPolicyId.SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY.value,
            "method_version": SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION,
            "min_substrates": int(min_substrates),
            "ranking_direction": str(ranking_direction),
            "tie_policy": SSGSEA_TIE_POLICY,
            "tie_block_contribution_rule": SSGSEA_TIE_BLOCK_CONTRIBUTION_RULE,
            "rank_walk_rule": (
                "untied positions are walked in rank order; equal-valued positions "
                "are walked as tie blocks with equivalent treatment for all sites "
                "inside the block"
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
            "Equal-valued finite sites are not ordered by row position or lexical "
            "site label; mixed substrate/non-substrate tie blocks contribute by "
            "the explicit block expectation policy.",
            "Sparse or missing substrate support weakens interpretation.",
            "Permutation p-values, when requested, use seeded random substrate-set "
            "label permutations scored with the same tie-block policy and "
            "deterministic child RNG streams keyed by method, method version, "
            "profile ID, kinase, stream, and caller-supplied seed. Version 2 "
            "changes the seeded stream identity because the method version is "
            "part of the seed material.",
            "The enrichment score does not prove kinase activation or causal "
            "regulation; causal kinase activity claims require external validation.",
            "This is a validated PhosPy implementation and is not a PTM-SEA "
            "parity claim.",
        ),
        output_scale=(
            "Profile-by-kinase rank-walk substrate-supported kinase score matrix "
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
    "KSEA_ZSCORE_ACTIVITY_POLICY_VERSION",
    "SSGSEA_PERMUTATION_RNG_SEED_MATERIAL",
    "SSGSEA_PERMUTATION_RNG_SEED_POLICY",
    "SSGSEA_PERMUTATION_RNG_SEED_POLICY_VERSION",
    "SSGSEA_SUBSTRATE_ENRICHMENT_ACTIVITY_POLICY_VERSION",
    "SSGSEA_TIE_BLOCK_CONTRIBUTION_RULE",
    "SSGSEA_TIE_POLICY",
]
