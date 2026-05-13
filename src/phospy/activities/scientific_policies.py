"""Scientific policy records for activity inference summaries."""

from __future__ import annotations

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.scoring.policy_models import ThresholdMode


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


def _threshold_operator_token(mode: ThresholdMode) -> str:
    if mode is ThresholdMode.GREATER_THAN:
        return ">"
    if mode is ThresholdMode.GREATER_THAN_OR_EQUAL:
        return ">="
    return ">"


__all__ = [
    "build_ksea_zscore_activity_policy",
    "build_simplified_weighted_substrate_activity_policy",
]
