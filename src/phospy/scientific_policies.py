"""Stable scientific-policy identifiers and serializable metadata records."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np

ScientificPolicyParameter = str | int | float | bool | None


class ScientificPolicyId(str, Enum):
    """Stable identifiers for scientific scoring and derivation behavior."""

    PROFILE_CORRELATION_SHIFTED_UNIT = "profile_correlation_shifted_unit_v1"
    MOTIF_PROFILE_RANK_FUSION = "motif_profile_rank_fusion_v1"
    SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY = "simplified_weighted_substrate_activity_v1"
    SIGNALOME_MODULE_CANDIDATE_SCORE = "signalome_module_candidate_score_v1"
    PROTEIN_MODULE_FROM_SITE_MEMBERSHIP = "protein_module_from_site_membership_v1"


@dataclass(frozen=True, slots=True)
class ScientificPolicyRecord:
    """Serializable metadata for one scientific scoring/derivation policy."""

    id: ScientificPolicyId
    name: str
    version: str
    description: str
    parameters: dict[str, ScientificPolicyParameter]
    assumptions: tuple[str, ...]
    output_scale: str | None = None
    quantitative_meaning: str | None = None

    def to_payload(self) -> dict[str, object]:
        return {
            "id": self.id.value,
            "name": self.name,
            "version": self.version,
            "description": self.description,
            "parameters": dict(self.parameters),
            "assumptions": list(self.assumptions),
            "output_scale": self.output_scale,
            "quantitative_meaning": self.quantitative_meaning,
        }

    @classmethod
    def from_payload(cls, payload: dict[str, object]) -> ScientificPolicyRecord:
        assumptions = payload.get("assumptions", ())
        if isinstance(assumptions, (list, tuple)):
            resolved_assumptions = tuple(str(value) for value in assumptions)
        else:
            resolved_assumptions = ()
        parameters_raw = payload.get("parameters", {})
        parameters: dict[str, ScientificPolicyParameter]
        if isinstance(parameters_raw, dict):
            parameters = {}
            for key, value in parameters_raw.items():
                if value is None or isinstance(value, (str, int, float, bool)):
                    parameters[str(key)] = value
                else:
                    parameters[str(key)] = str(value)
        else:
            parameters = {}
        output_scale = payload.get("output_scale")
        resolved_output_scale = None if output_scale is None else str(output_scale)
        quantitative_meaning = payload.get("quantitative_meaning")
        resolved_quantitative_meaning = (
            None if quantitative_meaning is None else str(quantitative_meaning)
        )
        return cls(
            id=ScientificPolicyId(str(payload.get("id"))),
            name=str(payload.get("name")),
            version=str(payload.get("version")),
            description=str(payload.get("description")),
            parameters=parameters,
            assumptions=resolved_assumptions,
            output_scale=resolved_output_scale,
            quantitative_meaning=resolved_quantitative_meaning,
        )


PROFILE_CORRELATION_SHIFTED_UNIT_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.PROFILE_CORRELATION_SHIFTED_UNIT,
    name="Profile Correlation Shifted Unit Support",
    version="1",
    description=(
        "Transforms Pearson correlation from [-1, 1] to [0, 1] using (r + 1) / 2."
    ),
    parameters={
        "transform": "(r + 1) / 2",
        "clip_to_unit_interval": True,
        "preserve_undefined_as_nan": True,
    },
    assumptions=(
        "Higher positive correlation indicates stronger support.",
        "Negative correlation is treated as lower support, not explicit "
        "inhibitory evidence.",
        "Undefined correlations remain missing (NaN).",
    ),
    output_scale=(
        "Relative support score in [0, 1] where larger means stronger positive "
        "profile agreement."
    ),
    quantitative_meaning="relative_support_score",
)


def build_motif_profile_rank_fusion_policy(
    *,
    allow_profile_only_fallback: bool,
    emit_weights: bool,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.MOTIF_PROFILE_RANK_FUSION,
        name="Motif/Profile Rank-Weighted Fusion",
        version="1",
        description=(
            "Combines motif-frequency and profile-correlation scores using "
            "rank-derived logarithmic weights."
        ),
        parameters={
            "motif_weight_formula": "log(rank(motif_size)+1) / total_weight",
            "profile_weight_formula": "log(rank(profile_size)+1) / total_weight",
            "allow_profile_only_fallback": bool(allow_profile_only_fallback),
            "emit_weights": bool(emit_weights),
        },
        assumptions=(
            "Motif-library size and quantified-substrate count proxy evidence "
            "strength.",
            "When motif evidence is missing for a kinase/site, profile evidence can "
            "be propagated.",
            "Outputs are relative support scores and are not calibrated probabilities.",
        ),
        output_scale="Relative downstream support score for kinase-site ranking.",
        quantitative_meaning="relative_support_score",
    )


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
                "mean phospho over predicted substrates with score > threshold"
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


def build_signalome_module_candidate_score_policy(
    *,
    requested_policy: str,
    candidate_scoring_policy: str,
    candidate_scoring_mode: str,
    max_exact_tree_sites: int | None,
    max_full_candidate_scoring_sites: int,
    candidate_scoring_evaluated: bool,
    candidate_scoring_skip_reason: str | None,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.SIGNALOME_MODULE_CANDIDATE_SCORE,
        name="Signalome Module Candidate Score",
        version="1",
        description=(
            "Ranks candidate module counts using within-cluster median "
            "correlation summaries over downstream kinase-score profiles."
        ),
        parameters={
            "requested_policy": str(requested_policy),
            "candidate_scoring_policy": str(candidate_scoring_policy),
            "candidate_scoring_mode": str(candidate_scoring_mode),
            "max_exact_tree_sites": max_exact_tree_sites,
            "max_full_candidate_scoring_sites": int(max_full_candidate_scoring_sites),
            "candidate_scoring_evaluated": bool(candidate_scoring_evaluated),
            "candidate_scoring_skip_reason": (
                None
                if candidate_scoring_skip_reason is None
                else str(candidate_scoring_skip_reason)
            ),
        },
        assumptions=(
            "Candidate quality is summarized by within-cluster correlation coherence.",
            "Degenerate or undefined profiles are excluded or tracked via diagnostics.",
            "Selected module count depends on thresholds and candidate-scoring policy.",
        ),
        output_scale=(
            "Candidate module-count support scores; higher values indicate stronger "
            "within-cluster profile coherence."
        ),
        quantitative_meaning="relative_module_candidate_support",
    )


PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.PROTEIN_MODULE_FROM_SITE_MEMBERSHIP,
    name="Protein Module From Site-Cluster Membership",
    version="1",
    description=(
        "Derives protein-level module IDs by grouping proteins with matching "
        "site-cluster membership patterns."
    ),
    parameters={
        "membership_representation": "binary site-cluster incidence vector",
        "module_id_assignment": "first-seen pattern order",
    },
    assumptions=(
        "Site-cluster membership captures protein-level signaling context.",
        "Proteins with identical site-cluster incidence vectors are grouped into "
        "the same module.",
    ),
    output_scale="Integer module IDs at the protein level.",
    quantitative_meaning="protein_module_membership_label",
)


def shift_correlation_to_unit_support(correlation: np.ndarray) -> np.ndarray:
    """Apply the shifted-unit profile support transform."""

    scores = (correlation + 1.0) / 2.0
    valid = np.isfinite(scores)
    scores[valid] = np.clip(scores[valid], 0.0, 1.0)
    return scores


__all__ = [
    "PROFILE_CORRELATION_SHIFTED_UNIT_POLICY",
    "PROTEIN_MODULE_FROM_SITE_MEMBERSHIP_POLICY",
    "ScientificPolicyId",
    "ScientificPolicyRecord",
    "build_motif_profile_rank_fusion_policy",
    "build_signalome_module_candidate_score_policy",
    "build_simplified_weighted_substrate_activity_policy",
    "shift_correlation_to_unit_support",
]
