"""Scientific policy records for prediction scoring and candidate selection."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.science.scoring.policy_models import (
    ProfileSelfInclusionPolicy,
    ThresholdMode,
)


@dataclass(frozen=True, slots=True)
class KinaseProfileScoringPolicy:
    """Executable metadata policy for kinase profile scoring behavior."""

    profile_missing_value_strategy: str
    min_substrates_floor: int
    requested_min_substrates: int
    profile_self_inclusion_policy: ProfileSelfInclusionPolicy | str = (
        ProfileSelfInclusionPolicy.ALLOW
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "profile_self_inclusion_policy",
            ProfileSelfInclusionPolicy.parse(
                self.profile_self_inclusion_policy,
                field_name="kinase profile scoring policy profile_self_inclusion_policy",
            ),
        )

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_kinase_profile_scoring_policy(
            profile_missing_value_strategy=self.profile_missing_value_strategy,
            min_substrates_floor=self.min_substrates_floor,
            requested_min_substrates=self.requested_min_substrates,
            profile_self_inclusion_policy=self.profile_self_inclusion_policy,
        )


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


PROFILE_CORRELATION_SHIFTED_UNIT_POLICY = ScientificPolicyRecord(
    id=ScientificPolicyId.PROFILE_CORRELATION_SHIFTED_UNIT,
    name="profile_correlation_v1",
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
        name="rank_weighted_motif_profile_fusion_v1",
        version="1",
        description=(
            "PhosR-inspired rank-weighted scoring that combines "
            "motif-frequency and profile-correlation scores using "
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
            "This is a PhosPy-specific scoring policy, not an exact PhosR "
            "implementation and not intended to provide numerical parity with PhosR.",
        ),
        output_scale="Relative downstream support score for kinase-site ranking.",
        quantitative_meaning="relative_support_score",
    )


def build_kinase_library_motif_scoring_policy(
    *,
    scoring_mode: str,
    resource_source_name: str | None,
    resource_source_version: str | None,
    resource_score_scale: str | None,
    workflow_score_scale: str,
    sequence_window: Mapping[str, object] | None = None,
) -> ScientificPolicyRecord:
    return ScientificPolicyRecord(
        id=ScientificPolicyId.KINASE_LIBRARY_MOTIF_SCORING,
        name="kinase_library_motif_scoring_v1",
        version="1",
        description=(
            "Scores phosphosite sequence windows against explicit Kinase "
            "Library-style position-specific kinase motif matrices."
        ),
        parameters={
            "scoring_mode": str(scoring_mode),
            "resource_source_name": resource_source_name,
            "resource_source_version": resource_source_version,
            "resource_score_scale": resource_score_scale,
            "workflow_score_scale": str(workflow_score_scale),
            "sequence_window": dict(sequence_window or {}),
        },
        assumptions=(
            "Motif matrices encode provider-scale relative residue-position support.",
            "Higher motif scores indicate stronger kinase motif compatibility.",
            "Workflow-normalized scores are relative within-run support values and "
            "are not calibrated probabilities.",
        ),
        output_scale=str(workflow_score_scale),
        quantitative_meaning="relative_motif_support_score",
    )


def build_kinase_profile_scoring_policy(
    *,
    profile_missing_value_strategy: str,
    min_substrates_floor: int,
    requested_min_substrates: int,
    profile_self_inclusion_policy: ProfileSelfInclusionPolicy | str = (
        ProfileSelfInclusionPolicy.ALLOW
    ),
) -> ScientificPolicyRecord:
    resolved_policy = ProfileSelfInclusionPolicy.parse(
        profile_self_inclusion_policy,
        field_name="kinase profile scoring policy profile_self_inclusion_policy",
    )
    leave_one_out_enabled = resolved_policy is ProfileSelfInclusionPolicy.LEAVE_ONE_OUT
    self_inclusion_allowed = resolved_policy is ProfileSelfInclusionPolicy.ALLOW
    return ScientificPolicyRecord(
        id=ScientificPolicyId.KINASE_PROFILE_SCORING,
        name="Kinase Profile Scoring Policy",
        version="1",
        description=(
            "Builds kinase reference profiles from quantified substrates and scores "
            "sites against those profiles with shifted Pearson correlation support."
        ),
        parameters={
            "profile_missing_value_strategy": str(profile_missing_value_strategy),
            "profile_self_inclusion_policy": resolved_policy.value,
            "self_inclusion_behavior": (
                "self_inclusion" if self_inclusion_allowed else "leave_one_out"
            ),
            "self_inclusion_allowed": bool(self_inclusion_allowed),
            "leave_one_out_enabled": bool(leave_one_out_enabled),
            "min_substrates_floor": int(min_substrates_floor),
            "requested_min_substrates": int(requested_min_substrates),
        },
        assumptions=(
            (
                "Profiles can include the same substrate site that is later "
                "scored when that site is present in the kinase profile "
                "definition."
                if self_inclusion_allowed
                else "Known substrate profile scores are recomputed without the "
                "scored site when leave-one-out support remains available."
            ),
            (
                "Leave-one-out profile recomputation is not applied in this policy."
                if self_inclusion_allowed
                else "Leave-one-out cells are left missing when excluding the "
                "scored substrate drops profile support below the configured "
                "minimum."
            ),
            "Profile missing-value strategy affects profile medians and can change "
            "site-level downstream support.",
        ),
        output_scale=(
            "Relative downstream support score in [0, 1] after shifted-correlation "
            "transformation."
        ),
        quantitative_meaning="relative_support_score",
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


__all__ = [
    "CandidateSubstrateSelectionPolicy",
    "KinaseProfileScoringPolicy",
    "PROFILE_CORRELATION_SHIFTED_UNIT_POLICY",
    "build_candidate_substrate_selection_policy",
    "build_kinase_profile_scoring_policy",
    "build_kinase_library_motif_scoring_policy",
    "build_motif_profile_rank_fusion_policy",
]
