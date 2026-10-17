"""Scientific policy records for kinase scoring."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.policies import PolicyEnum
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)


class ProfileSelfInclusionPolicy(PolicyEnum):
    ALLOW = "allow"
    LEAVE_ONE_OUT = "leave_one_out"


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


__all__ = [
    "KinaseLibraryMotifScoringPolicy",
    "KinaseProfileScoringPolicy",
    "PROFILE_CORRELATION_SHIFTED_UNIT_POLICY",
    "ProfileSelfInclusionPolicy",
    "build_kinase_library_motif_scoring_policy",
    "build_kinase_profile_scoring_policy",
    "build_motif_profile_rank_fusion_policy",
]
