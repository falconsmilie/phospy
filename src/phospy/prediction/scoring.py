"""Scoring helpers shared by kinase workflow execution."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
import pandas as pd

from phospy.policy_models import DownstreamScoreSource
from phospy.scientific_policies import (
    ScientificPolicyId,
    ScientificPolicyRecord,
    build_motif_profile_rank_fusion_policy,
)

DOWNSTREAM_SCORE_SOURCE_PROFILE = DownstreamScoreSource.PROFILE_SCORES.value
DOWNSTREAM_SCORE_SOURCE_RANK_WEIGHTED_FUSION = (
    DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES.value
)


@dataclass(frozen=True, slots=True)
class DownstreamScoreSelectionPolicy:
    """Versioned policy describing downstream score-lane selection semantics."""

    name: str
    version: str
    parameters: Mapping[str, object]
    description: str

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
            id=ScientificPolicyId.SIGNALOME_DOWNSTREAM_SCORE_SELECTION,
            name=self.name,
            version=self.version,
            description=self.description,
            parameters=self.parameters,
            assumptions=(
                "Downstream score-lane selection changes which upstream signal is "
                "interpreted by signalome stages.",
            ),
            output_scale="Selected downstream kinase-site score matrix.",
            quantitative_meaning="relative_downstream_support",
        )


SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY = (
    DownstreamScoreSelectionPolicy(
        name="signalome_downstream_score_rank_weighted_preferred_v1",
        version="1",
        parameters={
            "preferred_source": DOWNSTREAM_SCORE_SOURCE_RANK_WEIGHTED_FUSION,
            "fallback_source": DOWNSTREAM_SCORE_SOURCE_PROFILE,
        },
        description=(
            "Select rank-weighted fusion scores when available; otherwise use "
            "profile-only scores."
        ),
    )
)


def select_downstream_score_matrix(
    *,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame | None,
) -> tuple[pd.DataFrame, DownstreamScoreSource]:
    """Resolve the authoritative downstream prediction score matrix."""

    if rank_weighted_fusion_scores is not None:
        return (
            rank_weighted_fusion_scores,
            DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES,
        )
    return profile_scores, DownstreamScoreSource.PROFILE_SCORES


def resolve_downstream_score_matrix(
    *,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame | None,
) -> tuple[pd.DataFrame, DownstreamScoreSource, DownstreamScoreSelectionPolicy]:
    """Resolve downstream score matrix/source and the active selection policy."""

    selected, source = select_downstream_score_matrix(
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted_fusion_scores,
    )
    return (
        selected,
        source,
        SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY,
    )


def fuse_profile_and_motif_scores_by_rank_weight(
    *,
    motif_scores: pd.DataFrame,
    profile_scores: pd.DataFrame,
    motif_sizes: pd.Series,
    profile_sizes: pd.Series,
    allow_profile_only_fallback: bool = True,
    emit_weights: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Compatibility wrapper around the rank-fusion policy object."""

    policy = MotifProfileRankFusionPolicy(
        allow_profile_only_fallback=allow_profile_only_fallback,
        emit_weights=emit_weights,
    )
    return policy.fuse(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        motif_sizes=motif_sizes,
        profile_sizes=profile_sizes,
    )


@dataclass(frozen=True, slots=True)
class MotifProfileRankFusionPolicy:
    """Executable rank-fusion policy for motif and profile evidence."""

    allow_profile_only_fallback: bool = True
    emit_weights: bool = True

    @property
    def record(self) -> ScientificPolicyRecord:
        return build_motif_profile_rank_fusion_policy(
            allow_profile_only_fallback=bool(self.allow_profile_only_fallback),
            emit_weights=bool(self.emit_weights),
        )

    def fuse(
        self,
        *,
        motif_scores: pd.DataFrame,
        profile_scores: pd.DataFrame,
        motif_sizes: pd.Series,
        profile_sizes: pd.Series,
    ) -> tuple[pd.DataFrame, pd.DataFrame | None]:
        return _fuse_profile_and_motif_scores_by_rank_weight(
            motif_scores=motif_scores,
            profile_scores=profile_scores,
            motif_sizes=motif_sizes,
            profile_sizes=profile_sizes,
            allow_profile_only_fallback=self.allow_profile_only_fallback,
            emit_weights=self.emit_weights,
        )


def _fuse_profile_and_motif_scores_by_rank_weight(
    *,
    motif_scores: pd.DataFrame,
    profile_scores: pd.DataFrame,
    motif_sizes: pd.Series,
    profile_sizes: pd.Series,
    allow_profile_only_fallback: bool,
    emit_weights: bool,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Fuse motif-frequency and profile-correlation scores by rank-derived weights.

    Fuse motif-frequency scores and profile-correlation scores using rank-derived
    weights from motif library size and quantified substrate count.

    This is a deterministic score-fusion policy. It is not an enrichment
    statistic, classifier probability, or calibrated kinase activity estimate.
    """

    profile_kinases = set(profile_scores.columns)
    overlap = [kinase for kinase in motif_scores.columns if kinase in profile_kinases]

    if not overlap:
        if not allow_profile_only_fallback:
            raise ValueError(
                "no overlapping kinases between motif_scores and profile_scores"
            )
        weights = (
            _profile_only_weights(tuple(str(k) for k in profile_scores.columns))
            if emit_weights
            else None
        )
        return profile_scores.copy(), weights

    _require_index_members(motif_sizes, overlap, name="motif_sizes")
    _require_index_members(profile_sizes, overlap, name="profile_sizes")

    motif_rank_weight = np.log(motif_sizes.loc[overlap].rank(method="average") + 1.0)
    profile_rank_weight = np.log(
        profile_sizes.loc[overlap].rank(method="average") + 1.0
    )
    total_weight = motif_rank_weight + profile_rank_weight
    rank_weighted_fusion_scores = (
        motif_scores.loc[:, overlap].multiply(motif_rank_weight, axis=1)
        + profile_scores.loc[:, overlap].multiply(profile_rank_weight, axis=1)
    ).divide(total_weight, axis=1)
    # Motif columns can become all-NaN after min-max scaling when no sequence
    # contrast exists. Preserve profile evidence in that case instead of
    # discarding the kinase column for downstream prediction.
    motif_missing_profile_present = (
        motif_scores.loc[:, overlap].isna() & profile_scores.loc[:, overlap].notna()
    )
    rank_weighted_fusion_scores = rank_weighted_fusion_scores.where(
        ~motif_missing_profile_present,
        profile_scores.loc[:, overlap],
    )
    weights: pd.DataFrame | None
    if emit_weights:
        weights = pd.DataFrame(
            {
                "motif_weight": motif_rank_weight / total_weight,
                "profile_weight": profile_rank_weight / total_weight,
                "motif_rank_weight": motif_rank_weight,
                "profile_rank_weight": profile_rank_weight,
            },
            index=pd.Index(overlap),
        )
    else:
        weights = None

    if allow_profile_only_fallback:
        profile_only = [
            kinase for kinase in profile_scores.columns if kinase not in set(overlap)
        ]
        if profile_only:
            rank_weighted_fusion_scores = pd.concat(
                [rank_weighted_fusion_scores, profile_scores.loc[:, profile_only]],
                axis=1,
            )
            if weights is not None:
                weights = pd.concat(
                    [weights, _profile_only_weights(profile_only)],
                    axis=0,
                )

    if weights is not None:
        weights.index.name = "kinase"
    return rank_weighted_fusion_scores, weights


def _require_index_members(
    series: pd.Series,
    members: Sequence[str],
    *,
    name: str,
) -> None:
    missing = [member for member in members if member not in series.index]
    if missing:
        raise ValueError(f"{name} is missing entries for: {', '.join(missing)}")


def _profile_only_weights(kinases: Sequence[str]) -> pd.DataFrame:
    weights = pd.DataFrame(
        {
            "motif_weight": 0.0,
            "profile_weight": 1.0,
            "motif_rank_weight": 0.0,
            "profile_rank_weight": 1.0,
        },
        index=pd.Index(kinases),
    )
    weights.index.name = "kinase"
    return weights


__all__ = [
    "DOWNSTREAM_SCORE_SOURCE_RANK_WEIGHTED_FUSION",
    "DOWNSTREAM_SCORE_SOURCE_PROFILE",
    "DownstreamScoreSelectionPolicy",
    "MotifProfileRankFusionPolicy",
    "SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY",
    "fuse_profile_and_motif_scores_by_rank_weight",
    "resolve_downstream_score_matrix",
    "select_downstream_score_matrix",
]
