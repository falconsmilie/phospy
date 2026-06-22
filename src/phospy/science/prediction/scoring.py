"""Scoring helpers shared by kinase workflow execution."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
import pandas as pd

from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.science.prediction.scientific_policies import (
    build_motif_profile_rank_fusion_policy,
)
from phospy.science.scoring.policy_models import DownstreamScoreSource

DOWNSTREAM_SCORE_SOURCE_PROFILE = DownstreamScoreSource.PROFILE_SCORES.value
DOWNSTREAM_SCORE_SOURCE_RANK_WEIGHTED_FUSION = (
    DownstreamScoreSource.RANK_WEIGHTED_FUSION_SCORES.value
)
KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE = "fused_motif_profile_evidence"
KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT = (
    "profile_only_motif_missing_or_constant"
)
KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP = "profile_only_no_motif_overlap"
KINASE_SCORE_SOURCE_UNAVAILABLE_NO_SCORE = "unavailable_no_score"
KINASE_SCORE_SOURCE_VALUES: tuple[str, ...] = (
    KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE,
    KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT,
    KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP,
    KINASE_SCORE_SOURCE_UNAVAILABLE_NO_SCORE,
)
KINASE_SCORE_SOURCE_SUMMARY_COLUMNS: tuple[str, ...] = (
    "fused_motif_profile_evidence_count",
    "profile_only_motif_missing_or_constant_count",
    "profile_only_no_motif_overlap_count",
    "unavailable_no_score_count",
    "sites_with_score_count",
    "total_sites_count",
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
            "Select PhosR-inspired rank-weighted fusion scores when available; "
            "otherwise use profile-only scores."
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


def build_kinase_score_source_diagnostics(
    *,
    motif_scores: pd.DataFrame,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Label per-cell downstream score evidence and summarize per kinase."""

    _validate_score_source_inputs(
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted_fusion_scores,
    )

    source_matrix = pd.DataFrame(
        KINASE_SCORE_SOURCE_UNAVAILABLE_NO_SCORE,
        index=rank_weighted_fusion_scores.index.copy(),
        columns=rank_weighted_fusion_scores.columns.copy(),
        dtype=object,
    )
    motif_columns = set(motif_scores.columns.astype(str))
    for kinase in rank_weighted_fusion_scores.columns.astype(str):
        profile_column = _column_series(profile_scores, kinase)
        fused_column = _column_series(rank_weighted_fusion_scores, kinase)
        available_profile = profile_column.notna()
        available_fused = fused_column.notna()
        if kinase in motif_columns:
            motif_column = _column_series(motif_scores, kinase).reindex(
                rank_weighted_fusion_scores.index
            )
            fused_mask = motif_column.notna() & available_profile & available_fused
            profile_only_missing_motif_mask = (
                motif_column.isna() & available_profile & available_fused
            )
            source_matrix.loc[fused_mask, kinase] = (
                KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE
            )
            source_matrix.loc[profile_only_missing_motif_mask, kinase] = (
                KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT
            )
            continue
        no_overlap_mask = available_profile & available_fused
        source_matrix.loc[no_overlap_mask, kinase] = (
            KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP
        )

    score_source_summary = build_kinase_score_source_summary(
        motif_scores=motif_scores,
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted_fusion_scores,
    )
    return source_matrix, score_source_summary


def build_kinase_score_source_summary(
    *,
    motif_scores: pd.DataFrame,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize per-kinase downstream score evidence without cell matrix."""

    _validate_score_source_inputs(
        profile_scores=profile_scores,
        rank_weighted_fusion_scores=rank_weighted_fusion_scores,
    )

    total_sites = int(rank_weighted_fusion_scores.shape[0])
    motif_columns = set(motif_scores.columns.astype(str))
    kinases = rank_weighted_fusion_scores.columns.astype(str).tolist()
    fused_counts: list[int] = []
    profile_only_missing_motif_counts: list[int] = []
    profile_only_no_overlap_counts: list[int] = []
    unavailable_counts: list[int] = []

    for kinase in kinases:
        profile_column = _column_series(profile_scores, kinase)
        fused_column = _column_series(rank_weighted_fusion_scores, kinase)
        available_profile = profile_column.notna().to_numpy(dtype=bool, copy=False)
        available_fused = fused_column.notna().to_numpy(dtype=bool, copy=False)
        available_both = available_profile & available_fused

        if kinase in motif_columns:
            motif_column = _column_series(motif_scores, kinase).reindex(
                rank_weighted_fusion_scores.index
            )
            motif_available = motif_column.notna().to_numpy(dtype=bool, copy=False)
            fused_count = int((motif_available & available_both).sum())
            profile_only_missing_motif_count = int(
                ((~motif_available) & available_both).sum()
            )
            profile_only_no_overlap_count = 0
        else:
            fused_count = 0
            profile_only_missing_motif_count = 0
            profile_only_no_overlap_count = int(available_both.sum())

        sites_with_score_count = (
            fused_count
            + profile_only_missing_motif_count
            + profile_only_no_overlap_count
        )
        unavailable_count = int(total_sites - sites_with_score_count)

        fused_counts.append(fused_count)
        profile_only_missing_motif_counts.append(profile_only_missing_motif_count)
        profile_only_no_overlap_counts.append(profile_only_no_overlap_count)
        unavailable_counts.append(unavailable_count)

    score_source_summary = pd.DataFrame(
        data={
            "fused_motif_profile_evidence_count": np.asarray(
                fused_counts,
                dtype="int64",
            ),
            "profile_only_motif_missing_or_constant_count": np.asarray(
                profile_only_missing_motif_counts,
                dtype="int64",
            ),
            "profile_only_no_motif_overlap_count": np.asarray(
                profile_only_no_overlap_counts,
                dtype="int64",
            ),
            "unavailable_no_score_count": np.asarray(
                unavailable_counts,
                dtype="int64",
            ),
        },
        index=rank_weighted_fusion_scores.columns.copy(),
    )
    score_source_summary.loc[:, "sites_with_score_count"] = (
        score_source_summary.loc[:, "fused_motif_profile_evidence_count"]
        + score_source_summary.loc[:, "profile_only_motif_missing_or_constant_count"]
        + score_source_summary.loc[:, "profile_only_no_motif_overlap_count"]
    ).astype("int64")
    score_source_summary.loc[:, "total_sites_count"] = int(
        rank_weighted_fusion_scores.shape[0]
    )
    score_source_summary.index.name = "kinase"
    score_source_summary = score_source_summary.loc[
        :, list(KINASE_SCORE_SOURCE_SUMMARY_COLUMNS)
    ]
    return score_source_summary


def _validate_score_source_inputs(
    *,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame,
) -> None:
    if not rank_weighted_fusion_scores.index.equals(profile_scores.index):
        raise ValueError(
            "rank_weighted_fusion_scores.index must match profile_scores.index"
        )
    missing_profile_kinases = [
        kinase
        for kinase in rank_weighted_fusion_scores.columns
        if kinase not in profile_scores.columns
    ]
    if missing_profile_kinases:
        raise ValueError(
            "rank_weighted_fusion_scores contains kinases missing from "
            f"profile_scores: {', '.join(map(str, missing_profile_kinases))}"
        )


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


def _column_series(frame: pd.DataFrame, column_name: str) -> pd.Series:
    values = frame.loc[:, column_name]
    if isinstance(values, pd.DataFrame):
        if values.shape[1] != 1:
            raise ValueError(
                f"expected one column for {column_name!r}; got {values.shape[1]}"
            )
        return values[values.columns[0]]
    return values


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

    This deterministic PhosR-inspired score-fusion policy is PhosPy-specific.
    It is not an exact PhosR implementation, numerical compatibility mode,
    enrichment statistic, classifier probability, or calibrated kinase activity
    estimate.
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
    "KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE",
    "KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT",
    "KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP",
    "KINASE_SCORE_SOURCE_SUMMARY_COLUMNS",
    "KINASE_SCORE_SOURCE_UNAVAILABLE_NO_SCORE",
    "KINASE_SCORE_SOURCE_VALUES",
    "DOWNSTREAM_SCORE_SOURCE_RANK_WEIGHTED_FUSION",
    "DOWNSTREAM_SCORE_SOURCE_PROFILE",
    "DownstreamScoreSelectionPolicy",
    "MotifProfileRankFusionPolicy",
    "SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY",
    "build_kinase_score_source_diagnostics",
    "fuse_profile_and_motif_scores_by_rank_weight",
    "resolve_downstream_score_matrix",
    "select_downstream_score_matrix",
]
