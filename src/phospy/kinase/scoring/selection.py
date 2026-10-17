"""Scoring helpers shared by kinase workflow execution."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from phospy.kinase.result_boundary import (
        DownstreamScoreSelectionPolicy as _DownstreamScoreSelectionPolicy,
    )
    from phospy.kinase.result_boundary import (
        DownstreamScoreSource as _DownstreamScoreSource,
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


def select_downstream_score_matrix(
    *,
    profile_scores: pd.DataFrame,
    rank_weighted_fusion_scores: pd.DataFrame | None,
) -> tuple[pd.DataFrame, _DownstreamScoreSource]:
    """Resolve the authoritative downstream prediction score matrix."""

    from phospy.kinase.result_boundary import DownstreamScoreSource

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
) -> tuple[pd.DataFrame, _DownstreamScoreSource, _DownstreamScoreSelectionPolicy]:
    """Resolve downstream score matrix/source and the active selection policy."""

    from phospy.kinase.result_boundary import (
        SIGNALOME_DOWNSTREAM_SCORE_RANK_WEIGHTED_PREFERRED_POLICY,
    )

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


__all__ = [
    "KINASE_SCORE_SOURCE_FUSED_MOTIF_PROFILE_EVIDENCE",
    "KINASE_SCORE_SOURCE_PROFILE_ONLY_MOTIF_MISSING_OR_CONSTANT",
    "KINASE_SCORE_SOURCE_PROFILE_ONLY_NO_MOTIF_OVERLAP",
    "KINASE_SCORE_SOURCE_SUMMARY_COLUMNS",
    "KINASE_SCORE_SOURCE_UNAVAILABLE_NO_SCORE",
    "KINASE_SCORE_SOURCE_VALUES",
    "build_kinase_score_source_diagnostics",
    "build_kinase_score_source_summary",
    "resolve_downstream_score_matrix",
    "select_downstream_score_matrix",
]
