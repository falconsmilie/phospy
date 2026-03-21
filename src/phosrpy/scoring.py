from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .profiles import (
    AggregationMethod,
    KinaseProfileResult,
    build_kinase_substrate_profiles,
)


@dataclass(slots=True)
class KinaseScoringResult:
    profile_scores: pd.DataFrame
    combined_scores: pd.DataFrame | None = None
    weights: pd.DataFrame | None = None


@dataclass(slots=True)
class KinaseSubstrateScoreResult:
    profile_result: KinaseProfileResult
    profile_scores: pd.DataFrame
    motif_scores: pd.DataFrame | None = None
    combined_scores: pd.DataFrame | None = None
    ks_activity_matrix: pd.DataFrame | None = None
    weights: pd.DataFrame | None = None
    motif_sizes: pd.Series | None = None

    def to_scoring_result(self) -> KinaseScoringResult:
        return KinaseScoringResult(
            profile_scores=self.profile_scores,
            combined_scores=self.combined_scores,
            weights=self.weights,
        )


class KinaseScorer:
    """Score phosphosites against kinase-substrate activity profiles.

    This class provides a native Python seam for profile-based kinase scoring.
    Motif scoring remains separate, but profile construction, profile scoring,
    and motif/profile combination can now be orchestrated through the module's
    ``kinase_substrate_score()`` helper.
    """

    def __init__(self, kinase_profiles: pd.DataFrame) -> None:
        self.kinase_profiles = kinase_profiles.copy()

    @classmethod
    def from_substrate_map(
        cls,
        substrate_map: Mapping[str, Sequence[str]],
        phospho_matrix: pd.DataFrame,
        min_substrates: int = 1,
        aggregation: AggregationMethod = "median",
    ) -> KinaseScorer:
        profile_result = build_kinase_substrate_profiles(
            substrate_map=substrate_map,
            phospho_matrix=phospho_matrix,
            min_substrates=min_substrates,
            aggregation=aggregation,
        )
        return cls.from_profile_result(profile_result)

    @classmethod
    def from_profile_result(
        cls,
        profile_result: KinaseProfileResult,
    ) -> KinaseScorer:
        return cls(profile_result.profile_matrix)

    @classmethod
    def from_profile_dict(
        cls,
        kinase_profiles: Mapping[str, pd.Series | Sequence[float] | np.ndarray],
        sample_names: Sequence[str] | None = None,
    ) -> KinaseScorer:
        profile_frame = _build_profile_frame(
            kinase_profiles=kinase_profiles,
            sample_names=sample_names,
        )
        return cls(profile_frame)

    def score_phosphosite_profiles(self, phospho_matrix: pd.DataFrame) -> pd.DataFrame:
        """Return profile-based kinase scores on the 0..1 scale.

        The implementation mirrors the compatibility helper used in the R
        fixture scripts: Pearson correlations are computed between each
        phosphosite profile row and each kinase activity profile row, then
        rescaled from ``[-1, 1]`` to ``[0, 1]``.
        """

        sample_cols = list(phospho_matrix.columns)
        if set(sample_cols) != set(self.kinase_profiles.columns):
            msg = (
                "phospho_matrix columns must match kinase profile columns exactly; "
                "the order may differ, but the sets must be equal"
            )
            raise ValueError(msg)

        phospho_values = phospho_matrix.loc[:, sample_cols].to_numpy(dtype=float)
        profile_values = self.kinase_profiles.loc[:, sample_cols].to_numpy(dtype=float)

        correlation_matrix = _rowwise_correlation_matrix(
            left=phospho_values,
            right=profile_values,
        )
        score_matrix = (correlation_matrix + 1.0) / 2.0

        valid = ~np.isnan(score_matrix)
        score_matrix[valid] = np.clip(score_matrix[valid], 0.0, 1.0)

        return pd.DataFrame(
            score_matrix,
            index=phospho_matrix.index.copy(),
            columns=self.kinase_profiles.index.copy(),
        )

    def score(
        self,
        phospho_matrix: pd.DataFrame,
        motif_scores: pd.DataFrame | None = None,
        motif_sizes: pd.Series | None = None,
        profile_sizes: pd.Series | None = None,
        allow_profile_only_fallback: bool = False,
    ) -> KinaseScoringResult:
        profile_scores = self.score_phosphosite_profiles(phospho_matrix)

        if motif_scores is None:
            return KinaseScoringResult(profile_scores=profile_scores)

        if motif_sizes is None or profile_sizes is None:
            msg = (
                "motif_sizes and profile_sizes are required when motif_scores are "
                "provided"
            )
            raise ValueError(msg)

        combined_scores, weights = combine_profile_and_motif_scores(
            motif_scores=motif_scores,
            profile_scores=profile_scores,
            motif_sizes=motif_sizes,
            profile_sizes=profile_sizes,
            allow_profile_only_fallback=allow_profile_only_fallback,
        )
        return KinaseScoringResult(
            profile_scores=profile_scores,
            combined_scores=combined_scores,
            weights=weights,
        )


def kinase_substrate_score(
    substrate_map: Mapping[str, Sequence[str]],
    phospho_matrix: pd.DataFrame,
    motif_scores: pd.DataFrame | None = None,
    motif_sizes: pd.Series | None = None,
    min_motif_size: int = 1,
    min_substrates: int = 1,
    aggregation: AggregationMethod = "median",
    allow_profile_only_fallback: bool = False,
) -> KinaseSubstrateScoreResult:
    """Build kinase profiles, score phosphosites, and combine motif evidence.

    This is a narrower Python-native analogue of PhosR's
    ``kinaseSubstrateScore()``. It builds kinase substrate profiles from the
    provided substrate map, scores all phosphosites against those profiles, and
    optionally combines those profile scores with a caller-supplied motif score
    matrix. Motif *generation* is not implemented here yet; callers provide the
    motif scores and motif sizes explicitly when they want a combined matrix.
    """

    _validate_positive_int(min_motif_size, name="min_motif_size")

    profile_result = build_kinase_substrate_profiles(
        substrate_map=substrate_map,
        phospho_matrix=phospho_matrix,
        min_substrates=min_substrates,
        aggregation=aggregation,
    )
    scorer = KinaseScorer.from_profile_result(profile_result)
    profile_scores = scorer.score_phosphosite_profiles(phospho_matrix)
    profile_sizes = profile_result.substrate_counts.loc[profile_scores.columns].astype(
        float
    )

    if motif_scores is None:
        combined_scores, weights = _profile_only_score_outputs(
            profile_scores=profile_scores,
            allow_profile_only_fallback=allow_profile_only_fallback,
        )
        ks_activity_matrix = profile_result.profile_matrix.loc[
            profile_scores.columns
        ].copy()
        return KinaseSubstrateScoreResult(
            profile_result=profile_result,
            profile_scores=profile_scores,
            combined_scores=combined_scores,
            ks_activity_matrix=ks_activity_matrix,
            weights=weights,
        )

    if motif_sizes is None:
        msg = "motif_sizes are required when motif_scores are provided"
        raise ValueError(msg)

    _require_index_members(motif_sizes, motif_scores.columns, name="motif_sizes")

    allowed_kinases = [
        kinase
        for kinase in motif_scores.columns
        if float(motif_sizes.loc[kinase]) >= float(min_motif_size)
    ]
    filtered_motif_scores = motif_scores.loc[:, allowed_kinases].copy()
    filtered_motif_sizes = motif_sizes.loc[allowed_kinases].astype(float).copy()

    combined_scores, weights = combine_profile_and_motif_scores(
        motif_scores=filtered_motif_scores,
        profile_scores=profile_scores,
        motif_sizes=filtered_motif_sizes,
        profile_sizes=profile_sizes,
        allow_profile_only_fallback=allow_profile_only_fallback,
    )

    ks_activity_columns = (
        list(combined_scores.columns)
        if combined_scores is not None
        else list(profile_scores.columns)
    )
    ks_activity_matrix = profile_result.profile_matrix.loc[ks_activity_columns].copy()

    return KinaseSubstrateScoreResult(
        profile_result=profile_result,
        profile_scores=profile_scores,
        motif_scores=filtered_motif_scores,
        combined_scores=combined_scores,
        ks_activity_matrix=ks_activity_matrix,
        weights=weights,
        motif_sizes=filtered_motif_sizes,
    )


def combine_profile_and_motif_scores(
    motif_scores: pd.DataFrame,
    profile_scores: pd.DataFrame,
    motif_sizes: pd.Series,
    profile_sizes: pd.Series,
    allow_profile_only_fallback: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine motif and profile score matrices using rank-derived weights."""

    overlap = [
        kinase
        for kinase in motif_scores.columns
        if kinase in set(profile_scores.columns)
    ]

    if not overlap:
        if not allow_profile_only_fallback:
            msg = "No overlapping kinases between motif and profile score matrices"
            raise ValueError(msg)

        weights = pd.DataFrame(
            {
                "motif_weight": 0.0,
                "profile_weight": 1.0,
                "motif_rank_weight": 0.0,
                "profile_rank_weight": 1.0,
            },
            index=profile_scores.columns.copy(),
        )
        weights.index.name = "kinase"
        return profile_scores.copy(), weights

    _require_index_members(motif_sizes, overlap, name="motif_sizes")
    _require_index_members(profile_sizes, overlap, name="profile_sizes")

    motif_rank_weight = np.log(motif_sizes.loc[overlap].rank(method="average") + 1.0)
    profile_rank_weight = np.log(
        profile_sizes.loc[overlap].rank(method="average") + 1.0
    )
    total_weight = motif_rank_weight + profile_rank_weight

    combined_scores = (
        motif_scores.loc[:, overlap].multiply(motif_rank_weight, axis=1)
        + profile_scores.loc[:, overlap].multiply(profile_rank_weight, axis=1)
    ).divide(total_weight, axis=1)

    weights = pd.DataFrame(
        {
            "motif_weight": motif_rank_weight / total_weight,
            "profile_weight": profile_rank_weight / total_weight,
            "motif_rank_weight": motif_rank_weight,
            "profile_rank_weight": profile_rank_weight,
        },
        index=overlap,
    )
    weights.index.name = "kinase"
    return combined_scores, weights


def _build_profile_frame(
    kinase_profiles: Mapping[str, pd.Series | Sequence[float] | np.ndarray],
    sample_names: Sequence[str] | None,
) -> pd.DataFrame:
    rows: dict[str, pd.Series] = {}

    for kinase, profile in kinase_profiles.items():
        if isinstance(profile, pd.Series):
            rows[kinase] = profile.astype(float)
            continue

        if sample_names is None:
            msg = "sample_names are required when kinase profiles are not pandas Series"
            raise ValueError(msg)

        rows[kinase] = pd.Series(profile, index=list(sample_names), dtype=float)

    profile_frame = pd.DataFrame.from_dict(rows, orient="index")
    profile_frame.index.name = "kinase"
    return profile_frame


def _profile_only_score_outputs(
    profile_scores: pd.DataFrame,
    allow_profile_only_fallback: bool,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    if not allow_profile_only_fallback:
        return None, None

    weights = pd.DataFrame(
        {
            "motif_weight": 0.0,
            "profile_weight": 1.0,
            "motif_rank_weight": 0.0,
            "profile_rank_weight": 1.0,
        },
        index=profile_scores.columns.copy(),
    )
    weights.index.name = "kinase"
    return profile_scores.copy(), weights


def _require_index_members(
    series: pd.Series, members: Sequence[str], name: str
) -> None:
    missing = [member for member in members if member not in series.index]
    if missing:
        msg = f"{name} is missing entries for: {', '.join(missing)}"
        raise ValueError(msg)


def _rowwise_correlation_matrix(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left_centered = left - left.mean(axis=1, keepdims=True)
    right_centered = right - right.mean(axis=1, keepdims=True)

    left_scale = np.linalg.norm(left_centered, axis=1)
    right_scale = np.linalg.norm(right_centered, axis=1)
    denominator = np.outer(left_scale, right_scale)

    with np.errstate(divide="ignore", invalid="ignore"):
        correlation = left_centered @ right_centered.T / denominator

    correlation[denominator == 0.0] = np.nan
    return correlation


def _validate_positive_int(value: int, name: str) -> None:
    if value < 1:
        raise ValueError(f"{name} must be at least 1")
