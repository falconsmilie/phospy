from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..validation.errors import InputCompatibilityError


@dataclass(slots=True)
class KinaseScoringResult:
    """Profile and combined kinase scoring tables for one scoring run."""

    profile_scores: pd.DataFrame
    combined_scores: pd.DataFrame | None = None
    weights: pd.DataFrame | None = None


class KinaseScorer:
    """Score phosphosites against kinase-substrate activity profiles.

    This class provides a first native Python scoring seam for profile-based
    kinase scoring. It intentionally stays narrower than PhosR's full
    ``kinaseSubstrateScore()`` workflow: motif scoring and prediction remain
    separate concerns.
    """

    def __init__(
        self,
        kinase_profiles: pd.DataFrame,
        correlation_batch_size: int | None = 1024,
    ) -> None:
        self.kinase_profiles = kinase_profiles.copy()
        self.correlation_batch_size = _validate_correlation_batch_size(
            correlation_batch_size
        )

    @classmethod
    def from_profile_dict(
        cls,
        kinase_profiles: Mapping[str, pd.Series | Sequence[float] | np.ndarray],
        sample_names: Sequence[str] | None = None,
        correlation_batch_size: int | None = 1024,
    ) -> KinaseScorer:
        profile_frame = _build_profile_frame(
            kinase_profiles=kinase_profiles,
            sample_names=sample_names,
        )
        return cls(
            profile_frame,
            correlation_batch_size=correlation_batch_size,
        )

    def score_phosphosite_profiles(
        self,
        phospho_matrix: pd.DataFrame,
        *,
        correlation_batch_size: int | None = None,
    ) -> pd.DataFrame:
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
            raise InputCompatibilityError(msg)

        phospho_values = phospho_matrix.loc[:, sample_cols].to_numpy(dtype=float)
        profile_values = self.kinase_profiles.loc[:, sample_cols].to_numpy(dtype=float)

        correlation_matrix = _rowwise_correlation_matrix(
            left=phospho_values,
            right=profile_values,
            batch_size=_resolve_correlation_batch_size(
                requested_batch_size=correlation_batch_size,
                default_batch_size=self.correlation_batch_size,
            ),
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
        *,
        correlation_batch_size: int | None = None,
    ) -> KinaseScoringResult:
        profile_scores = self.score_phosphosite_profiles(
            phospho_matrix,
            correlation_batch_size=correlation_batch_size,
        )

        if motif_scores is None:
            return KinaseScoringResult(profile_scores=profile_scores)

        if motif_sizes is None or profile_sizes is None:
            msg = (
                "motif_sizes and profile_sizes are required when motif_scores are "
                "provided"
            )
            raise InputCompatibilityError(msg)

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


def combine_profile_and_motif_scores(
    motif_scores: pd.DataFrame,
    profile_scores: pd.DataFrame,
    motif_sizes: pd.Series,
    profile_sizes: pd.Series,
    allow_profile_only_fallback: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine motif and profile scores using rank-derived weights.

    Overlapping kinases are combined using motif and profile rank-derived weights.
    When ``allow_profile_only_fallback`` is ``True``, kinases present only in the
    profile score matrix are preserved with ``motif_weight=0`` and
    ``profile_weight=1`` instead of being dropped during partial overlap.
    """

    profile_kinases = set(profile_scores.columns)
    overlap = [kinase for kinase in motif_scores.columns if kinase in profile_kinases]

    if not overlap:
        if not allow_profile_only_fallback:
            msg = "No overlapping kinases between motif and profile score matrices"
            raise InputCompatibilityError(msg)

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

    if allow_profile_only_fallback:
        profile_only = [
            kinase for kinase in profile_scores.columns if kinase not in set(overlap)
        ]
        if profile_only:
            combined_scores = pd.concat(
                [combined_scores, profile_scores.loc[:, profile_only]],
                axis=1,
            )
            fallback_weights = pd.DataFrame(
                {
                    "motif_weight": 0.0,
                    "profile_weight": 1.0,
                    "motif_rank_weight": 0.0,
                    "profile_rank_weight": 1.0,
                },
                index=profile_only,
            )
            weights = pd.concat([weights, fallback_weights], axis=0)

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
            raise InputCompatibilityError(msg)

        rows[kinase] = pd.Series(profile, index=list(sample_names), dtype=float)

    profile_frame = pd.DataFrame.from_dict(rows, orient="index")
    profile_frame.index.name = "kinase"
    return profile_frame


def _require_index_members(
    series: pd.Series, members: Sequence[str], name: str
) -> None:
    missing = [member for member in members if member not in series.index]
    if missing:
        msg = f"{name} is missing entries for: {', '.join(missing)}"
        raise InputCompatibilityError(msg)


def _rowwise_correlation_matrix(
    left: np.ndarray,
    right: np.ndarray,
    batch_size: int | None = None,
) -> np.ndarray:
    resolved_batch_size = _validate_correlation_batch_size(batch_size)
    if resolved_batch_size is None:
        resolved_batch_size = left.shape[0] or 1

    right_centered = right - right.mean(axis=1, keepdims=True)
    right_scale = np.linalg.norm(right_centered, axis=1)

    correlation = np.empty((left.shape[0], right.shape[0]), dtype=float)

    for start in range(0, left.shape[0], resolved_batch_size):
        stop = min(start + resolved_batch_size, left.shape[0])
        left_chunk = left[start:stop]
        left_centered = left_chunk - left_chunk.mean(axis=1, keepdims=True)
        left_scale = np.linalg.norm(left_centered, axis=1)
        denominator = np.outer(left_scale, right_scale)

        with np.errstate(divide="ignore", invalid="ignore"):
            chunk_correlation = left_centered @ right_centered.T / denominator

        chunk_correlation[denominator == 0.0] = np.nan
        correlation[start:stop, :] = chunk_correlation

    return correlation


def _resolve_correlation_batch_size(
    requested_batch_size: int | None,
    default_batch_size: int | None,
) -> int | None:
    if requested_batch_size is None:
        return default_batch_size
    return _validate_correlation_batch_size(requested_batch_size)


def _validate_correlation_batch_size(batch_size: int | None) -> int | None:
    if batch_size is None:
        return None
    if batch_size < 1:
        msg = "correlation_batch_size must be at least 1 when provided"
        raise InputCompatibilityError(msg)
    return batch_size
