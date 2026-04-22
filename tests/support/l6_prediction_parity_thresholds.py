from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class RankingGateThresholds:
    mean_spearman_rank_corr_floor: float
    mean_top20_overlap_floor: float
    mean_top30_overlap_floor: float
    good_top10_count_floor: int


@dataclass(frozen=True, slots=True)
class ScoreDivergenceThresholds:
    score_correlation_floor: float


# 1.5.0 release-governance bars for repaired like-for-like ranking surfaces.
L6_PREDICTION_MATRIX_RANKING_GATES = RankingGateThresholds(
    mean_spearman_rank_corr_floor=0.99,
    mean_top20_overlap_floor=0.95,
    mean_top30_overlap_floor=0.95,
    good_top10_count_floor=24,
)

L6_RANKED_TOPK_EXPORT_RANKING_GATES = RankingGateThresholds(
    mean_spearman_rank_corr_floor=0.99,
    mean_top20_overlap_floor=0.95,
    mean_top30_overlap_floor=0.95,
    good_top10_count_floor=24,
)

L6_CROSS_POLICY_PREDICTION_MATRIX_DIVERGENCE_GATES = ScoreDivergenceThresholds(
    score_correlation_floor=0.95,
)
