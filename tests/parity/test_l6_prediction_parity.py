from __future__ import annotations

import pytest

from tests.support.l6_prediction_parity_metrics import (
    collect_l6_prediction_parity_metrics,
)
from tests.support.parity_reporting import (
    format_fraction,
    format_percent,
    format_shape,
    record_parity_metrics,
)

pytestmark = pytest.mark.parity

PROFILE_MAE_CEILING = 1e-12
PROFILE_PEARSON_FLOOR = 0.999999
PROFILE_SPEARMAN_FLOOR = 0.999999
COMBINED_MAE_CEILING = 1e-12
COMBINED_PEARSON_FLOOR = 0.999999
COMBINED_SPEARMAN_FLOOR = 0.999999
WEIGHTS_MAE_CEILING = 1e-12
WEIGHTS_MAX_ABS_DIFF_CEILING = 1e-11
PREDICTION_MATRIX_MAE_CEILING = 1e-12
PREDICTION_MATRIX_PEARSON_FLOOR = 0.999999
PREDICTION_MATRIX_SPEARMAN_FLOOR = 0.999999
CANDIDATE_OVERLAP_PRECISION_FLOOR = 0.999999
CANDIDATE_OVERLAP_RECALL_FLOOR = 0.999999
CANDIDATE_OVERLAP_F1_FLOOR = 0.999999
# Hard regression gates for promoted-reference ranking agreement.
# These bars intentionally mirror the legacy release thresholds now that the
# ranking comparison surface is like-for-like (predMat-derived ranking compared
# to promoted predMat reference ranking, and ranked top-k export compared to
# promoted top-k reference ranking).
PRED_MAT_RANK_SPEARMAN_FLOOR = 0.96
PRED_MAT_TOP20_OVERLAP_FLOOR = 0.85
PRED_MAT_TOP30_OVERLAP_FLOOR = 0.88
PRED_MAT_GOOD_TOP10_COUNT_FLOOR = 20
TOPK_EXPORT_RANK_SPEARMAN_FLOOR = 0.96
TOPK_EXPORT_TOP20_OVERLAP_FLOOR = 0.85
TOPK_EXPORT_TOP30_OVERLAP_FLOOR = 0.88
TOPK_EXPORT_GOOD_TOP10_COUNT_FLOOR = 20
CROSS_POLICY_PRED_MAT_CORRELATION_FLOOR = 0.95


def test_l6_scoring_tables_parity_matches_promoted_reference_surfaces() -> None:
    metrics = collect_l6_prediction_parity_metrics()

    assert metrics.profile.shared_site_count == metrics.profile.observed_shape[0], (
        "profile parity regressed on promoted-reference shared site coverage"
    )
    assert metrics.profile.expected_shape[0] == metrics.profile.observed_shape[0]
    assert metrics.profile.shared_kinase_count == metrics.profile.observed_shape[1]
    assert metrics.profile.expected_shape[1] >= metrics.profile.observed_shape[1]
    assert metrics.profile.mean_abs_diff <= PROFILE_MAE_CEILING, (
        "L6 profile scoring parity regressed against promoted reference surfaces"
    )
    assert metrics.profile.mean_pearson_corr >= PROFILE_PEARSON_FLOOR
    assert metrics.profile.mean_spearman_corr >= PROFILE_SPEARMAN_FLOOR

    assert metrics.combined.shared_site_count == metrics.combined.observed_shape[0]
    assert metrics.combined.expected_shape[0] == metrics.combined.observed_shape[0]
    assert metrics.combined.shared_kinase_count >= 15
    assert metrics.combined.shared_kinase_count <= metrics.combined.observed_shape[1]
    assert metrics.combined.expected_shape[1] >= metrics.combined.observed_shape[1]
    assert metrics.combined.mean_abs_diff <= COMBINED_MAE_CEILING, (
        "L6 combined-score parity regressed against promoted reference surfaces"
    )
    assert metrics.combined.mean_pearson_corr >= COMBINED_PEARSON_FLOOR
    assert metrics.combined.mean_spearman_corr >= COMBINED_SPEARMAN_FLOOR

    assert metrics.weights.shared_site_count >= 15
    assert metrics.weights.shared_site_count <= metrics.weights.observed_shape[0]
    assert metrics.weights.shared_kinase_count == metrics.weights.observed_shape[1]
    assert metrics.weights.mean_abs_diff <= WEIGHTS_MAE_CEILING
    assert metrics.weights.max_abs_diff <= WEIGHTS_MAX_ABS_DIFF_CEILING


def test_l6_prediction_matrix_numeric_parity_matches_promoted_reference_surfaces() -> (
    None
):
    metrics = collect_l6_prediction_parity_metrics()
    parity = metrics.prediction_matrix

    assert parity.shared_site_count == parity.observed_shape[0], (
        "prediction-matrix parity regressed on promoted-reference shared site coverage"
    )
    assert parity.expected_shape[0] == parity.observed_shape[0]
    assert parity.shared_kinase_count == parity.observed_shape[1]
    assert parity.expected_shape[1] >= parity.observed_shape[1]
    assert parity.mean_abs_diff <= PREDICTION_MATRIX_MAE_CEILING, (
        "prediction-matrix parity regressed on promoted-reference MAE"
    )
    assert parity.mean_pearson_corr >= PREDICTION_MATRIX_PEARSON_FLOOR
    assert parity.mean_spearman_corr >= PREDICTION_MATRIX_SPEARMAN_FLOOR


def test_l6_prediction_matrix_ranking_parity_matches_promoted_reference_surfaces() -> (
    None
):
    ranking = collect_l6_prediction_parity_metrics().prediction_matrix_ranking

    assert ranking.kinases_compared > 0
    assert ranking.mean_spearman_rank_corr >= PRED_MAT_RANK_SPEARMAN_FLOOR, (
        "prediction-matrix ranking parity regressed on per-kinase Spearman agreement"
    )
    assert ranking.mean_top20_overlap >= PRED_MAT_TOP20_OVERLAP_FLOOR, (
        "prediction-matrix ranking parity regressed on top-20 overlap"
    )
    assert ranking.mean_top30_overlap >= PRED_MAT_TOP30_OVERLAP_FLOOR
    assert ranking.good_top10_count >= PRED_MAT_GOOD_TOP10_COUNT_FLOOR, (
        "prediction-matrix ranking parity regressed on top-10 agreement coverage"
    )
    assert ranking.top_rank_total == ranking.kinases_compared


def test_l6_candidate_selection_parity_matches_promoted_reference_surfaces() -> None:
    metrics = collect_l6_prediction_parity_metrics()
    candidates = metrics.candidates

    assert candidates.observed_rows > 0
    assert candidates.expected_rows > 0
    assert candidates.overlap_rows > 0
    assert candidates.overlap_precision >= CANDIDATE_OVERLAP_PRECISION_FLOOR, (
        "candidate-set parity regressed on promoted-reference precision"
    )
    assert candidates.overlap_recall >= CANDIDATE_OVERLAP_RECALL_FLOOR, (
        "candidate-set parity regressed on promoted-reference recall"
    )
    assert candidates.overlap_f1 >= CANDIDATE_OVERLAP_F1_FLOOR, (
        "candidate-set parity regressed on promoted-reference F1 overlap"
    )
    assert candidates.shared_kinase_count > 0


def test_l6_ranked_topk_export_parity_matches_promoted_reference_surfaces() -> None:
    metrics = collect_l6_prediction_parity_metrics()
    ranking = metrics.ranked_topk_export

    assert ranking.kinases_compared > 0
    assert ranking.mean_spearman_rank_corr >= TOPK_EXPORT_RANK_SPEARMAN_FLOOR, (
        "ranked top-k export parity regressed on per-kinase Spearman agreement"
    )
    assert ranking.mean_top20_overlap >= TOPK_EXPORT_TOP20_OVERLAP_FLOOR, (
        "top-k export parity regressed on top-20 overlap"
    )
    assert ranking.mean_top30_overlap >= TOPK_EXPORT_TOP30_OVERLAP_FLOOR
    assert ranking.good_top10_count >= TOPK_EXPORT_GOOD_TOP10_COUNT_FLOOR
    assert ranking.top_rank_total == ranking.kinases_compared


def test_l6_cross_policy_prediction_matrix_divergence_stable_vs_r_parity() -> None:
    metrics = collect_l6_prediction_parity_metrics()
    divergence = metrics.policy_divergence
    pred_mat = divergence.prediction_matrix_ranking

    assert (
        divergence.prediction_matrix_score_corr
        >= CROSS_POLICY_PRED_MAT_CORRELATION_FLOOR
    ), (
        "cross-policy prediction-matrix divergence regressed on score correlation "
        "(stable vs r_parity)"
    )
    assert pred_mat.kinases_compared > 0
    assert pred_mat.top_rank_total == pred_mat.kinases_compared


def test_l6_cross_policy_ranked_topk_divergence_stable_vs_r_parity() -> None:
    topk = collect_l6_prediction_parity_metrics().policy_divergence.ranked_topk_export

    assert topk.kinases_compared > 0
    assert topk.top_rank_total == topk.kinases_compared


@pytest.mark.parity_diagnostic
def test_l6_prediction_parity_reporting_is_surface_explicit(
    request: pytest.FixtureRequest,
) -> None:
    metrics = collect_l6_prediction_parity_metrics()
    divergence = metrics.policy_divergence

    record_parity_metrics(
        request.config,
        family="l6_prediction",
        metrics=[
            (
                "contract: rewrite-vs-promoted-reference | scoring-table surfaces",
                "asserted",
            ),
            (
                "profile score shape (rewrite vs promoted reference)",
                format_shape(*metrics.profile.observed_shape),
            ),
            (
                "profile score promoted reference shape",
                format_shape(*metrics.profile.expected_shape),
            ),
            (
                "profile score mean abs diff (rewrite vs promoted reference)",
                metrics.profile.mean_abs_diff,
            ),
            (
                "profile score mean Pearson correlation (rewrite vs promoted reference)",
                format_percent(metrics.profile.mean_pearson_corr),
            ),
            (
                "profile score mean Spearman correlation (rewrite vs promoted reference)",
                format_percent(metrics.profile.mean_spearman_corr),
            ),
            (
                "combined score shape (rewrite vs promoted reference)",
                format_shape(*metrics.combined.observed_shape),
            ),
            (
                "combined score promoted reference shape",
                format_shape(*metrics.combined.expected_shape),
            ),
            (
                "combined score mean abs diff (rewrite vs promoted reference)",
                metrics.combined.mean_abs_diff,
            ),
            (
                "combined score mean Pearson correlation (rewrite vs promoted reference)",
                format_percent(metrics.combined.mean_pearson_corr),
            ),
            (
                "combined score mean Spearman correlation (rewrite vs promoted reference)",
                format_percent(metrics.combined.mean_spearman_corr),
            ),
            (
                "weight table shape (rewrite vs promoted reference)",
                format_shape(*metrics.weights.observed_shape),
            ),
            (
                "weight table promoted reference shape",
                format_shape(*metrics.weights.expected_shape),
            ),
            (
                "weight table mean abs diff (rewrite vs promoted reference)",
                metrics.weights.mean_abs_diff,
            ),
            (
                "weight table max abs diff (rewrite vs promoted reference)",
                metrics.weights.max_abs_diff,
            ),
            (
                "contract: rewrite-vs-promoted-reference | prediction-matrix numeric surface",
                "asserted",
            ),
            (
                "prediction matrix shape (rewrite vs promoted reference)",
                format_shape(*metrics.prediction_matrix.observed_shape),
            ),
            (
                "prediction matrix promoted reference shape",
                format_shape(*metrics.prediction_matrix.expected_shape),
            ),
            (
                "prediction matrix mean abs diff (rewrite vs promoted reference)",
                metrics.prediction_matrix.mean_abs_diff,
            ),
            (
                "prediction matrix mean Pearson correlation (rewrite vs promoted reference)",
                format_percent(metrics.prediction_matrix.mean_pearson_corr),
            ),
            (
                "prediction matrix mean Spearman correlation (rewrite vs promoted reference)",
                format_percent(metrics.prediction_matrix.mean_spearman_corr),
            ),
            (
                "contract: rewrite-vs-promoted-reference | candidate-selection surface",
                "asserted",
            ),
            (
                "candidate rows (rewrite candidate-set surface)",
                metrics.candidates.observed_rows,
            ),
            (
                "candidate rows (promoted reference candidate-set surface)",
                metrics.candidates.expected_rows,
            ),
            (
                "candidate overlap precision (rewrite vs promoted reference candidate-set surface)",
                format_percent(metrics.candidates.overlap_precision),
            ),
            (
                "candidate overlap recall (rewrite vs promoted reference candidate-set surface)",
                format_percent(metrics.candidates.overlap_recall),
            ),
            (
                "candidate overlap F1 (rewrite vs promoted reference candidate-set surface)",
                format_percent(metrics.candidates.overlap_f1),
            ),
            (
                "contract: rewrite-vs-promoted-reference | prediction-matrix ranking surface",
                "asserted",
            ),
            (
                "prediction-matrix ranking mean Spearman (rewrite vs promoted reference)",
                format_percent(
                    metrics.prediction_matrix_ranking.mean_spearman_rank_corr
                ),
            ),
            (
                "prediction-matrix ranking top-10 overlap (rewrite vs promoted reference)",
                format_percent(metrics.prediction_matrix_ranking.mean_top10_overlap),
            ),
            (
                "prediction-matrix ranking top-20 overlap (rewrite vs promoted reference)",
                format_percent(metrics.prediction_matrix_ranking.mean_top20_overlap),
            ),
            (
                "prediction-matrix ranking top-30 overlap (rewrite vs promoted reference)",
                format_percent(metrics.prediction_matrix_ranking.mean_top30_overlap),
            ),
            (
                "prediction-matrix top-rank matches (rewrite vs promoted reference)",
                format_fraction(
                    metrics.prediction_matrix_ranking.top_rank_matches,
                    metrics.prediction_matrix_ranking.top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "contract: rewrite-vs-promoted-reference | ranked top-k export surface",
                "asserted",
            ),
            (
                "top-k export ranking mean Spearman (rewrite vs promoted reference)",
                format_percent(metrics.ranked_topk_export.mean_spearman_rank_corr),
            ),
            (
                "top-k export ranking top-10 overlap (rewrite vs promoted reference)",
                format_percent(metrics.ranked_topk_export.mean_top10_overlap),
            ),
            (
                "top-k export ranking top-20 overlap (rewrite vs promoted reference)",
                format_percent(metrics.ranked_topk_export.mean_top20_overlap),
            ),
            (
                "top-k export ranking top-30 overlap (rewrite vs promoted reference)",
                format_percent(metrics.ranked_topk_export.mean_top30_overlap),
            ),
            (
                "top-k export top-rank matches (rewrite vs promoted reference)",
                format_fraction(
                    metrics.ranked_topk_export.top_rank_matches,
                    metrics.ranked_topk_export.top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "contract: cross-policy divergence | prediction-matrix surface",
                "stable (default) vs r_parity",
            ),
            (
                "prediction matrix score correlation (cross-policy divergence: stable vs r_parity)",
                format_percent(divergence.prediction_matrix_score_corr),
            ),
            (
                "prediction matrix score mean abs diff (cross-policy divergence: stable vs r_parity)",
                divergence.prediction_matrix_score_mae,
            ),
            (
                "prediction-matrix ranking mean Spearman (cross-policy divergence: stable vs r_parity)",
                format_percent(
                    divergence.prediction_matrix_ranking.mean_spearman_rank_corr
                ),
            ),
            (
                "prediction-matrix ranking top-10 overlap (cross-policy divergence: stable vs r_parity)",
                format_percent(divergence.prediction_matrix_ranking.mean_top10_overlap),
            ),
            (
                "prediction-matrix ranking top-20 overlap (cross-policy divergence: stable vs r_parity)",
                format_percent(divergence.prediction_matrix_ranking.mean_top20_overlap),
            ),
            (
                "prediction-matrix ranking top-30 overlap (cross-policy divergence: stable vs r_parity)",
                format_percent(divergence.prediction_matrix_ranking.mean_top30_overlap),
            ),
            (
                "prediction-matrix top-rank matches (cross-policy divergence: stable vs r_parity)",
                format_fraction(
                    divergence.prediction_matrix_ranking.top_rank_matches,
                    divergence.prediction_matrix_ranking.top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "contract: cross-policy divergence | ranked top-k export surface",
                "stable (default) vs r_parity",
            ),
            (
                "top-k export ranking mean Spearman (cross-policy divergence: stable vs r_parity)",
                format_percent(divergence.ranked_topk_export.mean_spearman_rank_corr),
            ),
            (
                "top-k export ranking top-10 overlap (cross-policy divergence: stable vs r_parity)",
                format_percent(divergence.ranked_topk_export.mean_top10_overlap),
            ),
            (
                "top-k export ranking top-20 overlap (cross-policy divergence: stable vs r_parity)",
                format_percent(divergence.ranked_topk_export.mean_top20_overlap),
            ),
            (
                "top-k export ranking top-30 overlap (cross-policy divergence: stable vs r_parity)",
                format_percent(divergence.ranked_topk_export.mean_top30_overlap),
            ),
            (
                "top-k export top-rank matches (cross-policy divergence: stable vs r_parity)",
                format_fraction(
                    divergence.ranked_topk_export.top_rank_matches,
                    divergence.ranked_topk_export.top_rank_total,
                    include_percent=True,
                ),
            ),
        ],
        notes=(
            "fixture lane: tests/fixtures/rewrite_parity/r_reference_l6_prediction/",
            "rewrite-vs-promoted-reference contracts: prediction matrix, candidate set, and ranked top-k export are asserted independently",
            "cross-policy contract: stable (default) vs r_parity divergence is asserted independently from promoted-reference parity",
        ),
    )
