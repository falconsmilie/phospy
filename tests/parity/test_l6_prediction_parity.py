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
CANDIDATE_OVERLAP_PRECISION_FLOOR = 0.999999
CANDIDATE_OVERLAP_RECALL_FLOOR = 0.999999
DEFAULT_PRED_MAT_RANK_SPEARMAN_FLOOR = 0.96
DEFAULT_PRED_MAT_TOP20_OVERLAP_FLOOR = 0.75
DEFAULT_PRED_MAT_TOP30_OVERLAP_FLOOR = 0.65
DEFAULT_PRED_MAT_GOOD_TOP10_COUNT_FLOOR = 20
R_PARITY_PRED_MAT_RANK_SPEARMAN_FLOOR = 0.94
R_PARITY_PRED_MAT_TOP10_OVERLAP_FLOOR = 0.74
R_PARITY_PRED_MAT_TOP20_OVERLAP_FLOOR = 0.72
R_PARITY_PRED_MAT_TOP30_OVERLAP_FLOOR = 0.65
R_PARITY_PRED_MAT_TOP_RANK_MATCH_FLOOR = 20
DEFAULT_TOPK_EXPORT_TOP20_OVERLAP_FLOOR = 0.75
DEFAULT_TOPK_EXPORT_TOP30_OVERLAP_FLOOR = 0.65
DEFAULT_TOPK_EXPORT_GOOD_TOP10_COUNT_FLOOR = 20
R_PARITY_TOPK_EXPORT_TOP10_OVERLAP_FLOOR = 0.74
R_PARITY_TOPK_EXPORT_TOP20_OVERLAP_FLOOR = 0.72
R_PARITY_TOPK_EXPORT_TOP30_OVERLAP_FLOOR = 0.65
R_PARITY_TOPK_EXPORT_TOP_RANK_MATCH_FLOOR = 20
CROSS_POLICY_CORRELATION_FLOOR = 0.95


def test_l6_full_prediction_and_scoring_parity_against_promoted_reference_tables(
    request: pytest.FixtureRequest,
) -> None:
    metrics = collect_l6_prediction_parity_metrics()

    assert metrics.profile.shared_site_count == metrics.profile.observed_shape[0]
    assert metrics.profile.expected_shape[0] == metrics.profile.observed_shape[0]
    assert metrics.profile.shared_kinase_count == metrics.profile.observed_shape[1]
    assert metrics.profile.expected_shape[1] >= metrics.profile.observed_shape[1]
    assert metrics.profile.mean_abs_diff <= PROFILE_MAE_CEILING
    assert metrics.profile.mean_pearson_corr >= PROFILE_PEARSON_FLOOR
    assert metrics.profile.mean_spearman_corr >= PROFILE_SPEARMAN_FLOOR

    assert metrics.combined.shared_site_count == metrics.combined.observed_shape[0]
    assert metrics.combined.expected_shape[0] == metrics.combined.observed_shape[0]
    assert metrics.combined.shared_kinase_count >= 15
    assert metrics.combined.shared_kinase_count <= metrics.combined.observed_shape[1]
    assert metrics.combined.expected_shape[1] >= metrics.combined.observed_shape[1]
    assert metrics.combined.mean_abs_diff <= COMBINED_MAE_CEILING
    assert metrics.combined.mean_pearson_corr >= COMBINED_PEARSON_FLOOR
    assert metrics.combined.mean_spearman_corr >= COMBINED_SPEARMAN_FLOOR

    assert metrics.weights.shared_site_count >= 15
    assert metrics.weights.shared_site_count <= metrics.weights.observed_shape[0]
    assert metrics.weights.shared_kinase_count == metrics.weights.observed_shape[1]
    assert metrics.weights.mean_abs_diff <= WEIGHTS_MAE_CEILING
    assert metrics.weights.max_abs_diff <= WEIGHTS_MAX_ABS_DIFF_CEILING

    assert metrics.candidates.observed_rows > 0
    assert metrics.candidates.expected_rows > 0
    assert metrics.candidates.overlap_rows > 0
    assert metrics.candidates.overlap_precision >= CANDIDATE_OVERLAP_PRECISION_FLOOR
    assert metrics.candidates.overlap_recall >= CANDIDATE_OVERLAP_RECALL_FLOOR
    assert metrics.candidates.shared_kinase_count > 0

    stable_pred_mat = metrics.stable_pred_mat_ranking
    r_parity_pred_mat = metrics.r_parity_pred_mat_ranking
    stable_topk_export = metrics.stable_topk_export_ranking
    r_parity_topk_export = metrics.r_parity_topk_export_ranking
    assert stable_pred_mat.kinases_compared > 0
    assert r_parity_pred_mat.kinases_compared > 0
    assert stable_topk_export.kinases_compared > 0
    assert r_parity_topk_export.kinases_compared > 0
    assert stable_pred_mat.kinases_compared == r_parity_pred_mat.kinases_compared
    assert stable_topk_export.kinases_compared == r_parity_topk_export.kinases_compared

    assert (
        stable_pred_mat.mean_spearman_rank_corr >= DEFAULT_PRED_MAT_RANK_SPEARMAN_FLOOR
    )
    assert stable_pred_mat.mean_top20_overlap >= DEFAULT_PRED_MAT_TOP20_OVERLAP_FLOOR
    assert stable_pred_mat.mean_top30_overlap >= DEFAULT_PRED_MAT_TOP30_OVERLAP_FLOOR
    assert stable_pred_mat.good_top10_count >= DEFAULT_PRED_MAT_GOOD_TOP10_COUNT_FLOOR
    assert stable_pred_mat.top_rank_total == stable_pred_mat.kinases_compared

    assert (
        r_parity_pred_mat.mean_spearman_rank_corr
        >= R_PARITY_PRED_MAT_RANK_SPEARMAN_FLOOR
    )
    assert r_parity_pred_mat.mean_top10_overlap >= R_PARITY_PRED_MAT_TOP10_OVERLAP_FLOOR
    assert r_parity_pred_mat.mean_top20_overlap >= R_PARITY_PRED_MAT_TOP20_OVERLAP_FLOOR
    assert r_parity_pred_mat.mean_top30_overlap >= R_PARITY_PRED_MAT_TOP30_OVERLAP_FLOOR
    assert r_parity_pred_mat.top_rank_matches >= R_PARITY_PRED_MAT_TOP_RANK_MATCH_FLOOR
    assert r_parity_pred_mat.top_rank_total == r_parity_pred_mat.kinases_compared

    assert (
        stable_topk_export.mean_top20_overlap >= DEFAULT_TOPK_EXPORT_TOP20_OVERLAP_FLOOR
    )
    assert (
        stable_topk_export.mean_top30_overlap >= DEFAULT_TOPK_EXPORT_TOP30_OVERLAP_FLOOR
    )
    assert (
        stable_topk_export.good_top10_count
        >= DEFAULT_TOPK_EXPORT_GOOD_TOP10_COUNT_FLOOR
    )
    assert stable_topk_export.top_rank_total == stable_topk_export.kinases_compared

    assert (
        r_parity_topk_export.mean_top10_overlap
        >= R_PARITY_TOPK_EXPORT_TOP10_OVERLAP_FLOOR
    )
    assert (
        r_parity_topk_export.mean_top20_overlap
        >= R_PARITY_TOPK_EXPORT_TOP20_OVERLAP_FLOOR
    )
    assert (
        r_parity_topk_export.mean_top30_overlap
        >= R_PARITY_TOPK_EXPORT_TOP30_OVERLAP_FLOOR
    )
    assert (
        r_parity_topk_export.top_rank_matches
        >= R_PARITY_TOPK_EXPORT_TOP_RANK_MATCH_FLOOR
    )
    assert r_parity_topk_export.top_rank_total == r_parity_topk_export.kinases_compared

    assert metrics.cross_policy_prediction_corr >= CROSS_POLICY_CORRELATION_FLOOR

    record_parity_metrics(
        request.config,
        family="l6_prediction",
        metrics=[
            ("profile score shape", format_shape(*metrics.profile.observed_shape)),
            (
                "profile score donor shape",
                format_shape(*metrics.profile.expected_shape),
            ),
            ("profile score mean abs diff", metrics.profile.mean_abs_diff),
            (
                "profile score mean Pearson correlation",
                format_percent(metrics.profile.mean_pearson_corr),
            ),
            (
                "profile score mean Spearman correlation",
                format_percent(metrics.profile.mean_spearman_corr),
            ),
            ("combined score shape", format_shape(*metrics.combined.observed_shape)),
            (
                "combined score donor shape",
                format_shape(*metrics.combined.expected_shape),
            ),
            ("combined score mean abs diff", metrics.combined.mean_abs_diff),
            (
                "combined score mean Pearson correlation",
                format_percent(metrics.combined.mean_pearson_corr),
            ),
            (
                "combined score mean Spearman correlation",
                format_percent(metrics.combined.mean_spearman_corr),
            ),
            ("weight table shape", format_shape(*metrics.weights.observed_shape)),
            ("weight table donor shape", format_shape(*metrics.weights.expected_shape)),
            ("weight table mean abs diff", metrics.weights.mean_abs_diff),
            ("weight table max abs diff", metrics.weights.max_abs_diff),
            (
                "candidate rows (rewrite, candidate-set surface)",
                metrics.candidates.observed_rows,
            ),
            (
                "candidate rows (donor, candidate-set surface)",
                metrics.candidates.expected_rows,
            ),
            (
                "candidate overlap precision (candidate-set surface)",
                format_percent(metrics.candidates.overlap_precision),
            ),
            (
                "candidate overlap recall (candidate-set surface)",
                format_percent(metrics.candidates.overlap_recall),
            ),
            (
                "candidate overlap F1 (candidate-set surface)",
                format_percent(metrics.candidates.overlap_f1),
            ),
            (
                "stable rank mean Spearman (prediction matrix surface)",
                format_percent(metrics.stable_pred_mat_ranking.mean_spearman_rank_corr),
            ),
            (
                "stable rank top-10 overlap (prediction matrix surface)",
                format_percent(metrics.stable_pred_mat_ranking.mean_top10_overlap),
            ),
            (
                "stable rank top-20 overlap (prediction matrix surface)",
                format_percent(metrics.stable_pred_mat_ranking.mean_top20_overlap),
            ),
            (
                "stable rank top-30 overlap (prediction matrix surface)",
                format_percent(metrics.stable_pred_mat_ranking.mean_top30_overlap),
            ),
            (
                "stable top-rank matches (prediction matrix surface)",
                format_fraction(
                    metrics.stable_pred_mat_ranking.top_rank_matches,
                    metrics.stable_pred_mat_ranking.top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity rank mean Spearman (prediction matrix surface)",
                format_percent(
                    metrics.r_parity_pred_mat_ranking.mean_spearman_rank_corr
                ),
            ),
            (
                "r_parity rank top-10 overlap (prediction matrix surface)",
                format_percent(metrics.r_parity_pred_mat_ranking.mean_top10_overlap),
            ),
            (
                "r_parity rank top-20 overlap (prediction matrix surface)",
                format_percent(metrics.r_parity_pred_mat_ranking.mean_top20_overlap),
            ),
            (
                "r_parity rank top-30 overlap (prediction matrix surface)",
                format_percent(metrics.r_parity_pred_mat_ranking.mean_top30_overlap),
            ),
            (
                "r_parity top-rank matches (prediction matrix surface)",
                format_fraction(
                    metrics.r_parity_pred_mat_ranking.top_rank_matches,
                    metrics.r_parity_pred_mat_ranking.top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "stable rank mean Spearman (top-k export surface)",
                format_percent(
                    metrics.stable_topk_export_ranking.mean_spearman_rank_corr
                ),
            ),
            (
                "stable rank top-10 overlap (top-k export surface)",
                format_percent(metrics.stable_topk_export_ranking.mean_top10_overlap),
            ),
            (
                "stable rank top-20 overlap (top-k export surface)",
                format_percent(metrics.stable_topk_export_ranking.mean_top20_overlap),
            ),
            (
                "stable rank top-30 overlap (top-k export surface)",
                format_percent(metrics.stable_topk_export_ranking.mean_top30_overlap),
            ),
            (
                "stable top-rank matches (top-k export surface)",
                format_fraction(
                    metrics.stable_topk_export_ranking.top_rank_matches,
                    metrics.stable_topk_export_ranking.top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity rank mean Spearman (top-k export surface)",
                format_percent(
                    metrics.r_parity_topk_export_ranking.mean_spearman_rank_corr
                ),
            ),
            (
                "r_parity rank top-10 overlap (top-k export surface)",
                format_percent(metrics.r_parity_topk_export_ranking.mean_top10_overlap),
            ),
            (
                "r_parity rank top-20 overlap (top-k export surface)",
                format_percent(metrics.r_parity_topk_export_ranking.mean_top20_overlap),
            ),
            (
                "r_parity rank top-30 overlap (top-k export surface)",
                format_percent(metrics.r_parity_topk_export_ranking.mean_top30_overlap),
            ),
            (
                "r_parity top-rank matches (top-k export surface)",
                format_fraction(
                    metrics.r_parity_topk_export_ranking.top_rank_matches,
                    metrics.r_parity_topk_export_ranking.top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "stable vs r_parity prediction correlation",
                format_percent(metrics.cross_policy_prediction_corr),
            ),
            (
                "stable vs r_parity prediction mean abs diff",
                metrics.cross_policy_prediction_mae,
            ),
            (
                "stable vs r_parity top-10 overlap",
                format_percent(metrics.cross_policy_mean_top10_overlap),
            ),
            (
                "stable vs r_parity top-20 overlap",
                format_percent(metrics.cross_policy_mean_top20_overlap),
            ),
            (
                "stable vs r_parity top-30 overlap",
                format_percent(metrics.cross_policy_mean_top30_overlap),
            ),
        ],
        notes=(
            "fixture lane: tests/fixtures/rewrite_parity/r_reference_l6_prediction/",
            "policy labels: adaptive_policy=stable (default lane), adaptive_policy=r_parity",
            "prediction-matrix ranking surface: pred_mat (rewrite) vs predMat.csv (fixture)",
            (
                "top-k export ranking surface: substrate_list (rewrite) vs "
                "native_prediction_top30.csv (fixture)"
            ),
            (
                "candidate-set surface: substrate_list kinase/site pairs vs "
                "native_candidate_substrates.csv"
            ),
        ),
    )
