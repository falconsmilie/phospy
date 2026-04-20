from __future__ import annotations

import pytest

from tests.support.adaptive_parity_metrics import (
    collect_adaptive_policy_comparison_metrics,
)
from tests.support.parity_reporting import (
    format_fraction,
    format_percent,
    format_shape,
    record_parity_metrics,
)

pytestmark = pytest.mark.parity


STABLE_DONOR_CORRELATION_FLOOR = 0.999
STABLE_DONOR_MAE_CEILING = 0.01
STABLE_DONOR_TOP_RANK_MATCH_FLOOR = 7
STABLE_DONOR_TOP10_OVERLAP_FLOOR = 0.99
STABLE_DONOR_TOP20_OVERLAP_FLOOR = 0.99
STABLE_DONOR_TOP30_OVERLAP_FLOOR = 0.99
R_PARITY_DONOR_CORRELATION_FLOOR = 0.999
R_PARITY_DONOR_MAE_CEILING = 0.10
R_PARITY_DONOR_TOP_SET_OVERLAP_FLOOR = 5
MODE_COMPARISON_CORRELATION_FLOOR = 0.999
MODE_COMPARISON_MAE_CEILING = 0.08
MODE_COMPARISON_TOP_SET_OVERLAP_FLOOR = 7
MODE_COMPARISON_TOP10_OVERLAP_FLOOR = 0.99
MODE_COMPARISON_TOP20_OVERLAP_FLOOR = 0.99
MODE_COMPARISON_TOP30_OVERLAP_FLOOR = 0.99


def test_adaptive_ensemble_outputs_match_promoted_fixture_tolerances(
    request: pytest.FixtureRequest,
) -> None:
    metrics = collect_adaptive_policy_comparison_metrics()
    stable = metrics.stable
    r_parity = metrics.r_parity

    # This promoted fixture was captured from the default lane.
    # Keep strict donor-fit thresholds for stable, and bounded-comparison
    # thresholds for r_parity while preserving shared coverage across lanes.
    assert stable.candidate_count == stable.selected_trace_candidate_count
    assert r_parity.candidate_count == r_parity.selected_trace_candidate_count
    assert stable.candidate_count == r_parity.candidate_count
    assert stable.candidate_kinase_count == r_parity.candidate_kinase_count
    assert stable.prediction_shape == r_parity.prediction_shape
    assert stable.kinases_compared == r_parity.kinases_compared
    assert stable.donor_prediction_rows == r_parity.donor_prediction_rows

    assert stable.donor_prediction_corr >= STABLE_DONOR_CORRELATION_FLOOR
    assert stable.donor_prediction_mae <= STABLE_DONOR_MAE_CEILING
    assert stable.donor_top_rank_matches >= STABLE_DONOR_TOP_RANK_MATCH_FLOOR
    assert stable.donor_mean_top10_overlap >= STABLE_DONOR_TOP10_OVERLAP_FLOOR
    assert stable.donor_mean_top20_overlap >= STABLE_DONOR_TOP20_OVERLAP_FLOOR
    assert stable.donor_mean_top30_overlap >= STABLE_DONOR_TOP30_OVERLAP_FLOOR

    assert r_parity.donor_prediction_corr >= R_PARITY_DONOR_CORRELATION_FLOOR
    assert r_parity.donor_prediction_mae <= R_PARITY_DONOR_MAE_CEILING
    assert (
        r_parity.donor_top_set_overlap_matches >= R_PARITY_DONOR_TOP_SET_OVERLAP_FLOOR
    )

    assert metrics.cross_policy_prediction_corr >= MODE_COMPARISON_CORRELATION_FLOOR
    assert metrics.cross_policy_prediction_mae <= MODE_COMPARISON_MAE_CEILING
    assert (
        metrics.cross_policy_top_set_overlap_matches
        >= MODE_COMPARISON_TOP_SET_OVERLAP_FLOOR
    )
    assert metrics.cross_policy_top_rank_matches < metrics.cross_policy_top_rank_total
    assert (
        metrics.cross_policy_mean_top10_overlap >= MODE_COMPARISON_TOP10_OVERLAP_FLOOR
    )
    assert (
        metrics.cross_policy_mean_top20_overlap >= MODE_COMPARISON_TOP20_OVERLAP_FLOOR
    )
    assert (
        metrics.cross_policy_mean_top30_overlap >= MODE_COMPARISON_TOP30_OVERLAP_FLOOR
    )

    assert r_parity.donor_mean_top10_overlap >= stable.donor_mean_top10_overlap
    assert r_parity.donor_mean_top20_overlap >= stable.donor_mean_top20_overlap
    assert r_parity.donor_mean_top30_overlap >= stable.donor_mean_top30_overlap

    record_parity_metrics(
        request.config,
        family="adaptive_prediction",
        metrics=[
            ("adaptive_policy: stable (default)", "executed"),
            ("stable (default) candidate count", stable.candidate_count),
            (
                "stable (default) prediction shape",
                format_shape(*stable.prediction_shape),
            ),
            ("stable (default) kinases compared", stable.kinases_compared),
            (
                "stable (default) donor correlation",
                format_percent(stable.donor_prediction_corr),
            ),
            (
                "stable (default) donor mean abs diff",
                stable.donor_prediction_mae,
            ),
            (
                "stable (default) donor top-4 rank matches",
                format_fraction(
                    stable.donor_top_rank_matches,
                    stable.donor_top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "stable (default) donor top-10 overlap",
                format_percent(stable.donor_mean_top10_overlap),
            ),
            (
                "stable (default) donor top-20 overlap",
                format_percent(stable.donor_mean_top20_overlap),
            ),
            (
                "stable (default) donor top-30 overlap",
                format_percent(stable.donor_mean_top30_overlap),
            ),
            ("adaptive_policy: r_parity", "executed"),
            ("r_parity candidate count", r_parity.candidate_count),
            ("r_parity prediction shape", format_shape(*r_parity.prediction_shape)),
            ("r_parity kinases compared", r_parity.kinases_compared),
            (
                "r_parity donor correlation",
                format_percent(r_parity.donor_prediction_corr),
            ),
            ("r_parity donor mean abs diff", r_parity.donor_prediction_mae),
            (
                "r_parity donor top-4 rank matches",
                format_fraction(
                    r_parity.donor_top_rank_matches,
                    r_parity.donor_top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity donor top-4 set overlap",
                format_fraction(
                    r_parity.donor_top_set_overlap_matches,
                    r_parity.donor_top_set_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity donor top-10 overlap",
                format_percent(r_parity.donor_mean_top10_overlap),
            ),
            (
                "r_parity donor top-20 overlap",
                format_percent(r_parity.donor_mean_top20_overlap),
            ),
            (
                "r_parity donor top-30 overlap",
                format_percent(r_parity.donor_mean_top30_overlap),
            ),
            (
                "adaptive policy mode comparison",
                "stable (default) vs r_parity",
            ),
            (
                "mode comparison donor correlation delta (r_parity - stable)",
                r_parity.donor_prediction_corr - stable.donor_prediction_corr,
            ),
            (
                "mode comparison donor MAE delta (r_parity - stable)",
                r_parity.donor_prediction_mae - stable.donor_prediction_mae,
            ),
            (
                "mode comparison prediction correlation (stable vs r_parity)",
                format_percent(metrics.cross_policy_prediction_corr),
            ),
            (
                "mode comparison prediction mean abs diff (stable vs r_parity)",
                metrics.cross_policy_prediction_mae,
            ),
            (
                "mode comparison top-4 rank matches (stable vs r_parity)",
                format_fraction(
                    metrics.cross_policy_top_rank_matches,
                    metrics.cross_policy_top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "mode comparison top-4 set overlap (stable vs r_parity)",
                format_fraction(
                    metrics.cross_policy_top_set_overlap_matches,
                    metrics.cross_policy_top_set_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "mode comparison top-10 overlap (stable vs r_parity)",
                format_percent(metrics.cross_policy_mean_top10_overlap),
            ),
            (
                "mode comparison top-20 overlap (stable vs r_parity)",
                format_percent(metrics.cross_policy_mean_top20_overlap),
            ),
            (
                "mode comparison top-30 overlap (stable vs r_parity)",
                format_percent(metrics.cross_policy_mean_top30_overlap),
            ),
        ],
        notes=(
            "display mapping: adaptive_policy=stable corresponds to the default lane",
            "fixture lane: tests/fixtures/rewrite_parity/adaptive_sampling_edge/PROVENANCE.md",
        ),
    )
