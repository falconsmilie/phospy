from __future__ import annotations

import pytest

from tests.support.adaptive_trace_parity_metrics import (
    collect_adaptive_trace_replay_metrics,
)
from tests.support.parity_reporting import (
    format_fraction,
    format_percent,
    record_parity_metrics,
)

pytestmark = pytest.mark.parity

INITIAL_OVERLAP_RATIO_FLOOR = 0.05
SAMPLE_OVERLAP_RATIO_FLOOR = 0.30
STABLE_DONOR_CORR_FLOOR = 0.84
STABLE_DONOR_MAE_CEILING = 0.12
R_PARITY_DONOR_CORR_FLOOR = 0.80
R_PARITY_DONOR_MAE_CEILING = 0.13
DONOR_TOP_SET_OVERLAP_FLOOR = 0.75
DONOR_TOP10_OVERLAP_FLOOR = 0.70
DONOR_TOP20_OVERLAP_FLOOR = 0.80
DONOR_TOP30_OVERLAP_FLOOR = 0.80
CROSS_POLICY_CORRELATION_FLOOR = 0.97
CROSS_POLICY_MAE_CEILING = 0.05
CROSS_POLICY_TOP_SET_OVERLAP_FLOOR = 0.60


def _ratio(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return 0.0
    return float(numerator / denominator)


def test_adaptive_replay_trace_parity_matches_promoted_trace_surfaces(
    request: pytest.FixtureRequest,
) -> None:
    metrics = collect_adaptive_trace_replay_metrics()
    stable = metrics.stable
    r_parity = metrics.r_parity

    stable_initial_ratio = _ratio(
        stable.initial_overlap_matches, stable.initial_overlap_total
    )
    stable_sample_ratio = _ratio(
        stable.sample_overlap_matches, stable.sample_overlap_total
    )
    r_parity_initial_ratio = _ratio(
        r_parity.initial_overlap_matches,
        r_parity.initial_overlap_total,
    )
    r_parity_sample_ratio = _ratio(
        r_parity.sample_overlap_matches,
        r_parity.sample_overlap_total,
    )
    stable_top_set_ratio = _ratio(
        stable.donor_top_set_overlap_matches,
        stable.donor_top_set_overlap_total,
    )
    r_parity_top_set_ratio = _ratio(
        r_parity.donor_top_set_overlap_matches,
        r_parity.donor_top_set_overlap_total,
    )
    cross_top_set_ratio = _ratio(
        metrics.cross_policy_top_set_overlap_matches,
        metrics.cross_policy_top_set_overlap_total,
    )

    assert stable.candidate_count == r_parity.candidate_count
    assert stable.candidate_kinase_count == r_parity.candidate_kinase_count
    assert stable.selected_trace_kinase_count == r_parity.selected_trace_kinase_count

    assert stable_initial_ratio >= INITIAL_OVERLAP_RATIO_FLOOR
    assert stable_sample_ratio >= SAMPLE_OVERLAP_RATIO_FLOOR
    assert r_parity_initial_ratio >= INITIAL_OVERLAP_RATIO_FLOOR
    assert r_parity_sample_ratio >= SAMPLE_OVERLAP_RATIO_FLOOR

    assert stable.donor_final_prediction_corr >= STABLE_DONOR_CORR_FLOOR
    assert stable.donor_final_prediction_mae <= STABLE_DONOR_MAE_CEILING
    assert r_parity.donor_final_prediction_corr >= R_PARITY_DONOR_CORR_FLOOR
    assert r_parity.donor_final_prediction_mae <= R_PARITY_DONOR_MAE_CEILING

    assert stable_top_set_ratio >= DONOR_TOP_SET_OVERLAP_FLOOR
    assert stable.donor_mean_top10_overlap >= DONOR_TOP10_OVERLAP_FLOOR
    assert stable.donor_mean_top20_overlap >= DONOR_TOP20_OVERLAP_FLOOR
    assert stable.donor_mean_top30_overlap >= DONOR_TOP30_OVERLAP_FLOOR

    assert r_parity_top_set_ratio >= DONOR_TOP_SET_OVERLAP_FLOOR
    assert r_parity.donor_mean_top10_overlap >= DONOR_TOP10_OVERLAP_FLOOR
    assert r_parity.donor_mean_top20_overlap >= DONOR_TOP20_OVERLAP_FLOOR
    assert r_parity.donor_mean_top30_overlap >= DONOR_TOP30_OVERLAP_FLOOR

    assert stable.deterministic_under_seed
    assert r_parity.deterministic_under_seed

    assert metrics.cross_policy_prediction_corr >= CROSS_POLICY_CORRELATION_FLOOR
    assert metrics.cross_policy_prediction_mae <= CROSS_POLICY_MAE_CEILING
    assert cross_top_set_ratio >= CROSS_POLICY_TOP_SET_OVERLAP_FLOOR

    record_parity_metrics(
        request.config,
        family="adaptive_replay",
        metrics=[
            ("adaptive_policy: stable (default)", "executed"),
            ("stable candidate count", stable.candidate_count),
            ("stable trace kinases", stable.selected_trace_kinase_count),
            (
                "stable initial replay overlap",
                format_fraction(
                    stable.initial_overlap_matches,
                    stable.initial_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "stable iteration-sample overlap",
                format_fraction(
                    stable.sample_overlap_matches,
                    stable.sample_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "stable donor final correlation",
                format_percent(stable.donor_final_prediction_corr),
            ),
            ("stable donor final mean abs diff", stable.donor_final_prediction_mae),
            (
                "stable donor top-rank matches",
                format_fraction(
                    stable.donor_top_rank_matches,
                    stable.donor_top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "stable donor top-set overlap",
                format_fraction(
                    stable.donor_top_set_overlap_matches,
                    stable.donor_top_set_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "stable donor top-10 overlap",
                format_percent(stable.donor_mean_top10_overlap),
            ),
            (
                "stable donor top-20 overlap",
                format_percent(stable.donor_mean_top20_overlap),
            ),
            (
                "stable donor top-30 overlap",
                format_percent(stable.donor_mean_top30_overlap),
            ),
            (
                "stable replay deterministic under fixed seed",
                stable.deterministic_under_seed,
            ),
            ("adaptive_policy: r_parity", "executed"),
            ("r_parity candidate count", r_parity.candidate_count),
            ("r_parity trace kinases", r_parity.selected_trace_kinase_count),
            (
                "r_parity initial replay overlap",
                format_fraction(
                    r_parity.initial_overlap_matches,
                    r_parity.initial_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity iteration-sample overlap",
                format_fraction(
                    r_parity.sample_overlap_matches,
                    r_parity.sample_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity donor final correlation",
                format_percent(r_parity.donor_final_prediction_corr),
            ),
            ("r_parity donor final mean abs diff", r_parity.donor_final_prediction_mae),
            (
                "r_parity donor top-rank matches",
                format_fraction(
                    r_parity.donor_top_rank_matches,
                    r_parity.donor_top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity donor top-set overlap",
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
                "r_parity replay deterministic under fixed seed",
                r_parity.deterministic_under_seed,
            ),
            (
                "stable vs r_parity replay prediction correlation",
                format_percent(metrics.cross_policy_prediction_corr),
            ),
            (
                "stable vs r_parity replay prediction mean abs diff",
                metrics.cross_policy_prediction_mae,
            ),
            (
                "stable vs r_parity replay top-set overlap",
                format_fraction(
                    metrics.cross_policy_top_set_overlap_matches,
                    metrics.cross_policy_top_set_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "stable vs r_parity replay top-10 overlap",
                format_percent(metrics.cross_policy_mean_top10_overlap),
            ),
            (
                "stable vs r_parity replay top-20 overlap",
                format_percent(metrics.cross_policy_mean_top20_overlap),
            ),
            (
                "stable vs r_parity replay top-30 overlap",
                format_percent(metrics.cross_policy_mean_top30_overlap),
            ),
        ],
        notes=(
            "fixture lane: tests/fixtures/rewrite_parity/adaptive_sampling_replay/",
            "comparison labels: adaptive_policy=stable (default), adaptive_policy=r_parity",
        ),
    )
