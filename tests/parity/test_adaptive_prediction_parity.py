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
CROSS_POLICY_DIVERGENCE_CORRELATION_FLOOR = 0.999
CROSS_POLICY_DIVERGENCE_MAE_CEILING = 0.08
CROSS_POLICY_DIVERGENCE_TOP_SET_OVERLAP_FLOOR = 7
CROSS_POLICY_DIVERGENCE_TOP10_OVERLAP_FLOOR = 0.99
CROSS_POLICY_DIVERGENCE_TOP20_OVERLAP_FLOOR = 0.99
CROSS_POLICY_DIVERGENCE_TOP30_OVERLAP_FLOOR = 0.99


def test_adaptive_prediction_surface_coverage_is_consistent_across_policies() -> None:
    metrics = collect_adaptive_policy_comparison_metrics()
    stable = metrics.stable
    r_parity = metrics.r_parity

    assert stable.candidate_count == stable.selected_trace_candidate_count
    assert r_parity.candidate_count == r_parity.selected_trace_candidate_count
    assert stable.candidate_count == r_parity.candidate_count
    assert stable.candidate_kinase_count == r_parity.candidate_kinase_count
    assert stable.prediction_shape == r_parity.prediction_shape
    assert stable.kinases_compared == r_parity.kinases_compared
    assert stable.donor_prediction_rows == r_parity.donor_prediction_rows


def test_adaptive_prediction_donor_parity_stable_policy() -> None:
    stable = collect_adaptive_policy_comparison_metrics().stable

    assert stable.donor_prediction_corr >= STABLE_DONOR_CORRELATION_FLOOR, (
        "stable donor-vs-rewrite parity regressed on prediction correlation"
    )
    assert stable.donor_prediction_mae <= STABLE_DONOR_MAE_CEILING, (
        "stable donor-vs-rewrite parity regressed on prediction MAE"
    )
    assert stable.donor_top_rank_matches >= STABLE_DONOR_TOP_RANK_MATCH_FLOOR
    assert stable.donor_mean_top10_overlap >= STABLE_DONOR_TOP10_OVERLAP_FLOOR
    assert stable.donor_mean_top20_overlap >= STABLE_DONOR_TOP20_OVERLAP_FLOOR
    assert stable.donor_mean_top30_overlap >= STABLE_DONOR_TOP30_OVERLAP_FLOOR


def test_adaptive_prediction_donor_parity_r_parity_policy() -> None:
    r_parity = collect_adaptive_policy_comparison_metrics().r_parity

    assert r_parity.donor_prediction_corr >= R_PARITY_DONOR_CORRELATION_FLOOR, (
        "r_parity donor-vs-rewrite parity regressed on prediction correlation"
    )
    assert r_parity.donor_prediction_mae <= R_PARITY_DONOR_MAE_CEILING, (
        "r_parity donor-vs-rewrite parity regressed on prediction MAE"
    )
    assert (
        r_parity.donor_top_set_overlap_matches >= R_PARITY_DONOR_TOP_SET_OVERLAP_FLOOR
    )


def test_adaptive_prediction_cross_policy_divergence_stable_vs_r_parity() -> None:
    metrics = collect_adaptive_policy_comparison_metrics()

    assert (
        metrics.cross_policy_prediction_corr
        >= CROSS_POLICY_DIVERGENCE_CORRELATION_FLOOR
    ), "cross-policy divergence regressed on prediction correlation"
    assert metrics.cross_policy_prediction_mae <= CROSS_POLICY_DIVERGENCE_MAE_CEILING
    assert (
        metrics.cross_policy_top_set_overlap_matches
        >= CROSS_POLICY_DIVERGENCE_TOP_SET_OVERLAP_FLOOR
    )
    assert metrics.cross_policy_top_rank_matches < metrics.cross_policy_top_rank_total
    assert (
        metrics.cross_policy_mean_top10_overlap
        >= CROSS_POLICY_DIVERGENCE_TOP10_OVERLAP_FLOOR
    )
    assert (
        metrics.cross_policy_mean_top20_overlap
        >= CROSS_POLICY_DIVERGENCE_TOP20_OVERLAP_FLOOR
    )
    assert (
        metrics.cross_policy_mean_top30_overlap
        >= CROSS_POLICY_DIVERGENCE_TOP30_OVERLAP_FLOOR
    )


def test_adaptive_prediction_parity_reporting_is_contract_explicit(
    request: pytest.FixtureRequest,
) -> None:
    metrics = collect_adaptive_policy_comparison_metrics()
    stable = metrics.stable
    r_parity = metrics.r_parity

    record_parity_metrics(
        request.config,
        family="adaptive_prediction",
        metrics=[
            ("adaptive_policy: stable (default)", "executed"),
            (
                "stable prediction shape (donor vs rewrite)",
                format_shape(*stable.prediction_shape),
            ),
            ("stable candidate count (donor vs rewrite)", stable.candidate_count),
            ("stable kinases compared (donor vs rewrite)", stable.kinases_compared),
            (
                "stable prediction correlation (donor vs rewrite)",
                format_percent(stable.donor_prediction_corr),
            ),
            (
                "stable prediction mean abs diff (donor vs rewrite)",
                stable.donor_prediction_mae,
            ),
            (
                "stable top-4 rank matches (donor vs rewrite)",
                format_fraction(
                    stable.donor_top_rank_matches,
                    stable.donor_top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "stable top-10 overlap (donor vs rewrite)",
                format_percent(stable.donor_mean_top10_overlap),
            ),
            (
                "stable top-20 overlap (donor vs rewrite)",
                format_percent(stable.donor_mean_top20_overlap),
            ),
            (
                "stable top-30 overlap (donor vs rewrite)",
                format_percent(stable.donor_mean_top30_overlap),
            ),
            ("adaptive_policy: r_parity", "executed"),
            (
                "r_parity prediction shape (donor vs rewrite)",
                format_shape(*r_parity.prediction_shape),
            ),
            ("r_parity candidate count (donor vs rewrite)", r_parity.candidate_count),
            ("r_parity kinases compared (donor vs rewrite)", r_parity.kinases_compared),
            (
                "r_parity prediction correlation (donor vs rewrite)",
                format_percent(r_parity.donor_prediction_corr),
            ),
            (
                "r_parity prediction mean abs diff (donor vs rewrite)",
                r_parity.donor_prediction_mae,
            ),
            (
                "r_parity top-4 rank matches (donor vs rewrite)",
                format_fraction(
                    r_parity.donor_top_rank_matches,
                    r_parity.donor_top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity top-4 set overlap (donor vs rewrite)",
                format_fraction(
                    r_parity.donor_top_set_overlap_matches,
                    r_parity.donor_top_set_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity top-10 overlap (donor vs rewrite)",
                format_percent(r_parity.donor_mean_top10_overlap),
            ),
            (
                "r_parity top-20 overlap (donor vs rewrite)",
                format_percent(r_parity.donor_mean_top20_overlap),
            ),
            (
                "r_parity top-30 overlap (donor vs rewrite)",
                format_percent(r_parity.donor_mean_top30_overlap),
            ),
            ("cross-policy divergence surface", "stable (default) vs r_parity"),
            (
                "cross-policy prediction correlation (stable vs r_parity)",
                format_percent(metrics.cross_policy_prediction_corr),
            ),
            (
                "cross-policy prediction mean abs diff (stable vs r_parity)",
                metrics.cross_policy_prediction_mae,
            ),
            (
                "cross-policy top-4 rank matches (stable vs r_parity)",
                format_fraction(
                    metrics.cross_policy_top_rank_matches,
                    metrics.cross_policy_top_rank_total,
                    include_percent=True,
                ),
            ),
            (
                "cross-policy top-4 set overlap (stable vs r_parity)",
                format_fraction(
                    metrics.cross_policy_top_set_overlap_matches,
                    metrics.cross_policy_top_set_overlap_total,
                    include_percent=True,
                ),
            ),
            (
                "cross-policy top-10 overlap (stable vs r_parity)",
                format_percent(metrics.cross_policy_mean_top10_overlap),
            ),
            (
                "cross-policy top-20 overlap (stable vs r_parity)",
                format_percent(metrics.cross_policy_mean_top20_overlap),
            ),
            (
                "cross-policy top-30 overlap (stable vs r_parity)",
                format_percent(metrics.cross_policy_mean_top30_overlap),
            ),
        ],
        notes=(
            "display mapping: adaptive_policy=stable corresponds to the default lane",
            "fixture lane: tests/fixtures/rewrite_parity/adaptive_sampling_edge/PROVENANCE.md",
            "donor-vs-rewrite contracts are reported separately from cross-policy divergence",
        ),
    )
