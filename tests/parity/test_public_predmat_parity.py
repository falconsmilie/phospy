from __future__ import annotations

import pytest

from tests.support.parity_reporting import (
    format_fraction,
    format_percent,
    format_shape,
    record_parity_metrics,
)
from tests.support.public_predmat_parity_metrics import (
    collect_public_predmat_benchmark_metrics,
    collect_public_predmat_order_invariance_metrics,
    load_public_predmat_contract,
)

pytestmark = pytest.mark.parity

CROSS_POLICY_CORRELATION_FLOOR = 0.80
CROSS_POLICY_MAE_CEILING = 0.12


def test_public_predmat_workflow_matches_rewrite_committed_benchmarks(
    request: pytest.FixtureRequest,
) -> None:
    contract = load_public_predmat_contract()
    metrics = collect_public_predmat_benchmark_metrics()
    stable = metrics.stable
    r_parity = metrics.r_parity

    for lane in (stable, r_parity):
        assert lane.observed_shape == lane.expected_shape
        assert lane.row_identity_match
        assert lane.column_identity_match
        assert lane.mean_abs_diff == pytest.approx(0.0, abs=1e-12)
        assert lane.max_abs_diff == pytest.approx(0.0, abs=1e-12)
        assert lane.dominant_matches == lane.dominant_total
        assert lane.deterministic_under_seed
        assert lane.donor_corr >= 0.99
        assert lane.donor_mae <= 0.30

    assert metrics.cross_policy_corr >= CROSS_POLICY_CORRELATION_FLOOR
    assert metrics.cross_policy_mae <= CROSS_POLICY_MAE_CEILING

    stable_contract = contract["benchmarks"]["stable"]
    r_parity_contract = contract["benchmarks"]["r_parity"]
    assert stable.dominant_total == len(stable_contract["dominant_kinase_by_site"])
    assert r_parity.dominant_total == len(r_parity_contract["dominant_kinase_by_site"])

    record_parity_metrics(
        request.config,
        family="public_predmat",
        metrics=[
            ("adaptive_policy: stable (default)", "executed"),
            ("stable predMat shape", format_shape(*stable.observed_shape)),
            ("stable row identity match", stable.row_identity_match),
            ("stable column identity match", stable.column_identity_match),
            (
                "stable dominant kinase matches",
                format_fraction(
                    stable.dominant_matches, stable.dominant_total, include_percent=True
                ),
            ),
            ("stable deterministic under fixed seed", stable.deterministic_under_seed),
            ("stable donor correlation", format_percent(stable.donor_corr)),
            ("stable donor mean abs diff", stable.donor_mae),
            ("adaptive_policy: r_parity", "executed"),
            ("r_parity predMat shape", format_shape(*r_parity.observed_shape)),
            ("r_parity row identity match", r_parity.row_identity_match),
            ("r_parity column identity match", r_parity.column_identity_match),
            (
                "r_parity dominant kinase matches",
                format_fraction(
                    r_parity.dominant_matches,
                    r_parity.dominant_total,
                    include_percent=True,
                ),
            ),
            (
                "r_parity deterministic under fixed seed",
                r_parity.deterministic_under_seed,
            ),
            ("r_parity donor correlation", format_percent(r_parity.donor_corr)),
            ("r_parity donor mean abs diff", r_parity.donor_mae),
            (
                "stable vs r_parity prediction correlation",
                format_percent(metrics.cross_policy_corr),
            ),
            ("stable vs r_parity prediction mean abs diff", metrics.cross_policy_mae),
        ],
        notes=(
            "fixture lane: tests/fixtures/public_workflow_reference/predmat_rewrite_*.csv",
            "rewrite public contract uses adaptive_policy (svm_mode is archival naming)",
        ),
    )


def test_public_predmat_workflow_stable_lane_is_order_invariant_end_to_end(
    request: pytest.FixtureRequest,
) -> None:
    metrics = collect_public_predmat_order_invariance_metrics()

    assert metrics.normalized_equal
    assert metrics.deterministic_under_seed
    assert metrics.dominant_matches == metrics.dominant_total

    record_parity_metrics(
        request.config,
        family="order_invariance",
        metrics=[
            ("lane", "public predMat stable (default)"),
            ("normalized predMat equality", metrics.normalized_equal),
            ("predMat output shape", format_shape(*metrics.output_shape)),
            (
                "dominant kinase matches",
                format_fraction(
                    metrics.dominant_matches,
                    metrics.dominant_total,
                    include_percent=True,
                ),
            ),
            ("deterministic under fixed seed", metrics.deterministic_under_seed),
        ],
        notes=(
            "order perturbations: reference substrate map order + site-sequence map order",
            "assertion uses row/column normalized equality",
        ),
    )
