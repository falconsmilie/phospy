from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pandas as pd
import pytest

import phospy.science.differential.duplicate_correlation as duplicate_correlation_module
from phospy.errors import PhosPyInputError
from phospy.science.differential.duplicate_correlation import (
    _fisher_trimmed_mean_consensus,
    estimate_duplicate_correlation_reml_consensus,
)
from phospy.science.differential.models.duplicate_correlation import (
    DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
    DUPLICATE_CORRELATION_TRIM_FRACTION,
    DuplicateCorrelationBoundary,
    DuplicateCorrelationConsensusResult,
    DuplicateCorrelationFailureReason,
    DuplicateCorrelationFeatureEstimate,
    DuplicateCorrelationFeatureStatus,
    DuplicateCorrelationReasonCount,
)

ROOT = Path(__file__).resolve().parents[4]
FIXTURE_ROOT = (
    ROOT
    / "tests"
    / "fixtures"
    / "rewrite_parity"
    / "differential_duplicate_correlation"
)

# Central tolerance for the residual-space two-variance-component REML route
# used by limma/statmod black-box fixtures. Most rows compare at floating-point
# noise; the small envelope covers QR/SVD basis differences on near-constant
# rows while still detecting final p-value-relevant estimator drift.
REML_REFERENCE_CORRELATION_ABSOLUTE_TOLERANCE = 1.0e-6
REML_REFERENCE_FISHER_Z_ABSOLUTE_TOLERANCE = 1.0e-6
STRICT_CORRELATION_ABSOLUTE_TOLERANCE = 1.0e-8

FIXTURE_IDS = (
    "fixture_a_complete_pairs",
    "fixture_b_three_observation_blocks",
    "fixture_c_incomplete_unequal_blocks",
    "fixture_d_feature_level_failures",
)

SUCCESS_STATUSES = {
    DuplicateCorrelationFeatureStatus.ESTIMATED,
    DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED,
}


@pytest.mark.parametrize("fixture_id", FIXTURE_IDS)
def test_reml_feature_correlations_match_limma_reference_fixtures(
    fixture_id: str,
) -> None:
    matrix, design, blocks, feature_ids = _fixture_inputs(fixture_id)
    expected = pd.read_csv(
        FIXTURE_ROOT / fixture_id / "feature_correlations.csv",
        keep_default_na=False,
    )

    result = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=feature_ids,
    )

    estimates = {estimate.feature_id: estimate for estimate in result.feature_estimates}
    for expected_row in expected.itertuples(index=False):
        estimate = estimates[str(expected_row.feature_id)]
        if expected_row.status == "estimated":
            assert estimate.status in SUCCESS_STATUSES
            assert estimate.correlation is not None
            assert math.isclose(
                estimate.correlation,
                float(expected_row.correlation),
                rel_tol=0.0,
                abs_tol=REML_REFERENCE_CORRELATION_ABSOLUTE_TOLERANCE,
            )
            assert math.isclose(
                math.atanh(estimate.correlation),
                float(expected_row.atanh_correlation),
                rel_tol=0.0,
                abs_tol=REML_REFERENCE_FISHER_Z_ABSOLUTE_TOLERANCE,
            )
            assert estimate.design_rank == int(expected_row.design_rank)
        else:
            assert estimate.status not in SUCCESS_STATUSES
            assert estimate.correlation is None
            assert str(expected_row.atanh_correlation_missing_kind) in {"NA", "NaN"}
        assert estimate.observed_value_count == int(expected_row.observed_value_count)


@pytest.mark.parametrize("fixture_id", FIXTURE_IDS)
def test_reml_consensus_matches_limma_reference_fixtures(fixture_id: str) -> None:
    matrix, design, blocks, feature_ids = _fixture_inputs(fixture_id)
    summary = _summary(fixture_id)

    result = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=feature_ids,
    )

    assert result.success
    assert result.consensus_correlation is not None
    assert math.isclose(
        result.consensus_correlation,
        float(summary["consensus_correlation"]),
        rel_tol=0.0,
        abs_tol=REML_REFERENCE_CORRELATION_ABSOLUTE_TOLERANCE,
    )
    assert math.isclose(
        math.atanh(result.consensus_correlation),
        float(summary["consensus_atanh_correlation"]),
        rel_tol=0.0,
        abs_tol=REML_REFERENCE_FISHER_Z_ABSOLUTE_TOLERANCE,
    )
    assert result.attempted_feature_count == int(summary["feature_count"])
    assert result.eligible_feature_count == int(summary["feature_count"])
    assert result.estimated_feature_count == int(
        summary["estimated_feature_correlation_count"]
    )
    assert result.failed_feature_count == int(
        summary["missing_feature_correlation_count"]
    )
    assert result.trim_fraction == DUPLICATE_CORRELATION_TRIM_FRACTION
    assert result.block_structure is not None
    assert result.block_structure.sample_count == int(summary["sample_count"])
    assert result.block_structure.block_count == int(summary["block_count"])
    assert result.block_structure.repeated_block_count == int(
        summary["repeated_block_count"]
    )
    assert result.block_structure.singleton_block_count == int(
        summary["singleton_block_count"]
    )
    assert result.block_structure.minimum_block_size == int(summary["min_block_size"])
    assert result.block_structure.maximum_block_size == int(summary["max_block_size"])
    assert result.design_rank == _design_rank_from_fixture(fixture_id)
    assert result.sample_count == int(summary["sample_count"])


def test_feature_level_failures_are_typed_and_counted_by_reason() -> None:
    fixture_id = "fixture_d_feature_level_failures"
    matrix, design, blocks, feature_ids = _fixture_inputs(fixture_id)

    result = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=feature_ids,
    )

    statuses = {
        estimate.feature_id: estimate.status for estimate in result.feature_estimates
    }
    assert statuses == {
        "D_valid_null": DuplicateCorrelationFeatureStatus.ESTIMATED,
        "D_valid_effect": DuplicateCorrelationFeatureStatus.ESTIMATED,
        "D_missing_still_estimable": DuplicateCorrelationFeatureStatus.ESTIMATED,
        "D_constant_all": (
            DuplicateCorrelationFeatureStatus.ZERO_OR_UNUSABLE_RESIDUAL_VARIATION
        ),
        "D_near_constant": DuplicateCorrelationFeatureStatus.ESTIMATED,
        "D_rank_loss_only_A": (
            DuplicateCorrelationFeatureStatus.LOST_FIXED_EFFECT_ESTIMABILITY
        ),
        "D_insufficient_one": (
            DuplicateCorrelationFeatureStatus.INSUFFICIENT_FINITE_OBSERVATIONS
        ),
        "D_all_missing": (
            DuplicateCorrelationFeatureStatus.INSUFFICIENT_FINITE_OBSERVATIONS
        ),
    }
    reason_counts = {item.reason: item.count for item in result.failure_reason_counts}
    assert reason_counts == {
        DuplicateCorrelationFailureReason.INSUFFICIENT_FINITE_OBSERVATIONS: 2,
        DuplicateCorrelationFailureReason.LOST_FIXED_EFFECT_ESTIMABILITY: 1,
        DuplicateCorrelationFailureReason.ZERO_OR_UNUSABLE_RESIDUAL_VARIATION: 1,
    }
    assert result.failed_feature_count == 4
    assert all(
        estimate.correlation is None
        for estimate in result.feature_estimates
        if estimate.status not in SUCCESS_STATUSES
    )


def test_boundary_convergence_is_recorded_but_still_contributes() -> None:
    matrix, design, blocks, feature_ids = _fixture_inputs("fixture_a_complete_pairs")
    result = estimate_duplicate_correlation_reml_consensus(
        matrix[:1, :],
        design,
        blocks,
        feature_ids=feature_ids[:1],
    )

    assert result.estimated_feature_count == 1
    assert {estimate.status for estimate in result.feature_estimates} == {
        DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED
    }
    assert {estimate.boundary for estimate in result.feature_estimates} == {
        DuplicateCorrelationBoundary.UPPER
    }
    assert result.convergence_summary is not None
    assert result.convergence_summary.boundary_feature_count == 1
    assert result.boundary_summary is not None
    assert result.boundary_summary.upper_boundary_feature_count == 1


def test_negative_correlation_can_approach_pair_block_pd_lower_bound() -> None:
    result = estimate_duplicate_correlation_reml_consensus(
        np.array([[1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]]),
        np.ones((8, 1), dtype=np.float64),
        ("b1", "b1", "b2", "b2", "b3", "b3", "b4", "b4"),
        feature_ids=("strong_negative",),
    )

    estimate = result.feature_estimates[0]
    assert result.success
    assert estimate.status is DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED
    assert estimate.boundary is DuplicateCorrelationBoundary.LOWER
    assert estimate.correlation is not None
    assert estimate.correlation == pytest.approx(-0.99)


def test_raw_correlation_above_upper_cap_is_clamped_deterministically(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw_fit = _variance_fit_for_raw_correlation(0.995)
    assert (
        raw_fit.correlation
        > duplicate_correlation_module.DUPLICATE_CORRELATION_UPPER_CORRELATION_CAP
    )

    def high_correlation_fit(**_: object) -> object:
        return raw_fit

    monkeypatch.setattr(
        duplicate_correlation_module,
        "_fit_reml_variance_components",
        high_correlation_fit,
    )
    matrix = np.array([[1.0, 2.0, 2.5, 1.5, 1.25, 2.25, 2.75, 1.75]])
    design = np.ones((8, 1), dtype=np.float64)
    blocks = ("b1", "b1", "b2", "b2", "b3", "b3", "b4", "b4")

    first = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=("upper_clamped",),
    )
    second = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=("upper_clamped",),
    )

    estimate = first.feature_estimates[0]
    assert first.success
    assert first.to_payload() == second.to_payload()
    assert estimate.status is DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED
    assert estimate.boundary is DuplicateCorrelationBoundary.UPPER
    assert estimate.correlation == pytest.approx(
        duplicate_correlation_module.DUPLICATE_CORRELATION_UPPER_CORRELATION_CAP
    )
    assert estimate.upper_correlation_bound == pytest.approx(
        duplicate_correlation_module.DUPLICATE_CORRELATION_UPPER_CORRELATION_CAP
    )
    assert first.convergence_summary is not None
    assert first.convergence_summary.boundary_feature_count == 1
    assert first.boundary_summary is not None
    assert first.boundary_summary.upper_boundary_feature_count == 1


def test_feature_missingness_changes_lower_clamp_and_boundary_diagnostics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw_fit = _variance_fit_for_raw_correlation(-0.75)

    def negative_correlation_fit(**_: object) -> object:
        return raw_fit

    monkeypatch.setattr(
        duplicate_correlation_module,
        "_fit_reml_variance_components",
        negative_correlation_fit,
    )
    result = estimate_duplicate_correlation_reml_consensus(
        np.array(
            [
                [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
                [1.0, 2.0, math.nan, 4.0, 5.0, math.nan, 7.0, 8.0, math.nan],
            ],
            dtype=np.float64,
        ),
        np.ones((9, 1), dtype=np.float64),
        ("b1", "b1", "b1", "b2", "b2", "b2", "b3", "b3", "b3"),
        feature_ids=("full_triples", "observed_pairs"),
    )

    estimates = {estimate.feature_id: estimate for estimate in result.feature_estimates}
    triple_lower = (-1.0 / 2.0) + (
        duplicate_correlation_module.DUPLICATE_CORRELATION_LOWER_CORRELATION_BOUND_OFFSET
    )
    pair_lower = -1.0 + (
        duplicate_correlation_module.DUPLICATE_CORRELATION_LOWER_CORRELATION_BOUND_OFFSET
    )

    full_triples = estimates["full_triples"]
    assert full_triples.status is DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED
    assert full_triples.boundary is DuplicateCorrelationBoundary.LOWER
    assert full_triples.lower_correlation_bound == pytest.approx(triple_lower)
    assert full_triples.correlation == pytest.approx(triple_lower)
    assert full_triples.correlation is not None
    assert full_triples.correlation > -1.0 / 2.0

    observed_pairs = estimates["observed_pairs"]
    assert observed_pairs.status is DuplicateCorrelationFeatureStatus.ESTIMATED
    assert observed_pairs.boundary is None
    assert observed_pairs.lower_correlation_bound == pytest.approx(pair_lower)
    assert observed_pairs.correlation == pytest.approx(raw_fit.correlation)
    assert observed_pairs.correlation is not None
    assert observed_pairs.correlation > -1.0

    assert result.boundary_summary is not None
    assert result.boundary_summary.lower_boundary_feature_count == 1
    assert result.boundary_summary.lower_correlation_bound == pytest.approx(
        (-1.0 / 2.0)
        + duplicate_correlation_module.DUPLICATE_CORRELATION_POSITIVE_DEFINITE_TOLERANCE
    )
    assert not result.success
    assert result.failure_reason is (
        DuplicateCorrelationFailureReason.INVALID_OR_NON_POSITIVE_DEFINITE_COVARIANCE
    )


def test_valid_negative_interior_raw_correlation_is_not_forced_to_zero(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    raw_fit = _variance_fit_for_raw_correlation(-0.25)

    def negative_correlation_fit(**_: object) -> object:
        return raw_fit

    monkeypatch.setattr(
        duplicate_correlation_module,
        "_fit_reml_variance_components",
        negative_correlation_fit,
    )

    result = estimate_duplicate_correlation_reml_consensus(
        np.array([[1.0, 2.0, 2.5, 1.5, 1.25, 2.25, 2.75, 1.75]]),
        np.ones((8, 1), dtype=np.float64),
        ("b1", "b1", "b2", "b2", "b3", "b3", "b4", "b4"),
        feature_ids=("negative_interior",),
    )

    estimate = result.feature_estimates[0]
    assert result.success
    assert estimate.status is DuplicateCorrelationFeatureStatus.ESTIMATED
    assert estimate.boundary is None
    assert estimate.correlation == pytest.approx(raw_fit.correlation)
    assert result.consensus_correlation == pytest.approx(raw_fit.correlation)


@pytest.mark.parametrize("maximum_block_size", (2, 3, 4, 8))
@pytest.mark.parametrize("which_clamp", ("lower", "upper"))
def test_gls_covariance_is_positive_definite_at_supported_clamps(
    maximum_block_size: int,
    which_clamp: str,
) -> None:
    lower = (
        -1.0 / float(maximum_block_size - 1)
        + duplicate_correlation_module.DUPLICATE_CORRELATION_LOWER_CORRELATION_BOUND_OFFSET
    )
    correlation = (
        lower
        if which_clamp == "lower"
        else duplicate_correlation_module.DUPLICATE_CORRELATION_UPPER_CORRELATION_CAP
    )
    one_block_covariance = np.full(
        (maximum_block_size, maximum_block_size),
        correlation,
        dtype=np.float64,
    )
    np.fill_diagonal(one_block_covariance, 1.0)
    assert float(np.min(np.linalg.eigvalsh(one_block_covariance))) > 0.0

    blocks = tuple(
        f"b{block_position}"
        for block_position in range(2)
        for _ in range(maximum_block_size)
    )
    sample_count = len(blocks)

    fit = duplicate_correlation_module.fit_compound_symmetry_gls(
        np.arange(1, sample_count + 1, dtype=np.float64)[np.newaxis, :],
        np.ones((sample_count, 1), dtype=np.float64),
        blocks,
        correlation,
        feature_ids=("supported_clamp",),
    )

    assert fit.consensus_correlation == pytest.approx(correlation)
    assert fit.feature_fit_statuses == ("fit",)


def test_exact_trim_count_semantics_for_small_and_large_feature_sets() -> None:
    small = tuple(math.tanh(value) for value in (-0.4, 0.0, 0.2, 0.3, 0.5, 0.7))
    small_consensus = _fisher_trimmed_mean_consensus(small)

    assert small_consensus.trim_count_each_tail == 0
    assert small_consensus.retained_count == 6
    assert math.isclose(
        small_consensus.correlation,
        math.tanh(math.fsum(math.atanh(value) for value in small) / 6.0),
    )

    large = tuple(math.tanh(float(index) / 10.0) for index in range(-10, 10))
    large_consensus = _fisher_trimmed_mean_consensus(tuple(reversed(large)))
    retained = sorted(math.atanh(value) for value in large)[3:17]

    assert large_consensus.trim_count_each_tail == 3
    assert large_consensus.retained_count == 14
    assert math.isclose(
        large_consensus.correlation,
        math.tanh(math.fsum(retained) / 14.0),
        rel_tol=0.0,
        abs_tol=STRICT_CORRELATION_ABSOLUTE_TOLERANCE,
    )


def test_zero_near_correlation_is_estimated_without_forcing_zero_fallback() -> None:
    matrix = np.array(
        [
            [
                -0.49687212,
                1.36149211,
                -1.35796556,
                -0.82130014,
                -0.68520838,
                0.18410005,
                0.18299914,
                -1.17608756,
            ]
        ]
    )
    design = np.ones((8, 1), dtype=np.float64)
    blocks = ("b1", "b1", "b2", "b2", "b3", "b3", "b4", "b4")

    result = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=("zero_near",),
    )

    estimate = result.feature_estimates[0]
    assert result.success
    assert estimate.status is DuplicateCorrelationFeatureStatus.ESTIMATED
    assert estimate.correlation is not None
    assert abs(estimate.correlation) <= 5.0e-4


def test_seeded_simulation_consensus_moves_with_known_correlation_strength() -> None:
    low = _simulate_duplicate_correlation_matrix(known_correlation=0.05)
    high = _simulate_duplicate_correlation_matrix(known_correlation=0.55)

    low_result = estimate_duplicate_correlation_reml_consensus(
        low[0],
        low[1],
        low[2],
        feature_ids=low[3],
    )
    high_result = estimate_duplicate_correlation_reml_consensus(
        high[0],
        high[1],
        high[2],
        feature_ids=high[3],
    )

    assert low_result.success
    assert high_result.success
    assert low_result.consensus_correlation is not None
    assert high_result.consensus_correlation is not None
    assert low_result.estimated_feature_count == 300
    assert high_result.estimated_feature_count == 300
    assert low_result.consensus_correlation < 0.15
    assert high_result.consensus_correlation > 0.45
    assert high_result.consensus_correlation > low_result.consensus_correlation + 0.30


def test_singleton_blocks_do_not_create_false_repeated_observations() -> None:
    result = estimate_duplicate_correlation_reml_consensus(
        np.array([[1.0, 2.0, 3.0, 4.0]]),
        np.ones((4, 1), dtype=np.float64),
        ("s1", "s2", "s3", "s4"),
        feature_ids=("all_singletons",),
    )

    assert not result.success
    assert result.failure_reason is DuplicateCorrelationFailureReason.NO_REPEATED_BLOCKS
    assert result.feature_estimates[0].status is (
        DuplicateCorrelationFeatureStatus.NO_REPEATED_OBSERVATIONS
    )
    assert result.block_structure is not None
    assert result.block_structure.repeated_block_count == 0
    assert result.block_structure.correlated_pair_count == 0


def test_renaming_blocks_without_changing_membership_leaves_estimates_unchanged() -> (
    None
):
    matrix, design, blocks, feature_ids = _fixture_inputs(
        "fixture_b_three_observation_blocks"
    )
    baseline = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=feature_ids,
    )
    renamed_by_block = {
        block: f"renamed_{position}"
        for position, block in enumerate(sorted(set(blocks)), start=1)
    }
    renamed_blocks = tuple(renamed_by_block[block] for block in blocks)

    renamed = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        renamed_blocks,
        feature_ids=feature_ids,
    )

    _assert_same_consensus_and_feature_correlations(baseline, renamed)


def test_permuting_samples_with_matrix_design_and_blocks_is_invariant() -> None:
    matrix, design, blocks, feature_ids = _fixture_inputs(
        "fixture_b_three_observation_blocks"
    )
    permutation = np.array(
        [
            5,
            1,
            9,
            0,
            7,
            13,
            3,
            11,
            15,
            2,
            17,
            4,
            19,
            6,
            21,
            8,
            23,
            10,
            12,
            14,
            16,
            18,
            20,
            22,
        ],
        dtype=np.int64,
    )
    baseline = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=feature_ids,
    )
    permuted = estimate_duplicate_correlation_reml_consensus(
        matrix[:, permutation],
        design[permutation, :],
        tuple(blocks[index] for index in permutation.tolist()),
        feature_ids=feature_ids,
    )

    _assert_same_consensus_and_feature_correlations(baseline, permuted)


def test_permuting_feature_order_leaves_consensus_unchanged() -> None:
    matrix, design, blocks, feature_ids = _fixture_inputs(
        "fixture_d_feature_level_failures"
    )
    baseline = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=feature_ids,
    )

    reordered = estimate_duplicate_correlation_reml_consensus(
        matrix[::-1, :],
        design,
        blocks,
        feature_ids=tuple(reversed(feature_ids)),
    )

    assert baseline.consensus_correlation is not None
    assert reordered.consensus_correlation is not None
    assert math.isclose(
        baseline.consensus_correlation,
        reordered.consensus_correlation,
        rel_tol=0.0,
        abs_tol=STRICT_CORRELATION_ABSOLUTE_TOLERANCE,
    )


def test_affine_rescaling_one_feature_does_not_change_its_correlation() -> None:
    matrix, design, blocks, feature_ids = _fixture_inputs(
        "fixture_b_three_observation_blocks"
    )
    baseline = estimate_duplicate_correlation_reml_consensus(
        matrix,
        design,
        blocks,
        feature_ids=feature_ids,
    )
    rescaled_matrix = matrix.copy()
    rescaled_matrix[0, :] = (rescaled_matrix[0, :] * 7.5) - 3.0

    rescaled = estimate_duplicate_correlation_reml_consensus(
        rescaled_matrix,
        design,
        blocks,
        feature_ids=feature_ids,
    )

    assert baseline.feature_estimates[0].correlation is not None
    assert rescaled.feature_estimates[0].correlation is not None
    assert math.isclose(
        baseline.feature_estimates[0].correlation,
        rescaled.feature_estimates[0].correlation,
        rel_tol=0.0,
        abs_tol=STRICT_CORRELATION_ABSOLUTE_TOLERANCE,
    )


def test_optimizer_non_convergence_is_a_typed_feature_outcome(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def non_converged_fit(**_: object) -> object:
        return duplicate_correlation_module._VarianceComponentFit(
            converged=False,
            residual_component=1.0,
            block_component=1.0,
            deviance=1.0,
            iteration_count=2,
        )

    monkeypatch.setattr(
        duplicate_correlation_module,
        "_fit_reml_variance_components",
        non_converged_fit,
    )
    matrix, design, blocks, feature_ids = _fixture_inputs("fixture_a_complete_pairs")

    result = estimate_duplicate_correlation_reml_consensus(
        matrix[:2, :],
        design,
        blocks,
        feature_ids=feature_ids[:2],
    )

    assert not result.success
    assert result.failure_reason is (
        DuplicateCorrelationFailureReason.NUMERICAL_NON_CONVERGENCE
    )
    assert {estimate.status for estimate in result.feature_estimates} == {
        DuplicateCorrelationFeatureStatus.OPTIMISATION_FAILED
    }
    assert {item.reason: item.count for item in result.failure_reason_counts} == {
        DuplicateCorrelationFailureReason.OPTIMISATION_FAILED: 2
    }


def test_non_finite_objective_is_distinct_from_optimizer_non_convergence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def non_finite_variance_fit(**_: object) -> object:
        return duplicate_correlation_module._VarianceComponentFit(
            converged=False,
            residual_component=math.nan,
            block_component=math.nan,
            deviance=math.inf,
            iteration_count=0,
        )

    monkeypatch.setattr(
        duplicate_correlation_module,
        "_fit_reml_variance_components",
        non_finite_variance_fit,
    )

    result = estimate_duplicate_correlation_reml_consensus(
        np.array([[1.0, 2.0, 2.0, 1.0, 1.5, 2.5, 2.5, 1.5]]),
        np.ones((8, 1), dtype=np.float64),
        ("b1", "b1", "b2", "b2", "b3", "b3", "b4", "b4"),
        feature_ids=("nonfinite_objective",),
    )

    estimate = result.feature_estimates[0]
    assert estimate.status is (
        DuplicateCorrelationFeatureStatus.NON_FINITE_OBJECTIVE_OR_ESTIMATE
    )
    assert estimate.failure_reason is (
        DuplicateCorrelationFailureReason.NON_FINITE_OBJECTIVE_OR_ESTIMATE
    )
    assert result.non_finite_feature_count == 1


def test_no_eligible_feature_correlations_fails_without_zero_fallback() -> None:
    design = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        dtype=np.float64,
    )

    result = estimate_duplicate_correlation_reml_consensus(
        np.array(
            [
                [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                [2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            ],
            dtype=np.float64,
        ),
        design,
        ("pair_1", "pair_1", "pair_2", "pair_2", "pair_3", "pair_3"),
        feature_ids=("constant_1", "constant_2"),
    )

    assert not result.success
    assert result.consensus_correlation is None
    assert result.failure_reason is (
        DuplicateCorrelationFailureReason.NO_FEATURE_WITH_ESTIMABLE_CORRELATION
    )
    assert result.estimated_feature_count == 0
    assert {estimate.status for estimate in result.feature_estimates} == {
        DuplicateCorrelationFeatureStatus.ZERO_OR_UNUSABLE_RESIDUAL_VARIATION
    }
    assert {item.reason: item.count for item in result.failure_reason_counts} == {
        DuplicateCorrelationFailureReason.ZERO_OR_UNUSABLE_RESIDUAL_VARIATION: 2
    }


def test_invalid_consensus_correlation_is_a_typed_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def invalid_consensus(correlations: tuple[float, ...]) -> object:
        return duplicate_correlation_module._TrimmedConsensus(
            correlation=0.995,
            trim_count_each_tail=0,
            retained_count=len(correlations),
        )

    monkeypatch.setattr(
        duplicate_correlation_module,
        "_fisher_trimmed_mean_consensus",
        invalid_consensus,
    )
    matrix, design, blocks, feature_ids = _fixture_inputs("fixture_a_complete_pairs")

    result = estimate_duplicate_correlation_reml_consensus(
        matrix[:2, :],
        design,
        blocks,
        feature_ids=feature_ids[:2],
    )

    assert not result.success
    assert result.consensus_correlation is None
    assert result.estimated_feature_count == 2
    assert result.failure_reason is (
        DuplicateCorrelationFailureReason.INVALID_OR_NON_POSITIVE_DEFINITE_COVARIANCE
    )


def test_invalid_inputs_fail_with_domain_specific_errors() -> None:
    with pytest.raises(PhosPyInputError, match="block_ids length"):
        estimate_duplicate_correlation_reml_consensus(
            np.array([[1.0, 2.0]]),
            np.ones((2, 1), dtype=np.float64),
            ("only_one_block_id",),
        )

    with pytest.raises(PhosPyInputError, match="rank deficient"):
        estimate_duplicate_correlation_reml_consensus(
            np.array([[1.0, 2.0, 3.0]]),
            np.ones((3, 2), dtype=np.float64),
            ("b1", "b1", "b2"),
        )

    with pytest.raises(PhosPyInputError, match="more than two residual degrees"):
        estimate_duplicate_correlation_reml_consensus(
            np.array([[1.0, 2.0, 3.0]]),
            np.array(
                [
                    [1.0, 0.0],
                    [0.0, 1.0],
                    [1.0, 0.0],
                ],
                dtype=np.float64,
            ),
            ("b1", "b1", "b2"),
        )


def test_duplicate_correlation_rejects_named_fixed_block_design_columns() -> None:
    with pytest.raises(PhosPyInputError, match="exclude fixed block dummy columns"):
        estimate_duplicate_correlation_reml_consensus(
            np.array([[1.0, 2.0, 3.0, 4.0]]),
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [1.0, 0.0, 1.0],
                    [0.0, 1.0, 1.0],
                ],
                dtype=np.float64,
            ),
            ("pair_1", "pair_1", "pair_2", "pair_2"),
            feature_ids=("site_1",),
            design_column_names=("A", "B", "block[pair_2]"),
        )


def test_status_and_count_models_reject_invalid_construction() -> None:
    with pytest.raises(PhosPyInputError, match="must identify the reached boundary"):
        DuplicateCorrelationFeatureEstimate(
            feature_id="site_1",
            status=DuplicateCorrelationFeatureStatus.BOUNDARY_CONVERGED,
            correlation=0.99,
            objective_value=1.0,
        )

    failed = DuplicateCorrelationFeatureEstimate(
        feature_id="site_1",
        status=DuplicateCorrelationFeatureStatus.INSUFFICIENT_FINITE_OBSERVATIONS,
        failure_reason=(
            DuplicateCorrelationFailureReason.INSUFFICIENT_FINITE_OBSERVATIONS
        ),
    )
    with pytest.raises(PhosPyInputError, match="failure_reason_counts"):
        DuplicateCorrelationConsensusResult(
            method=DUPLICATE_CORRELATION_METHOD_REML_FISHER_TRIMMED_MEAN,
            trim_fraction=DUPLICATE_CORRELATION_TRIM_FRACTION,
            success=False,
            consensus_correlation=None,
            eligible_feature_count=1,
            estimated_feature_count=0,
            failed_feature_count=1,
            non_finite_feature_count=0,
            failure_reason=(
                DuplicateCorrelationFailureReason.NO_FEATURE_WITH_ESTIMABLE_CORRELATION
            ),
            feature_estimates=(failed,),
            failure_reason_counts=(
                DuplicateCorrelationReasonCount(
                    reason=DuplicateCorrelationFailureReason.OPTIMISATION_FAILED,
                    count=1,
                ),
            ),
        )


def _fixture_inputs(
    fixture_id: str,
) -> tuple[
    npt.NDArray[np.float64], npt.NDArray[np.float64], tuple[str, ...], tuple[str, ...]
]:
    base = FIXTURE_ROOT / fixture_id
    matrix = pd.read_csv(base / "matrix.csv", index_col=0)
    design = pd.read_csv(base / "design.csv", index_col=0)
    blocks = pd.read_csv(base / "blocks.csv")["block_id"].astype(str)
    return (
        matrix.to_numpy(dtype=np.float64),
        design.to_numpy(dtype=np.float64),
        tuple(blocks.tolist()),
        tuple(str(value) for value in matrix.index.tolist()),
    )


def _summary(fixture_id: str) -> dict[str, str]:
    frame = pd.read_csv(FIXTURE_ROOT / fixture_id / "duplicate_correlation_summary.csv")
    return {str(row.field): str(row.value) for row in frame.itertuples(index=False)}


def _variance_fit_for_raw_correlation(correlation: float) -> object:
    return duplicate_correlation_module._VarianceComponentFit(
        converged=True,
        residual_component=1.0,
        block_component=correlation / (1.0 - correlation),
        deviance=1.0,
        iteration_count=1,
    )


def _simulate_duplicate_correlation_matrix(
    *,
    known_correlation: float,
) -> tuple[
    npt.NDArray[np.float64],
    npt.NDArray[np.float64],
    tuple[str, ...],
    tuple[str, ...],
]:
    rng = np.random.default_rng(20_260_819)
    feature_count = 300
    block_count = 10
    observations_per_block = 3
    sample_count = block_count * observations_per_block
    condition_positions = np.tile(
        np.arange(observations_per_block, dtype=np.int64),
        block_count,
    )
    block_positions = np.repeat(
        np.arange(block_count, dtype=np.int64),
        observations_per_block,
    )
    design = np.zeros((sample_count, observations_per_block), dtype=np.float64)
    design[np.arange(sample_count, dtype=np.int64), condition_positions] = 1.0
    feature_baseline = rng.normal(10.0, 1.0, size=(feature_count, 1))
    condition_effects = rng.normal(
        0.0,
        0.15,
        size=(feature_count, observations_per_block),
    )
    block_effects = rng.normal(
        0.0,
        math.sqrt(known_correlation),
        size=(feature_count, block_count),
    )
    residuals = rng.normal(
        0.0,
        math.sqrt(1.0 - known_correlation),
        size=(feature_count, sample_count),
    )
    matrix = np.asarray(
        feature_baseline
        + condition_effects[:, condition_positions]
        + block_effects[:, block_positions]
        + residuals,
        dtype=np.float64,
    )
    block_ids = tuple(f"block_{position}" for position in block_positions.tolist())
    feature_ids = tuple(f"simulated_{position}" for position in range(feature_count))
    return matrix, design, block_ids, feature_ids


def _design_rank_from_fixture(fixture_id: str) -> int:
    expected = pd.read_csv(
        FIXTURE_ROOT / fixture_id / "feature_correlations.csv",
        keep_default_na=False,
    )
    return int(expected["design_rank"].iloc[0])


def _assert_same_consensus_and_feature_correlations(
    first: DuplicateCorrelationConsensusResult,
    second: DuplicateCorrelationConsensusResult,
) -> None:
    assert first.consensus_correlation is not None
    assert second.consensus_correlation is not None
    assert math.isclose(
        first.consensus_correlation,
        second.consensus_correlation,
        rel_tol=0.0,
        abs_tol=STRICT_CORRELATION_ABSOLUTE_TOLERANCE,
    )
    first_correlations = {
        estimate.feature_id: estimate.correlation
        for estimate in first.feature_estimates
    }
    second_correlations = {
        estimate.feature_id: estimate.correlation
        for estimate in second.feature_estimates
    }
    assert first_correlations.keys() == second_correlations.keys()
    for feature_id, first_correlation in first_correlations.items():
        second_correlation = second_correlations[feature_id]
        assert first_correlation is not None
        assert second_correlation is not None
        assert math.isclose(
            first_correlation,
            second_correlation,
            rel_tol=0.0,
            abs_tol=STRICT_CORRELATION_ABSOLUTE_TOLERANCE,
        )
