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

# Central tolerance for independent Cholesky-profiled REML estimates versus
# pinned limma black-box fixtures. Boundary-adjacent positive correlations are
# capped at the fixture policy value of 0.99; the tolerance keeps non-boundary
# fixtures tight while allowing minor optimizer-path differences.
REML_REFERENCE_CORRELATION_ABSOLUTE_TOLERANCE = 8.0e-3
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
    expected = pd.read_csv(FIXTURE_ROOT / fixture_id / "feature_correlations.csv")

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
        else:
            assert estimate.status not in SUCCESS_STATUSES
            assert estimate.correlation is None


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


def test_boundary_convergence_is_recorded_but_still_contributes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def decreasing_objective(*, correlation: float, **_: object) -> float:
        return -correlation

    monkeypatch.setattr(
        duplicate_correlation_module,
        "_restricted_reml_objective",
        decreasing_objective,
    )

    result = estimate_duplicate_correlation_reml_consensus(
        np.array([[1.0, 2.0, 2.0, 1.0, 1.5, 2.5, 2.5, 1.5]]),
        np.ones((8, 1), dtype=np.float64),
        ("b1", "b1", "b2", "b2", "b3", "b3", "b4", "b4"),
        feature_ids=("boundary_feature",),
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
    assert -1.0 < estimate.correlation < -0.99


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


def test_optimizer_non_convergence_is_a_typed_feature_outcome() -> None:
    matrix, design, blocks, feature_ids = _fixture_inputs("fixture_a_complete_pairs")

    result = estimate_duplicate_correlation_reml_consensus(
        matrix[:2, :],
        design,
        blocks,
        feature_ids=feature_ids[:2],
        optimizer_max_iterations=1,
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
    def non_finite_objective(*, correlation: float, **_: object) -> float:
        return math.inf

    monkeypatch.setattr(
        duplicate_correlation_module,
        "_restricted_reml_objective",
        non_finite_objective,
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

    with pytest.raises(PhosPyInputError, match="at least two residual degrees"):
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


def _design_rank_from_fixture(fixture_id: str) -> int:
    expected = pd.read_csv(FIXTURE_ROOT / fixture_id / "feature_correlations.csv")
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
