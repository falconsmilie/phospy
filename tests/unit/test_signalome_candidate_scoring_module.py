from __future__ import annotations

import numpy as np
import pytest

import phospy.signalome.clustering.candidate_scoring as candidate_scoring_module
from phospy.signalome.clustering.candidate_scoring import (
    _CandidateClusterScoreResult,
    compute_candidate_cluster_scores,
    resolve_candidate_scoring_policy,
    summarize_profile_degeneracy,
)
from phospy.signalome.clustering.policies import (
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_TREE_ENGINE_EXACT,
)
from phospy.signalome.clustering.tree_building import (
    prepare_scoring_values_for_clustering,
)


def _values() -> np.ndarray:
    return np.asarray(
        [
            [1.0, 0.0, 0.5],
            [0.9, 0.1, 0.4],
            [0.0, 1.0, 0.5],
            [0.1, 0.9, 0.6],
        ],
        dtype=float,
    )


def test_resolve_candidate_scoring_policy_auto_switches_by_site_limit() -> None:
    assert (
        resolve_candidate_scoring_policy(
            scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
            candidate_scoring_policy=None,
            n_sites=10,
            max_full_candidate_scoring_sites=10,
        )
        == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    )
    assert (
        resolve_candidate_scoring_policy(
            scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
            candidate_scoring_policy=None,
            n_sites=11,
            max_full_candidate_scoring_sites=10,
        )
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )


def test_compute_candidate_cluster_scores_has_independent_full_and_sampled_paths() -> (
    None
):
    values = _values()
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)

    full_result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 3),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        max_exact_tree_sites=10,
        max_full_candidate_scoring_sites=10,
    )
    sampled_result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 3),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=10,
        max_full_candidate_scoring_sites=10,
    )

    assert isinstance(full_result, _CandidateClusterScoreResult)
    assert isinstance(sampled_result, _CandidateClusterScoreResult)
    assert full_result.candidate_scoring_mode == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    assert (
        sampled_result.candidate_scoring_mode
        == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    )
    assert full_result.candidate_scoring_evaluated is True
    assert sampled_result.candidate_scoring_evaluated is True
    assert full_result.candidate_scoring_sampling is None
    assert isinstance(sampled_result.candidate_scoring_sampling, dict)


@pytest.mark.parametrize(
    ("values", "expected_pair_count", "minimum_score"),
    [
        (np.asarray([[1.0, 2.0, 3.0]], dtype=float), 0, 0.0),
        (
            np.asarray(
                [
                    [1.0, 2.0, 3.0],
                    [1.0, 2.0, 3.1],
                ],
                dtype=float,
            ),
            1,
            0.99,
        ),
    ],
)
def test_candidate_scoring_handles_one_and_two_site_inputs(
    values: np.ndarray,
    expected_pair_count: int,
    minimum_score: float,
) -> None:
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)
    result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(1, 2),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )

    assert result.candidate_scoring_mode == SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    assert result.candidate_scoring_sampling is not None
    assert result.candidate_scoring_sampling["actual_sampled_pair_count"] == (
        expected_pair_count
    )
    assert result.candidate_scores[1].mean_median_correlation >= minimum_score


def test_candidate_scoring_excludes_degenerate_rows_and_scores_near_identical_profiles() -> (
    None
):
    values = np.asarray(
        [
            [1.0, 1.0, 1.0],  # constant
            [np.nan, np.nan, np.nan],  # all-missing
            [1.0, 1.0000001, 1.0000002],  # near-constant
            [0.0, 1.0, 2.0],  # informative
            [0.1, 1.1, 2.1],  # near-identical informative profile
        ],
        dtype=float,
    )
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)
    result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(1, 2),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )

    assert profile_degeneracy.zero_variance_count == 1
    assert profile_degeneracy.near_constant_count == 1
    assert profile_degeneracy.excluded_count == 3
    assert result.candidate_scores[1].min_median_correlation == pytest.approx(1.0)
    assert result.candidate_scores[1].mean_median_correlation == pytest.approx(1.0)


def test_candidate_scoring_sampled_and_full_paths_match_on_small_input() -> None:
    values = _values()
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)

    full_result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 4),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )
    sampled_result = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 4),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )

    assert set(full_result.candidate_scores) == set(sampled_result.candidate_scores)
    for cluster_count, full_score in full_result.candidate_scores.items():
        sampled_score = sampled_result.candidate_scores[cluster_count]
        assert sampled_score.min_median_correlation == pytest.approx(
            full_score.min_median_correlation
        )
        assert sampled_score.mean_median_correlation == pytest.approx(
            full_score.mean_median_correlation
        )


def test_candidate_scoring_sampled_path_is_deterministic_for_fixed_input() -> None:
    x = np.linspace(0.0, 1.0, 300, dtype=float)
    values = np.column_stack((x, x**2, np.sin(3.0 * x)))
    profile_degeneracy = summarize_profile_degeneracy(values)
    clustering_values = prepare_scoring_values_for_clustering(values)

    first = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 3),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=2000,
        max_full_candidate_scoring_sites=10,
    )
    second = compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=values,
        candidate_range=range(2, 3),
        profile_degeneracy=profile_degeneracy,
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=2000,
        max_full_candidate_scoring_sites=10,
    )

    assert first.candidate_scores[2] == second.candidate_scores[2]
    assert first.candidate_scoring_sampling == second.candidate_scoring_sampling


def test_sampled_candidate_scoring_matches_independent_hash_seed_draw_golden(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    values = np.asarray(
        [
            [1.0, 0.0, 2.0, 3.0],
            [2.0, 1.0, 0.0, 4.0],
            [3.0, 2.0, 1.0, 0.0],
            [4.0, 1.0, 3.0, 2.0],
            [5.0, 4.0, 2.0, 1.0],
        ],
        dtype=np.float64,
    )

    # This is a deliberately independent scalar formulation of the production
    # FNV-1a row hash and SplitMix64 mixer. Keeping it here pins the bit-level
    # sampling contract without using the helpers under test to derive expected
    # values.
    mask = (1 << 64) - 1

    def splitmix64(value: int) -> int:
        value = (value + 0x9E3779B97F4A7C15) & mask
        value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & mask
        value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & mask
        return (value ^ (value >> 31)) & mask

    row_hashes: list[int] = []
    for row in np.ascontiguousarray(values).view(np.uint64):
        row_hash = 1469598103934665603
        for word in row:
            row_hash = ((row_hash ^ int(word)) * 1099511628211) & mask
        row_hashes.append(splitmix64(row_hash))

    assert row_hashes == [
        2579086850102505202,
        10361442417904159495,
        12459094153805549769,
        18117290090067621064,
        14546578479251792247,
    ]
    ordered_hashes = sorted(row_hashes)
    seed_material = 0xD2B74407B1CE6E93 ^ len(row_hashes) ^ 3
    seed_material ^= int(
        np.bitwise_xor.reduce(np.asarray(ordered_hashes, dtype=np.uint64))
    )
    seed_material ^= int(sum(ordered_hashes) & mask)
    seed_material ^= ordered_hashes[len(ordered_hashes) // 2]
    seed = splitmix64(int(seed_material))
    assert seed == 16630141970914373421

    tie_breakers = [
        splitmix64(row_hash ^ 0xA0761D6478BD642F) for row_hash in row_hashes
    ]
    assert tie_breakers == [
        6819241095243295143,
        10756062657952495562,
        2542974076769639541,
        17174044461912123156,
        2812843493748259844,
    ]
    canonical_order = sorted(
        range(len(row_hashes)),
        key=lambda position: (row_hashes[position], tie_breakers[position]),
    )
    assert canonical_order == [0, 1, 2, 4, 3]
    sampled_offsets = np.random.default_rng(seed).choice(5, size=3, replace=False)
    assert sampled_offsets.tolist() == [3, 1, 0]
    sampled_positions = np.asarray(canonical_order, dtype=int)[sampled_offsets]
    assert sampled_positions.tolist() == [4, 1, 0]

    sampled_correlations = np.asarray(
        np.corrcoef(values[sampled_positions]),
        dtype=float,
    )
    sampled_pairs = sampled_correlations[np.triu_indices(3, k=1)]
    expected_score = float(np.median(sampled_pairs))
    assert sampled_pairs.tolist() == pytest.approx(
        [-0.3207134902949093, -0.8485281374238571, 0.529150262212918]
    )
    assert expected_score == pytest.approx(-0.3207134902949093)

    monkeypatch.setattr(
        candidate_scoring_module,
        "MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER",
        3,
    )
    result = compute_candidate_cluster_scores(
        clustering_values=prepare_scoring_values_for_clustering(values),
        correlation_values=values,
        candidate_range=range(1, 2),
        profile_degeneracy=summarize_profile_degeneracy(values),
        n_sites=values.shape[0],
        scoring_mode=SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        max_exact_tree_sites=1000,
        max_full_candidate_scoring_sites=1000,
    )

    score = result.candidate_scores[1]
    assert score.min_median_correlation == pytest.approx(expected_score)
    assert score.mean_median_correlation == pytest.approx(expected_score)
    assert result.candidate_scoring_sampling == {
        "sampling_cap": 3,
        "sampling_method": "deterministic_uniform_without_replacement",
        "deterministic_seed_policy": (
            "order_invariant_seed_from_row_hashes_and_sample_size"
        ),
        "actual_sampled_pair_count": 3,
        "per_cluster_sample_count_summary": {
            "min": 3,
            "max": 3,
            "mean": 3.0,
            "total": 3,
        },
    }
