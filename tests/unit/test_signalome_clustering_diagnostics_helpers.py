from __future__ import annotations

import pytest

from phospy.signalomes.clustering.diagnostics import (
    approximation_used_from_candidate_mode,
    build_candidate_scoring_sampling_provenance,
    build_module_selection_diagnostics,
)


def test_build_module_selection_diagnostics_normalizes_payload_types() -> None:
    diagnostics = build_module_selection_diagnostics(
        strategy="correlation_thresholds",
        selected_module_count=2.0,
        requested_module_count=3.0,
        threshold_used=0.5,
        max_clusters_evaluated=10.0,
        candidate_scores={2: object()},
        reason="ok",
        zero_variance_profile_count=1.0,
        near_constant_profile_count=2.0,
        excluded_from_correlation_count=3.0,
    )

    assert diagnostics.strategy == "correlation_thresholds"
    assert diagnostics.selected_module_count == 2
    assert diagnostics.requested_module_count == 3
    assert diagnostics.threshold_used == pytest.approx(0.5)
    assert diagnostics.max_clusters_evaluated == 10
    assert diagnostics.zero_variance_profile_count == 1
    assert diagnostics.near_constant_profile_count == 2
    assert diagnostics.excluded_from_correlation_count == 3
    assert 2 in diagnostics.candidate_scores


def test_build_candidate_scoring_sampling_provenance_has_stable_schema() -> None:
    provenance = build_candidate_scoring_sampling_provenance(
        max_sites_per_cluster=256,
        per_cluster_sample_counts=[32, 16, 8],
        actual_sampled_pair_count=1234,
        sampling_method="deterministic_uniform_without_replacement",
        deterministic_seed_policy="order_invariant_seed_from_row_hashes_and_sample_size",
    )

    assert provenance["sampling_cap"] == 256
    assert provenance["sampling_method"] == "deterministic_uniform_without_replacement"
    assert (
        provenance["deterministic_seed_policy"]
        == "order_invariant_seed_from_row_hashes_and_sample_size"
    )
    assert provenance["actual_sampled_pair_count"] == 1234
    summary = provenance["per_cluster_sample_count_summary"]
    assert summary == {"min": 8, "max": 32, "mean": pytest.approx(56 / 3), "total": 56}


def test_approximation_used_from_candidate_mode_requires_sampled_and_evaluated() -> (
    None
):
    assert approximation_used_from_candidate_mode(
        candidate_scoring_mode="sampled",
        candidate_scoring_evaluated=True,
    )
    assert not approximation_used_from_candidate_mode(
        candidate_scoring_mode="sampled",
        candidate_scoring_evaluated=False,
    )
    assert not approximation_used_from_candidate_mode(
        candidate_scoring_mode="full",
        candidate_scoring_evaluated=True,
    )
