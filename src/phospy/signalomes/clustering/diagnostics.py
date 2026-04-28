"""Helpers for clustering backend diagnostic metadata."""

from __future__ import annotations

import numpy as np

from phospy.signalomes.models import SignalomeModuleSelectionDiagnostics


def approximation_used_from_candidate_mode(
    *,
    candidate_scoring_mode: str,
    candidate_scoring_evaluated: bool,
) -> bool:
    """Return whether approximate candidate scoring was used."""

    return bool(candidate_scoring_evaluated and candidate_scoring_mode == "sampled")


def build_module_selection_diagnostics(
    *,
    strategy: str,
    selected_module_count: int,
    requested_module_count: int | None,
    threshold_used: float | None,
    max_clusters_evaluated: int,
    candidate_scores: dict[int, object],
    reason: str,
    zero_variance_profile_count: int,
    near_constant_profile_count: int,
    excluded_from_correlation_count: int,
) -> SignalomeModuleSelectionDiagnostics:
    """Build a normalized module-selection diagnostics payload."""

    return SignalomeModuleSelectionDiagnostics(
        strategy=str(strategy),
        selected_module_count=int(selected_module_count),
        requested_module_count=(
            None if requested_module_count is None else int(requested_module_count)
        ),
        threshold_used=(None if threshold_used is None else float(threshold_used)),
        max_clusters_evaluated=int(max_clusters_evaluated),
        candidate_scores=dict(candidate_scores),
        reason=str(reason),
        zero_variance_profile_count=int(zero_variance_profile_count),
        near_constant_profile_count=int(near_constant_profile_count),
        excluded_from_correlation_count=int(excluded_from_correlation_count),
    )


def build_candidate_scoring_sampling_provenance(
    *,
    max_sites_per_cluster: int,
    per_cluster_sample_counts: list[int],
    actual_sampled_pair_count: int,
    sampling_method: str,
    deterministic_seed_policy: str,
) -> dict[str, object]:
    """Build deterministic sampled candidate-scoring provenance metadata."""

    if per_cluster_sample_counts:
        sample_min = int(min(per_cluster_sample_counts))
        sample_max = int(max(per_cluster_sample_counts))
        sample_mean = float(np.mean(per_cluster_sample_counts))
        sample_total = int(sum(per_cluster_sample_counts))
    else:
        sample_min = 0
        sample_max = 0
        sample_mean = 0.0
        sample_total = 0

    return {
        "sampling_cap": int(max_sites_per_cluster),
        "sampling_method": str(sampling_method),
        "deterministic_seed_policy": str(deterministic_seed_policy),
        "actual_sampled_pair_count": int(actual_sampled_pair_count),
        "per_cluster_sample_count_summary": {
            "min": sample_min,
            "max": sample_max,
            "mean": sample_mean,
            "total": sample_total,
        },
    }


__all__ = [
    "approximation_used_from_candidate_mode",
    "build_candidate_scoring_sampling_provenance",
    "build_module_selection_diagnostics",
]
