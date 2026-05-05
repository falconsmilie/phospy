"""Helpers for clustering backend diagnostic metadata."""

from __future__ import annotations

import numpy as np
import pandas as pd  # pyright: ignore[reportMissingTypeStubs]

from phospy.signalomes.clustering.contracts import ClusterTreeEngine
from phospy.signalomes.clustering.diagnostic_schemas import (
    SignalomeBackendDiagnostics,
    SignalomeCandidateScoringSamplingDiagnostics,
    SignalomeClusteringLimitMetadata,
    SignalomeClusteringThresholdMetadata,
    SignalomeTreeEngineDiagnostics,
    build_backend_diagnostics,
    candidate_scoring_sampling_diagnostics_to_payload,
    validate_candidate_scoring_sampling_diagnostics,
)
from phospy.signalomes.models import (
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionDiagnostics,
    SignalomeModuleSelectionStrategy,
)


def approximation_used_from_candidate_mode(
    *,
    candidate_scoring_mode: str,
    candidate_scoring_evaluated: bool,
) -> bool:
    """Return whether approximate candidate scoring was used."""

    return bool(candidate_scoring_evaluated and candidate_scoring_mode == "sampled")


def build_module_selection_diagnostics(
    *,
    strategy: SignalomeModuleSelectionStrategy,
    selected_module_count: int,
    requested_module_count: int | None,
    threshold_used: float | None,
    max_clusters_evaluated: int,
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    reason: str,
    zero_variance_profile_count: int,
    near_constant_profile_count: int,
    excluded_from_correlation_count: int,
) -> SignalomeModuleSelectionDiagnostics:
    """Build a normalized module-selection diagnostics payload."""

    return SignalomeModuleSelectionDiagnostics(
        strategy=strategy,
        selected_module_count=int(selected_module_count),
        requested_module_count=(
            None if requested_module_count is None else int(requested_module_count)
        ),
        threshold_used=(None if threshold_used is None else float(threshold_used)),
        max_clusters_evaluated=int(max_clusters_evaluated),
        candidate_scores=candidate_scores.copy(),
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
) -> SignalomeCandidateScoringSamplingDiagnostics:
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

    payload = {
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
    return validate_candidate_scoring_sampling_diagnostics(
        payload,
        field_name="candidate_scoring_sampling",
    )


def candidate_scoring_sampling_provenance_to_dataframe(
    payload: SignalomeCandidateScoringSamplingDiagnostics,
) -> pd.DataFrame:
    """Convert sampled candidate-scoring diagnostics to a stable one-row DataFrame."""

    normalized = validate_candidate_scoring_sampling_diagnostics(
        payload,
        field_name="candidate_scoring_sampling",
    )
    summary = normalized["per_cluster_sample_count_summary"]
    return pd.DataFrame(
        [
            {
                "sampling_cap": int(normalized["sampling_cap"]),
                "sampling_method": str(normalized["sampling_method"]),
                "deterministic_seed_policy": str(
                    normalized["deterministic_seed_policy"]
                ),
                "actual_sampled_pair_count": int(
                    normalized["actual_sampled_pair_count"]
                ),
                "per_cluster_sample_count_min": int(summary["min"]),
                "per_cluster_sample_count_max": int(summary["max"]),
                "per_cluster_sample_count_mean": float(summary["mean"]),
                "per_cluster_sample_count_total": int(summary["total"]),
            }
        ]
    )


def candidate_scoring_sampling_provenance_to_payload(
    payload: SignalomeCandidateScoringSamplingDiagnostics,
) -> dict[str, object]:
    """Serialize sampled candidate-scoring diagnostics to a stable payload."""

    return candidate_scoring_sampling_diagnostics_to_payload(payload)


class SignalomeEngineDiagnosticsBuilder:
    """Build backend diagnostics and typed metadata for engine results."""

    def backend_diagnostics(
        self,
        *,
        clustering_engine: str,
        backend_version: str,
        tree_engine: ClusterTreeEngine,
        tree_engine_diagnostics: SignalomeTreeEngineDiagnostics,
        selected_module_count: int,
        input_site_count: int,
        exact_tree_path_used: bool,
    ) -> SignalomeBackendDiagnostics:
        return build_backend_diagnostics(
            backend_name=str(clustering_engine),
            backend_version=str(backend_version),
            tree_engine=str(tree_engine.name),
            tree_engine_version=str(tree_engine.version),
            tree_engine_diagnostics=tree_engine_diagnostics,
            selected_module_count=selected_module_count,
            input_site_count=input_site_count,
            exact_tree_path_used=exact_tree_path_used,
        )

    def threshold_metadata(
        self, *, primary_threshold: float, fallback_threshold: float
    ) -> SignalomeClusteringThresholdMetadata:
        return {
            "primary_threshold": float(primary_threshold),
            "fallback_threshold": float(fallback_threshold),
        }

    def limit_metadata(
        self,
        *,
        max_exact_tree_sites: int | None,
        max_full_candidate_scoring_sites: int,
        max_clusters: int,
    ) -> SignalomeClusteringLimitMetadata:
        return {
            "max_exact_tree_sites": (
                None if max_exact_tree_sites is None else int(max_exact_tree_sites)
            ),
            "max_full_candidate_scoring_sites": int(max_full_candidate_scoring_sites),
            "max_clusters": int(max_clusters),
        }


__all__ = [
    "SignalomeEngineDiagnosticsBuilder",
    "approximation_used_from_candidate_mode",
    "build_candidate_scoring_sampling_provenance",
    "build_module_selection_diagnostics",
    "candidate_scoring_sampling_provenance_to_dataframe",
    "candidate_scoring_sampling_provenance_to_payload",
]
