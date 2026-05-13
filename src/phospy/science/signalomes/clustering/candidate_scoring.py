"""Candidate module scoring policies and helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from phospy.science.signalomes.clustering.diagnostic_schemas import (
    SignalomeCandidateScoringSamplingDiagnostics,
)
from phospy.science.signalomes.clustering.diagnostics import (
    build_candidate_scoring_sampling_provenance,
)
from phospy.science.signalomes.clustering.policies import (
    MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
    NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE,
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY,
    SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE,
    SIGNALOME_CLUSTERING_SCORING_MODE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
    _CandidateScoringMode,
)
from phospy.science.signalomes.clustering.scale_guards import (
    raise_if_full_candidate_scoring_limit_exceeded,
    resolve_max_exact_tree_sites,
)
from phospy.science.signalomes.clustering.tree_building import (
    ClusterTreeOperations,
    build_exact_cluster_tree_with_guard,
    resolve_cluster_tree_operations,
)
from phospy.science.signalomes.models import SignalomeClusterCandidateScore


@dataclass(frozen=True, slots=True)
class _ProfileDegeneracySummary:
    zero_variance_count: int
    near_constant_count: int
    excluded_count: int
    excluded_mask: np.ndarray


@dataclass(frozen=True, slots=True)
class _CandidateClusterScoreResult:
    candidate_scores: dict[int, SignalomeClusterCandidateScore]
    candidate_labels: dict[int, np.ndarray]
    approximation_note: str
    candidate_scoring_mode: _CandidateScoringMode
    exact_cluster_tree_built: bool
    candidate_scoring_evaluated: bool
    candidate_scoring_skip_reason: str | None
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None


def summarize_profile_degeneracy(
    scoring_values: np.ndarray,
) -> _ProfileDegeneracySummary:
    """Classify profiles that cannot support robust Pearson correlations."""

    values = np.asarray(scoring_values, dtype=float)
    n_sites = int(values.shape[0])
    if n_sites == 0:
        return _ProfileDegeneracySummary(
            zero_variance_count=0,
            near_constant_count=0,
            excluded_count=0,
            excluded_mask=np.zeros(0, dtype=bool),
        )

    profile_variances = np.var(values, axis=1)
    finite_mask = np.isfinite(profile_variances)
    zero_variance_mask = finite_mask & (profile_variances == 0.0)
    near_constant_mask = (
        finite_mask
        & (profile_variances > 0.0)
        & (profile_variances <= NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE)
    )
    excluded_mask = (~finite_mask) | zero_variance_mask | near_constant_mask
    return _ProfileDegeneracySummary(
        zero_variance_count=int(zero_variance_mask.sum(dtype=int)),
        near_constant_count=int(near_constant_mask.sum(dtype=int)),
        excluded_count=int(excluded_mask.sum(dtype=int)),
        excluded_mask=excluded_mask,
    )


def build_correlation_exclusion_note(summary: _ProfileDegeneracySummary) -> str:
    """Build a reason suffix describing degenerate profile exclusion."""

    if summary.excluded_count <= 0:
        return ""
    profile_label = "profile" if summary.excluded_count == 1 else "profiles"
    detail_tokens: list[str] = []
    if summary.zero_variance_count > 0:
        detail_tokens.append(f"{summary.zero_variance_count} zero-variance")
    if summary.near_constant_count > 0:
        detail_tokens.append(f"{summary.near_constant_count} near-constant")
    detail_suffix = f" ({', '.join(detail_tokens)})" if detail_tokens else ""
    return (
        f" Excluded {summary.excluded_count} degenerate {profile_label} from "
        f"correlation scoring{detail_suffix}."
    )


def build_correlation_matrix_with_exclusions(
    scoring_values: np.ndarray,
    *,
    excluded_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Compute a row-wise Pearson correlation matrix while excluding bad rows."""

    values = np.asarray(scoring_values, dtype=float)
    n_sites = int(values.shape[0])
    correlations = np.full((n_sites, n_sites), np.nan, dtype=float)
    if n_sites == 0:
        return correlations

    if excluded_mask is None:
        excluded = np.zeros(n_sites, dtype=bool)
    else:
        excluded = np.asarray(excluded_mask, dtype=bool)
        if excluded.shape != (n_sites,):
            raise ValueError(
                "excluded_mask must be a boolean vector aligned with scoring_values rows"
            )

    included_positions = np.flatnonzero(~excluded)
    if included_positions.size == 0:
        return correlations
    if included_positions.size == 1:
        correlations[included_positions[0], included_positions[0]] = 1.0
        return correlations

    included_correlations = np.corrcoef(values[included_positions])
    included_correlations = np.asarray(included_correlations, dtype=float)
    if included_correlations.ndim == 0:
        included_correlations = np.asarray(
            [[float(included_correlations)]],
            dtype=float,
        )
    np.fill_diagonal(included_correlations, 1.0)
    included_correlations = np.clip(included_correlations, -1.0, 1.0)
    correlations[np.ix_(included_positions, included_positions)] = included_correlations
    return correlations


def cluster_median_correlation(
    site_correlations: np.ndarray,
    labels: np.ndarray,
    label: int,
) -> float:
    """Return median within-cluster correlation for one cluster label."""

    cluster_positions = np.flatnonzero(labels == label)
    if cluster_positions.size <= 1:
        return 0.0
    cluster_correlations = site_correlations[
        np.ix_(cluster_positions, cluster_positions)
    ]
    cluster_correlations = cluster_correlations.copy()
    np.fill_diagonal(cluster_correlations, np.nan)
    values = cluster_correlations[~np.isnan(cluster_correlations)]
    if values.size == 0:
        return 0.0
    return float(np.median(values))


def cluster_median_correlation_approximate(
    *,
    scoring_values: np.ndarray,
    labels: np.ndarray,
    label: int,
    max_sites_per_cluster: int,
) -> float:
    """Approximate cluster-local median correlation using deterministic sampling."""

    correlation, _sampled_site_count, _sampled_pair_count = (
        _cluster_median_correlation_approximate_with_sampling_diagnostics(
            scoring_values=scoring_values,
            labels=labels,
            label=label,
            max_sites_per_cluster=max_sites_per_cluster,
        )
    )
    return float(correlation)


def _cluster_median_correlation_approximate_with_sampling_diagnostics(
    *,
    scoring_values: np.ndarray,
    labels: np.ndarray,
    label: int,
    max_sites_per_cluster: int,
) -> tuple[float, int, int]:
    """Return approximate median correlation plus sampled-size/pair diagnostics."""

    cluster_positions = np.flatnonzero(labels == label)
    sampled_site_count = int(cluster_positions.size)
    if cluster_positions.size <= 1:
        return 0.0, sampled_site_count, 0

    if cluster_positions.size > max_sites_per_cluster:
        cluster_positions = _sample_cluster_positions_for_approximation(
            scoring_values=scoring_values,
            cluster_positions=cluster_positions,
            sample_size=max_sites_per_cluster,
        )
    sampled_site_count = int(cluster_positions.size)

    cluster_values = np.asarray(scoring_values, dtype=float)[cluster_positions]
    profile_degeneracy = summarize_profile_degeneracy(cluster_values)
    if cluster_values.shape[0] - profile_degeneracy.excluded_count <= 1:
        return 0.0, sampled_site_count, 0

    cluster_correlations = build_correlation_matrix_with_exclusions(
        cluster_values,
        excluded_mask=profile_degeneracy.excluded_mask,
    ).copy()
    np.fill_diagonal(cluster_correlations, np.nan)
    values = cluster_correlations[~np.isnan(cluster_correlations)]
    if values.size == 0:
        return 0.0, sampled_site_count, 0
    return (
        float(np.median(values)),
        sampled_site_count,
        int(values.size // 2),
    )


def _build_candidate_scoring_sampling_provenance(
    *,
    max_sites_per_cluster: int,
    per_cluster_sample_counts: list[int],
    actual_sampled_pair_count: int,
) -> SignalomeCandidateScoringSamplingDiagnostics:
    """Build deterministic sampled candidate-scoring provenance metadata."""

    return build_candidate_scoring_sampling_provenance(
        max_sites_per_cluster=max_sites_per_cluster,
        per_cluster_sample_counts=per_cluster_sample_counts,
        actual_sampled_pair_count=actual_sampled_pair_count,
        sampling_method=SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
        deterministic_seed_policy=SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY,
    )


def _sample_cluster_positions_for_approximation(
    *,
    scoring_values: np.ndarray,
    cluster_positions: np.ndarray,
    sample_size: int,
) -> np.ndarray:
    """Sample cluster positions using order-invariant deterministic seeding."""

    cluster_values = np.asarray(scoring_values, dtype=np.float64)[cluster_positions]
    row_hashes = _stable_row_hashes(cluster_values)
    seed = int(_build_order_invariant_sampling_seed(row_hashes, sample_size))
    random_generator = np.random.default_rng(seed)

    tie_breakers = _splitmix64(row_hashes ^ np.uint64(0xA0761D6478BD642F))
    canonical_order = np.lexsort((tie_breakers, row_hashes))
    sampled_offsets = random_generator.choice(
        cluster_positions.size,
        size=sample_size,
        replace=False,
    )
    return cluster_positions[canonical_order[sampled_offsets]]


def _stable_row_hashes(values: np.ndarray) -> np.ndarray:
    """Build stable 64-bit hashes for rows in a float matrix."""

    matrix = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    if matrix.ndim != 2:
        raise ValueError("values must be a 2D matrix")
    if matrix.shape[0] == 0:
        return np.zeros(0, dtype=np.uint64)

    words = matrix.view(np.uint64)
    row_hashes = np.full(words.shape[0], np.uint64(1469598103934665603))
    fnv_prime = np.uint64(1099511628211)
    for word_index in range(words.shape[1]):
        row_hashes ^= words[:, word_index]
        row_hashes *= fnv_prime
    return _splitmix64(row_hashes)


def _build_order_invariant_sampling_seed(
    row_hashes: np.ndarray,
    sample_size: int,
) -> np.uint64:
    """Build a deterministic seed from order-invariant row summaries."""

    sorted_hashes = np.sort(np.asarray(row_hashes, dtype=np.uint64), kind="mergesort")
    seed = np.uint64(0xD2B74407B1CE6E93)
    seed ^= np.uint64(sorted_hashes.size)
    seed ^= np.uint64(sample_size)
    if sorted_hashes.size > 0:
        seed ^= np.bitwise_xor.reduce(sorted_hashes)
        seed ^= np.sum(sorted_hashes, dtype=np.uint64)
        seed ^= sorted_hashes[sorted_hashes.size // 2]
    return np.uint64(_splitmix64(np.asarray([seed], dtype=np.uint64))[0])


def _splitmix64(values: np.ndarray) -> np.ndarray:
    """Vectorized SplitMix64 mixer."""

    mixed = np.asarray(values, dtype=np.uint64).copy()
    mixed = mixed + np.uint64(0x9E3779B97F4A7C15)
    mixed = (mixed ^ (mixed >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    mixed = (mixed ^ (mixed >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return mixed ^ (mixed >> np.uint64(31))


def resolve_candidate_scoring_policy(
    *,
    scoring_mode: SignalomeClusteringScoringMode,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None,
    n_sites: int,
    max_full_candidate_scoring_sites: int,
) -> SignalomeCandidateScoringPolicy:
    if candidate_scoring_policy is not None:
        if (
            scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_EXACT
            and candidate_scoring_policy != SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
        ):
            raise ValueError(
                "scoring_mode='exact' cannot be combined with "
                "candidate_scoring_policy='sampled'"
            )
        if (
            scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE
            and candidate_scoring_policy != SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
        ):
            raise ValueError(
                "scoring_mode='approximate' cannot be combined with "
                "candidate_scoring_policy='full'"
            )
        return candidate_scoring_policy

    if scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_EXACT:
        return SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    if scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE:
        return SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED
    if n_sites <= int(max_full_candidate_scoring_sites):
        return SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    return SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED


def compute_candidate_cluster_scores(
    *,
    clustering_values: np.ndarray,
    correlation_values: np.ndarray,
    candidate_range: range,
    profile_degeneracy: _ProfileDegeneracySummary,
    n_sites: int,
    scoring_mode: SignalomeClusteringScoringMode,
    tree_engine: SignalomeTreeEngine,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy,
    max_exact_tree_sites: int | None,
    max_full_candidate_scoring_sites: int,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> _CandidateClusterScoreResult:
    """Score candidate cluster counts using full or sampled correlation paths."""

    candidate_counts = [int(cluster_count) for cluster_count in candidate_range]
    if not candidate_counts:
        return _CandidateClusterScoreResult(
            candidate_scores={},
            candidate_labels={},
            approximation_note="",
            candidate_scoring_mode=SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
            exact_cluster_tree_built=False,
            candidate_scoring_evaluated=False,
            candidate_scoring_skip_reason=None,
            candidate_scoring_sampling=None,
        )

    resolved_max_exact_tree_sites = resolve_max_exact_tree_sites(max_exact_tree_sites)
    # Guard ordering policy:
    # - If full candidate-correlation scoring exceeds max_full_candidate_scoring_sites
    #   while exact-tree construction is still permitted, fail here before any
    #   exact-tree construction is attempted.
    # - If both max_full_candidate_scoring_sites and max_exact_tree_sites are
    #   exceeded, defer to the exact-tree guard below as the canonical first
    #   failure for that configuration.
    raise_if_full_candidate_scoring_limit_exceeded(
        n_sites=n_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        max_exact_tree_sites=resolved_max_exact_tree_sites,
        candidate_scoring_policy=candidate_scoring_policy,
    )

    cluster_tree = build_exact_cluster_tree_with_guard(
        clustering_values=clustering_values,
        n_sites=n_sites,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=resolved_max_exact_tree_sites,
        cluster_tree_operations=cluster_tree_operations,
    )
    tree_operations = resolve_cluster_tree_operations(cluster_tree_operations)
    candidate_labels = tree_operations.build_cluster_labels_from_tree(
        cluster_tree=cluster_tree,
        cluster_counts=candidate_counts,
    )
    exact_cluster_tree_built = n_sites > 1

    if candidate_scoring_policy == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL:
        site_correlations = build_correlation_matrix_with_exclusions(
            correlation_values,
            excluded_mask=profile_degeneracy.excluded_mask,
        )
        candidate_scores: dict[int, SignalomeClusterCandidateScore] = {}
        for cluster_count in candidate_counts:
            labels = candidate_labels[cluster_count]
            cluster_medians = [
                cluster_median_correlation(site_correlations, labels, label)
                for label in np.unique(labels)
            ]
            if not cluster_medians:
                continue
            candidate_scores[cluster_count] = SignalomeClusterCandidateScore(
                min_median_correlation=float(min(cluster_medians)),
                mean_median_correlation=float(np.mean(cluster_medians)),
            )
        return _CandidateClusterScoreResult(
            candidate_scores=candidate_scores,
            candidate_labels=candidate_labels,
            approximation_note="",
            candidate_scoring_mode=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            exact_cluster_tree_built=exact_cluster_tree_built,
            candidate_scoring_evaluated=True,
            candidate_scoring_skip_reason=None,
            candidate_scoring_sampling=None,
        )

    candidate_scores = {}
    per_cluster_sample_counts: list[int] = []
    actual_sampled_pair_count = 0
    for cluster_count in candidate_counts:
        labels = candidate_labels[cluster_count]
        cluster_medians: list[float] = []
        for label in np.unique(labels):
            (
                cluster_median,
                sampled_site_count,
                sampled_pair_count,
            ) = _cluster_median_correlation_approximate_with_sampling_diagnostics(
                scoring_values=correlation_values,
                labels=labels,
                label=int(label),
                max_sites_per_cluster=MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
            )
            cluster_medians.append(cluster_median)
            per_cluster_sample_counts.append(int(sampled_site_count))
            actual_sampled_pair_count += int(sampled_pair_count)
        if not cluster_medians:
            continue
        candidate_scores[cluster_count] = SignalomeClusterCandidateScore(
            min_median_correlation=float(min(cluster_medians)),
            mean_median_correlation=float(np.mean(cluster_medians)),
        )
    if scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE:
        approximation_note = (
            " Used sampled within-cluster correlation estimates (seeded, "
            "order-invariant sampling) because candidate scoring was set to "
            "sampled. This sampling applies to candidate module-count evaluation "
            "only; exact cluster-tree construction and final module assignment "
            "remain exact."
        )
    else:
        approximation_note = (
            " Used sampled within-cluster correlation estimates (seeded, "
            "order-invariant sampling) to avoid materializing a full site-by-site "
            "correlation matrix during candidate module-count evaluation. Exact "
            "cluster-tree construction and final module assignment remain exact."
        )
    return _CandidateClusterScoreResult(
        candidate_scores=candidate_scores,
        candidate_labels=candidate_labels,
        approximation_note=approximation_note,
        candidate_scoring_mode=SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        exact_cluster_tree_built=exact_cluster_tree_built,
        candidate_scoring_evaluated=True,
        candidate_scoring_skip_reason=None,
        candidate_scoring_sampling=_build_candidate_scoring_sampling_provenance(
            max_sites_per_cluster=MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
            per_cluster_sample_counts=per_cluster_sample_counts,
            actual_sampled_pair_count=actual_sampled_pair_count,
        ),
    )


__all__ = [
    "build_correlation_exclusion_note",
    "build_correlation_matrix_with_exclusions",
    "cluster_median_correlation",
    "cluster_median_correlation_approximate",
    "compute_candidate_cluster_scores",
    "resolve_candidate_scoring_policy",
    "summarize_profile_degeneracy",
]
