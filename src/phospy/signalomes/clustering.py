from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, linkage
from sklearn.cluster import AgglomerativeClustering

from ..internal.types import SignalomeModuleSelectionStrategy
from ..validation.values.enums import validate_module_selection_strategy
from ..validation.values.numeric import validate_fraction, validate_positive_int

__all__ = [
    "ClusterCandidateScore",
    "ClusterSitesResult",
    "DEFAULT_SIGNALOME_MODULE_SELECTION_POLICY",
    "SignalomeModuleSelectionDiagnostics",
    "SignalomeModuleSelectionPolicy",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "select_module_count",
    "select_module_count_with_diagnostics",
]

MAX_FULL_CORRELATION_SITE_COUNT = 2000
MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER = 256
NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE = 1e-12


@dataclass(frozen=True, slots=True)
class ClusterCandidateScore:
    """Cached score summary for one candidate module count."""

    min_median_correlation: float
    mean_median_correlation: float


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionPolicy:
    """Explicit policy for automatic signalome module-count selection.

    ``strategy`` controls how the automatic selector behaves:

    - ``"correlation_thresholds"`` applies the current PhosPy correlation-based
      heuristic using the configured primary and fallback thresholds.
    - ``"single_module"`` bypasses automatic selection and forces one module
      unless the caller explicitly requests a module count.
    """

    strategy: SignalomeModuleSelectionStrategy = "correlation_thresholds"
    primary_threshold: float = 0.5
    fallback_threshold: float = 0.1
    max_clusters: int = 10

    def __post_init__(self) -> None:
        validate_module_selection_strategy(self.strategy)
        validate_fraction(self.primary_threshold, name="primary_threshold")
        validate_fraction(self.fallback_threshold, name="fallback_threshold")
        validate_positive_int(self.max_clusters, name="max_clusters")

    @classmethod
    def from_value(cls, value: object) -> SignalomeModuleSelectionPolicy:
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, Mapping):
            return cls(**dict(value))
        msg = "module_selection_policy must be a SignalomeModuleSelectionPolicy or mapping"
        raise TypeError(msg)


DEFAULT_SIGNALOME_MODULE_SELECTION_POLICY = SignalomeModuleSelectionPolicy()


@dataclass(frozen=True, slots=True)
class SignalomeModuleSelectionDiagnostics:
    """Structured explanation of how a signalome module count was chosen."""

    strategy: SignalomeModuleSelectionStrategy
    selected_module_count: int
    requested_module_count: int | None
    threshold_used: float | None
    max_clusters_evaluated: int
    candidate_scores: dict[int, ClusterCandidateScore]
    reason: str
    zero_variance_profile_count: int = 0
    near_constant_profile_count: int = 0
    excluded_from_correlation_count: int = 0

    @property
    def used_automatic_selection(self) -> bool:
        return self.requested_module_count is None


@dataclass(frozen=True, slots=True)
class ClusterSitesResult:
    """Cluster labels plus module-selection diagnostics."""

    site_clusters: pd.Series
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics


@dataclass(frozen=True, slots=True)
class _ModuleSelectionComputation:
    diagnostics: SignalomeModuleSelectionDiagnostics
    candidate_labels: dict[int, np.ndarray]


@dataclass(frozen=True, slots=True)
class _ProfileDegeneracySummary:
    zero_variance_count: int
    near_constant_count: int
    excluded_count: int
    excluded_mask: np.ndarray


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> pd.Series:
    """Cluster phosphosites into signalome site clusters."""

    return cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        policy=policy,
    ).site_clusters


def cluster_sites_with_diagnostics(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> ClusterSitesResult:
    """Cluster phosphosites and capture how the module count was chosen."""

    scoring_values = scoring_matrix.to_numpy(dtype=float)
    selection = _compute_module_selection(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        policy=policy,
    )
    diagnostics = selection.diagnostics
    n_sites = scoring_values.shape[0]
    module_count = max(1, min(diagnostics.selected_module_count, n_sites))

    if module_count == 1:
        labels = np.ones(n_sites, dtype=int)
    else:
        cached_labels = selection.candidate_labels.get(module_count)
        if cached_labels is not None:
            labels = cached_labels + 1
        else:
            labels = fit_cluster_labels(scoring_values, module_count) + 1

    return ClusterSitesResult(
        site_clusters=pd.Series(
            labels,
            index=scoring_matrix.index,
            dtype=int,
            name="site_cluster",
        ),
        module_selection_diagnostics=diagnostics,
    )


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> int:
    """Choose a signalome module count from the scoring matrix."""

    return select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        policy=policy,
    ).selected_module_count


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> SignalomeModuleSelectionDiagnostics:
    """Choose a module count and explain why that count was selected."""

    return _compute_module_selection(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        policy=policy,
    ).diagnostics


def _compute_module_selection(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    policy: SignalomeModuleSelectionPolicy | None = None,
) -> _ModuleSelectionComputation:
    """Evaluate module-count selection and keep reusable clustering artefacts."""

    resolved_policy = SignalomeModuleSelectionPolicy.from_value(policy)
    scoring_array = np.asarray(scoring_values, dtype=float)
    n_sites = scoring_array.shape[0]
    profile_degeneracy = summarize_profile_degeneracy(scoring_array)
    zero_variance_count = profile_degeneracy.zero_variance_count
    near_constant_count = profile_degeneracy.near_constant_count
    if n_sites <= 1:
        return _ModuleSelectionComputation(
            diagnostics=SignalomeModuleSelectionDiagnostics(
                strategy=resolved_policy.strategy,
                selected_module_count=1,
                requested_module_count=requested_module_count,
                threshold_used=None,
                max_clusters_evaluated=1,
                candidate_scores={},
                reason="single phosphosite input only supports one signalome module",
                zero_variance_profile_count=zero_variance_count,
                near_constant_profile_count=near_constant_count,
                excluded_from_correlation_count=0,
            ),
            candidate_labels={},
        )

    if requested_module_count is not None:
        resolved_count = max(1, min(int(requested_module_count), n_sites))
        return _ModuleSelectionComputation(
            diagnostics=SignalomeModuleSelectionDiagnostics(
                strategy=resolved_policy.strategy,
                selected_module_count=resolved_count,
                requested_module_count=int(requested_module_count),
                threshold_used=None,
                max_clusters_evaluated=min(resolved_policy.max_clusters, n_sites),
                candidate_scores={},
                reason="module_count was provided explicitly by the caller",
                zero_variance_profile_count=zero_variance_count,
                near_constant_profile_count=near_constant_count,
                excluded_from_correlation_count=0,
            ),
            candidate_labels={},
        )

    if resolved_policy.strategy == "single_module":
        return _ModuleSelectionComputation(
            diagnostics=SignalomeModuleSelectionDiagnostics(
                strategy=resolved_policy.strategy,
                selected_module_count=1,
                requested_module_count=None,
                threshold_used=None,
                max_clusters_evaluated=1,
                candidate_scores={},
                reason="module_selection_strategy='single_module' forces one module",
                zero_variance_profile_count=zero_variance_count,
                near_constant_profile_count=near_constant_count,
                excluded_from_correlation_count=0,
            ),
            candidate_labels={},
        )

    max_clusters = min(resolved_policy.max_clusters, n_sites)
    if max_clusters < 2:
        return _ModuleSelectionComputation(
            diagnostics=SignalomeModuleSelectionDiagnostics(
                strategy=resolved_policy.strategy,
                selected_module_count=1,
                requested_module_count=None,
                threshold_used=None,
                max_clusters_evaluated=max_clusters,
                candidate_scores={},
                reason="fewer than two cluster counts are available for evaluation",
                zero_variance_profile_count=zero_variance_count,
                near_constant_profile_count=near_constant_count,
                excluded_from_correlation_count=0,
            ),
            candidate_labels={},
        )

    candidate_range = range(2, max_clusters + 1)
    approximation_note = ""
    correlation_exclusion_note = build_correlation_exclusion_note(profile_degeneracy)
    if n_sites - profile_degeneracy.excluded_count <= 1:
        return _ModuleSelectionComputation(
            diagnostics=SignalomeModuleSelectionDiagnostics(
                strategy=resolved_policy.strategy,
                selected_module_count=1,
                requested_module_count=None,
                threshold_used=None,
                max_clusters_evaluated=1,
                candidate_scores={},
                reason=(
                    "fewer than two non-degenerate phosphosite profiles remained "
                    "after filtering degenerate rows for correlation scoring"
                )
                + correlation_exclusion_note,
                zero_variance_profile_count=zero_variance_count,
                near_constant_profile_count=near_constant_count,
                excluded_from_correlation_count=profile_degeneracy.excluded_count,
            ),
            candidate_labels={},
        )
    if n_sites <= MAX_FULL_CORRELATION_SITE_COUNT:
        site_correlations = build_correlation_matrix_with_exclusions(
            scoring_array,
            excluded_mask=profile_degeneracy.excluded_mask,
        )
        candidate_scores, candidate_labels = score_cluster_candidates(
            scoring_values=scoring_array,
            site_correlations=site_correlations,
            cluster_range=candidate_range,
        )
    else:
        candidate_scores, candidate_labels = score_cluster_candidates_approximate(
            scoring_values=scoring_array,
            cluster_range=candidate_range,
            max_sites_per_cluster=MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
        )
        approximation_note = (
            " Used sampled within-cluster correlation estimates to avoid "
            "materializing a full site-by-site correlation matrix."
        )
    primary_candidates = filter_cluster_candidates(
        candidate_scores,
        threshold=resolved_policy.primary_threshold,
    )
    if primary_candidates:
        selected_count = max(
            primary_candidates.items(),
            key=lambda item: (item[1], -item[0]),
        )[0]
        return _ModuleSelectionComputation(
            diagnostics=SignalomeModuleSelectionDiagnostics(
                strategy=resolved_policy.strategy,
                selected_module_count=selected_count,
                requested_module_count=None,
                threshold_used=resolved_policy.primary_threshold,
                max_clusters_evaluated=max_clusters,
                candidate_scores=candidate_scores,
                reason=(
                    "selected the highest-scoring candidate that satisfied the "
                    "primary within-cluster correlation threshold"
                )
                + correlation_exclusion_note
                + approximation_note,
                zero_variance_profile_count=zero_variance_count,
                near_constant_profile_count=near_constant_count,
                excluded_from_correlation_count=profile_degeneracy.excluded_count,
            ),
            candidate_labels=candidate_labels,
        )

    fallback_candidates = filter_cluster_candidates(
        candidate_scores,
        threshold=resolved_policy.fallback_threshold,
    )
    if fallback_candidates:
        selected_count = max(
            fallback_candidates.items(),
            key=lambda item: (item[1], -item[0]),
        )[0]
        return _ModuleSelectionComputation(
            diagnostics=SignalomeModuleSelectionDiagnostics(
                strategy=resolved_policy.strategy,
                selected_module_count=selected_count,
                requested_module_count=None,
                threshold_used=resolved_policy.fallback_threshold,
                max_clusters_evaluated=max_clusters,
                candidate_scores=candidate_scores,
                reason=(
                    "no candidate satisfied the primary threshold; selected the "
                    "highest-scoring fallback candidate"
                )
                + correlation_exclusion_note
                + approximation_note,
                zero_variance_profile_count=zero_variance_count,
                near_constant_profile_count=near_constant_count,
                excluded_from_correlation_count=profile_degeneracy.excluded_count,
            ),
            candidate_labels=candidate_labels,
        )

    return _ModuleSelectionComputation(
        diagnostics=SignalomeModuleSelectionDiagnostics(
            strategy=resolved_policy.strategy,
            selected_module_count=1,
            requested_module_count=None,
            threshold_used=None,
            max_clusters_evaluated=max_clusters,
            candidate_scores=candidate_scores,
            reason=(
                "no candidate module count satisfied the configured correlation "
                "thresholds, so the workflow fell back to one module"
            )
            + correlation_exclusion_note
            + approximation_note,
            zero_variance_profile_count=zero_variance_count,
            near_constant_profile_count=near_constant_count,
            excluded_from_correlation_count=profile_degeneracy.excluded_count,
        ),
        candidate_labels=candidate_labels,
    )


def summarize_profile_degeneracy(
    scoring_values: np.ndarray,
) -> _ProfileDegeneracySummary:
    """Classify profiles that cannot support robust Pearson correlations."""

    n_sites = int(np.asarray(scoring_values, dtype=float).shape[0])
    if n_sites == 0:
        return _ProfileDegeneracySummary(
            zero_variance_count=0,
            near_constant_count=0,
            excluded_count=0,
            excluded_mask=np.zeros(0, dtype=bool),
        )

    profile_variances = np.var(np.asarray(scoring_values, dtype=float), axis=1)
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
            msg = "excluded_mask must be a boolean vector aligned with scoring_values rows"
            raise ValueError(msg)

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


def filter_cluster_candidates(
    candidate_scores: dict[int, ClusterCandidateScore],
    *,
    threshold: float,
) -> dict[int, float]:
    """Return candidate counts whose cluster medians satisfy a threshold."""

    return {
        cluster_count: score.mean_median_correlation
        for cluster_count, score in candidate_scores.items()
        if score.min_median_correlation >= threshold
    }


def score_cluster_candidates(
    *,
    scoring_values: np.ndarray,
    site_correlations: np.ndarray,
    cluster_range: Iterable[int],
) -> tuple[dict[int, ClusterCandidateScore], dict[int, np.ndarray]]:
    """Score candidate counts and return reusable labels from one Ward hierarchy."""

    cluster_counts = [int(cluster_count) for cluster_count in cluster_range]
    if not cluster_counts:
        return {}, {}

    linkage_matrix = build_cluster_tree(scoring_values)
    candidate_labels = build_cluster_labels_from_tree(
        linkage_matrix=linkage_matrix,
        cluster_counts=cluster_counts,
    )

    candidates: dict[int, ClusterCandidateScore] = {}
    for cluster_count in cluster_counts:
        labels = candidate_labels[cluster_count]
        cluster_medians = [
            cluster_median_correlation(site_correlations, labels, label)
            for label in np.unique(labels)
        ]
        if not cluster_medians:
            continue
        candidates[cluster_count] = ClusterCandidateScore(
            min_median_correlation=float(min(cluster_medians)),
            mean_median_correlation=float(np.mean(cluster_medians)),
        )
    return candidates, candidate_labels


def score_cluster_candidates_approximate(
    *,
    scoring_values: np.ndarray,
    cluster_range: Iterable[int],
    max_sites_per_cluster: int,
) -> tuple[dict[int, ClusterCandidateScore], dict[int, np.ndarray]]:
    """Score candidate counts using sampled cluster-local correlations.

    This avoids materializing a full site-by-site correlation matrix for very
    large phosphosite sets.
    """

    cluster_counts = [int(cluster_count) for cluster_count in cluster_range]
    if not cluster_counts:
        return {}, {}

    linkage_matrix = build_cluster_tree(scoring_values)
    candidate_labels = build_cluster_labels_from_tree(
        linkage_matrix=linkage_matrix,
        cluster_counts=cluster_counts,
    )

    candidates: dict[int, ClusterCandidateScore] = {}
    for cluster_count in cluster_counts:
        labels = candidate_labels[cluster_count]
        cluster_medians = [
            cluster_median_correlation_approximate(
                scoring_values=scoring_values,
                labels=labels,
                label=label,
                max_sites_per_cluster=max_sites_per_cluster,
            )
            for label in np.unique(labels)
        ]
        if not cluster_medians:
            continue
        candidates[cluster_count] = ClusterCandidateScore(
            min_median_correlation=float(min(cluster_medians)),
            mean_median_correlation=float(np.mean(cluster_medians)),
        )
    return candidates, candidate_labels


def build_cluster_tree(scoring_values: np.ndarray) -> np.ndarray:
    """Build one Ward hierarchical tree for candidate module evaluation."""

    return linkage(np.asarray(scoring_values, dtype=float), method="ward")


def build_cluster_labels_from_tree(
    *,
    linkage_matrix: np.ndarray,
    cluster_counts: Iterable[int],
) -> dict[int, np.ndarray]:
    """Cut one cached hierarchy into labels for each candidate module count."""

    cluster_count_list = [int(cluster_count) for cluster_count in cluster_counts]
    if not cluster_count_list:
        return {}

    cut_labels = cut_tree(linkage_matrix, n_clusters=cluster_count_list)
    if cut_labels.ndim == 1:
        cut_labels = cut_labels.reshape(-1, 1)

    return {
        cluster_count: cut_labels[:, position].astype(int, copy=False)
        for position, cluster_count in enumerate(cluster_count_list)
    }


def fit_cluster_labels(scoring_values: np.ndarray, cluster_count: int) -> np.ndarray:
    """Fit Ward agglomerative clustering once for one candidate count."""

    return (
        AgglomerativeClustering(
            n_clusters=cluster_count,
            linkage="ward",
        )
        .fit_predict(scoring_values)
        .astype(int)
    )


def cluster_median_correlation(
    site_correlations: np.ndarray,
    labels: np.ndarray,
    label: int,
) -> float:
    """Return the median within-cluster correlation for one cluster label."""

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
    """Approximate the within-cluster median correlation for one label."""

    cluster_positions = np.flatnonzero(labels == label)
    if cluster_positions.size <= 1:
        return 0.0

    if cluster_positions.size > max_sites_per_cluster:
        sampled_positions = np.linspace(
            0,
            cluster_positions.size - 1,
            num=max_sites_per_cluster,
            dtype=int,
        )
        cluster_positions = cluster_positions[sampled_positions]

    cluster_values = scoring_values[cluster_positions]
    profile_degeneracy = summarize_profile_degeneracy(cluster_values)
    if cluster_values.shape[0] - profile_degeneracy.excluded_count <= 1:
        return 0.0

    cluster_correlations = build_correlation_matrix_with_exclusions(
        cluster_values,
        excluded_mask=profile_degeneracy.excluded_mask,
    ).copy()
    np.fill_diagonal(cluster_correlations, np.nan)
    values = cluster_correlations[~np.isnan(cluster_correlations)]
    if values.size == 0:
        return 0.0
    return float(np.median(values))
