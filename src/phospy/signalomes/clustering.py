"""Signalome clustering and module-count selection helpers."""

from __future__ import annotations

import heapq
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from phospy.errors.workflows import SignalomeScaleError
from phospy.signalomes.constants import (
    MODULE_ID_COLUMN,
    PROTEIN_COLUMN,
    SITE_CLUSTER_COLUMN,
    SITE_ID_COLUMN,
)
from phospy.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionDiagnostics,
    SignalomeModuleSelectionStrategy,
)

# Performance contracts for module-count selection scoring:
# - At or below `MAX_FULL_CORRELATION_SITE_COUNT`, candidate scoring computes a
#   full site-by-site correlation matrix.
# - Above this threshold, candidate scoring uses sampled within-cluster
#   correlations with at most `MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER` sites
#   per cluster.
#
# These thresholds only control how module-selection scores are computed. They do
# not change the input scoring matrix, selected output table schema, or whether
# approximation use is surfaced in diagnostics (`diagnostics.reason`).
MAX_FULL_CORRELATION_SITE_COUNT = 2000
MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER = 256
NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE = 1e-12
SIGNALOME_CLUSTERING_SCORING_MODE_AUTO = "auto"
SIGNALOME_CLUSTERING_SCORING_MODE_EXACT = "exact"
SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE = "approximate"
SignalomeClusteringScoringMode = Literal["auto", "exact", "approximate"]
SIGNALOME_CLUSTER_TREE_BACKEND_EXACT = "exact"
SignalomeClusterTreeBackend = Literal["exact"]
SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL = "full"
SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED = "sampled"
SignalomeCandidateScoringBackend = Literal["full", "sampled"]
SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED = "not_evaluated"


@dataclass(frozen=True, slots=True)
class ClusterSitesResult:
    """Cluster labels and diagnostics for module-count selection."""

    site_clusters: pd.Series
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics
    cluster_tree_backend: str = SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
    candidate_scoring_mode: str = SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    exact_cluster_tree_built: bool = False


@dataclass(frozen=True, slots=True)
class _ModuleSelectionComputation:
    diagnostics: SignalomeModuleSelectionDiagnostics
    candidate_labels: dict[int, np.ndarray]
    cluster_tree_backend: str
    candidate_scoring_mode: str
    exact_cluster_tree_built: bool


@dataclass(frozen=True, slots=True)
class _ProfileDegeneracySummary:
    zero_variance_count: int
    near_constant_count: int
    excluded_count: int
    excluded_mask: np.ndarray


@dataclass(frozen=True, slots=True)
class _WardClusterTree:
    n_sites: int
    merges: tuple[tuple[int, int], ...]


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = (
        SIGNALOME_CLUSTERING_SCORING_MODE_AUTO
    ),
    cluster_tree_backend: SignalomeClusterTreeBackend = (
        SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
    ),
    candidate_scoring_backend: SignalomeCandidateScoringBackend | None = None,
    max_exact_cluster_tree_sites: int | None = None,
    max_full_correlation_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
) -> pd.Series:
    """Cluster phosphosites into site clusters."""

    return cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
        max_full_correlation_sites=max_full_correlation_sites,
    ).site_clusters


def cluster_sites_with_diagnostics(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = (
        SIGNALOME_CLUSTERING_SCORING_MODE_AUTO
    ),
    cluster_tree_backend: SignalomeClusterTreeBackend = (
        SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
    ),
    candidate_scoring_backend: SignalomeCandidateScoringBackend | None = None,
    max_exact_cluster_tree_sites: int | None = None,
    max_full_correlation_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
) -> ClusterSitesResult:
    """Cluster phosphosites and capture module-selection diagnostics."""

    scoring_values = np.asarray(scoring_matrix.to_numpy(dtype=float, copy=False))
    selection = _compute_module_selection(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
        max_full_correlation_sites=max_full_correlation_sites,
    )
    diagnostics = selection.diagnostics
    n_sites = int(scoring_values.shape[0])
    module_count = max(1, min(int(diagnostics.selected_module_count), n_sites))

    if module_count == 1:
        labels = np.ones(n_sites, dtype=int)
    else:
        cached_labels = selection.candidate_labels.get(module_count)
        if cached_labels is not None:
            labels = cached_labels.astype(int, copy=False) + 1
        else:
            labels = (
                fit_cluster_labels(
                    scoring_values=_prepare_scoring_values_for_clustering(
                        scoring_values
                    ),
                    cluster_count=module_count,
                    cluster_tree_backend=cluster_tree_backend,
                    candidate_scoring_backend=(
                        SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL
                        if candidate_scoring_backend is None
                        else candidate_scoring_backend
                    ),
                    max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
                )
                + 1
            )

    return ClusterSitesResult(
        site_clusters=pd.Series(
            labels,
            index=scoring_matrix.index.copy(),
            dtype=int,
            name=SITE_CLUSTER_COLUMN,
        ),
        module_selection_diagnostics=diagnostics,
        cluster_tree_backend=selection.cluster_tree_backend,
        candidate_scoring_mode=selection.candidate_scoring_mode,
        exact_cluster_tree_built=selection.exact_cluster_tree_built,
    )


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = (
        SIGNALOME_CLUSTERING_SCORING_MODE_AUTO
    ),
) -> int:
    """Select a module count from a scoring matrix."""

    return select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
    ).selected_module_count


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = (
        SIGNALOME_CLUSTERING_SCORING_MODE_AUTO
    ),
) -> SignalomeModuleSelectionDiagnostics:
    """Select a module count and return diagnostics."""

    array = (
        scoring_values.to_numpy(dtype=float, copy=False)
        if isinstance(scoring_values, pd.DataFrame)
        else np.asarray(scoring_values, dtype=float)
    )
    return _compute_module_selection(
        scoring_values=array,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
    ).diagnostics


def derive_protein_modules(
    *,
    site_clusters: pd.Series,
    site_to_protein: pd.Series,
) -> pd.Series:
    """Collapse site-level clusters into protein-level module assignments."""

    aligned_site_to_protein = site_to_protein.copy()
    aligned_site_to_protein.index = pd.Index(
        aligned_site_to_protein.index.astype(str),
        name=SITE_ID_COLUMN,
    )
    cluster_index = pd.Index(site_clusters.index.astype(str), name=SITE_ID_COLUMN)
    missing_sites = [
        site_id
        for site_id in cluster_index
        if site_id not in aligned_site_to_protein.index
    ]
    if missing_sites:
        preview = ", ".join(missing_sites[:3])
        suffix = "..." if len(missing_sites) > 3 else ""
        raise ValueError(
            f"site_to_protein is missing clustered site mappings: {preview}{suffix}"
        )
    aligned_site_to_protein = aligned_site_to_protein.loc[cluster_index].astype(str)

    proteins = pd.Index(aligned_site_to_protein.tolist(), dtype=object)
    membership = pd.crosstab(site_clusters, proteins)
    membership = (membership > 0).astype(int)

    assignments: dict[str, int] = {}
    pattern_to_module: dict[tuple[int, ...], int] = {}
    next_module_id = 1
    for protein in membership.columns:
        pattern = tuple(int(value) for value in membership.loc[:, protein].tolist())
        if pattern not in pattern_to_module:
            pattern_to_module[pattern] = next_module_id
            next_module_id += 1
        assignments[str(protein)] = pattern_to_module[pattern]

    protein_modules = pd.Series(assignments, dtype="int64", name=MODULE_ID_COLUMN)
    protein_modules.index.name = PROTEIN_COLUMN
    return protein_modules


def _compute_module_selection(
    *,
    scoring_values: np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = (
        SIGNALOME_CLUSTERING_SCORING_MODE_AUTO
    ),
    cluster_tree_backend: SignalomeClusterTreeBackend = (
        SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
    ),
    candidate_scoring_backend: SignalomeCandidateScoringBackend | None = None,
    max_exact_cluster_tree_sites: int | None = None,
    max_full_correlation_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
) -> _ModuleSelectionComputation:
    _validate_threshold(primary_threshold, field_name="primary_threshold")
    _validate_threshold(fallback_threshold, field_name="fallback_threshold")
    if max_clusters < 1:
        raise ValueError("max_clusters must be >= 1")
    if scoring_mode not in {
        SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
        SIGNALOME_CLUSTERING_SCORING_MODE_EXACT,
        SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE,
    }:
        raise ValueError("scoring_mode must be one of: auto, exact, approximate")
    if cluster_tree_backend != SIGNALOME_CLUSTER_TREE_BACKEND_EXACT:
        raise ValueError("cluster_tree_backend must be 'exact'")
    if candidate_scoring_backend not in {
        None,
        SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
        SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
    }:
        raise ValueError("candidate_scoring_backend must be one of: full, sampled")
    if (
        max_exact_cluster_tree_sites is not None
        and int(max_exact_cluster_tree_sites) < 1
    ):
        raise ValueError("max_exact_cluster_tree_sites must be >= 1")
    if max_full_correlation_sites < 1:
        raise ValueError("max_full_correlation_sites must be >= 1")

    scoring_array = np.asarray(scoring_values, dtype=float)
    if scoring_array.ndim != 2:
        raise ValueError("scoring_values must be a 2D array")
    n_sites = int(scoring_array.shape[0])
    profile_degeneracy = summarize_profile_degeneracy(scoring_array)
    correlation_exclusion_note = build_correlation_exclusion_note(profile_degeneracy)

    early_selection, resolved_max_clusters = _resolve_pre_scoring_module_selection(
        requested_module_count=requested_module_count,
        n_sites=n_sites,
        max_clusters=max_clusters,
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
    )
    if early_selection is not None:
        return early_selection
    resolved_candidate_scoring_backend = _resolve_candidate_scoring_backend(
        scoring_mode=scoring_mode,
        candidate_scoring_backend=candidate_scoring_backend,
        n_sites=n_sites,
        max_full_correlation_sites=max_full_correlation_sites,
    )

    clustering_values = _prepare_scoring_values_for_clustering(scoring_array)
    (
        candidate_scores,
        candidate_labels,
        approximation_note,
        candidate_scoring_mode,
        exact_cluster_tree_built,
    ) = _compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=scoring_array,
        candidate_range=range(2, resolved_max_clusters + 1),
        profile_degeneracy=profile_degeneracy,
        n_sites=n_sites,
        scoring_mode=scoring_mode,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=resolved_candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
        max_full_correlation_sites=max_full_correlation_sites,
    )

    primary_selection = _select_threshold_candidate(
        candidate_scores=candidate_scores,
        candidate_labels=candidate_labels,
        max_clusters=resolved_max_clusters,
        threshold=primary_threshold,
        requested_module_count=requested_module_count,
        reason=(
            "selected the highest-scoring candidate that satisfied the primary "
            "within-cluster correlation threshold"
        ),
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
        approximation_note=approximation_note,
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        cluster_tree_backend=cluster_tree_backend,
    )
    if primary_selection is not None:
        return primary_selection

    fallback_selection = _select_threshold_candidate(
        candidate_scores=candidate_scores,
        candidate_labels=candidate_labels,
        max_clusters=resolved_max_clusters,
        threshold=fallback_threshold,
        requested_module_count=requested_module_count,
        reason=(
            "no candidate satisfied the primary threshold; selected the "
            "highest-scoring fallback candidate"
        ),
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
        approximation_note=approximation_note,
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        cluster_tree_backend=cluster_tree_backend,
    )
    if fallback_selection is not None:
        return fallback_selection

    return _build_module_selection_result(
        strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        selected_module_count=1,
        requested_module_count=requested_module_count,
        threshold_used=None,
        max_clusters_evaluated=resolved_max_clusters,
        candidate_scores=candidate_scores,
        reason=(
            "no candidate module count satisfied the configured correlation "
            "thresholds, so the workflow fell back to one module"
        )
        + correlation_exclusion_note
        + approximation_note,
        profile_degeneracy=profile_degeneracy,
        excluded_from_correlation_count=profile_degeneracy.excluded_count,
        candidate_labels=candidate_labels,
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        cluster_tree_backend=cluster_tree_backend,
    )


def _build_module_selection_result(
    *,
    strategy: SignalomeModuleSelectionStrategy,
    selected_module_count: int,
    requested_module_count: int | None,
    threshold_used: float | None,
    max_clusters_evaluated: int,
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    reason: str,
    profile_degeneracy: _ProfileDegeneracySummary,
    excluded_from_correlation_count: int,
    candidate_labels: dict[int, np.ndarray],
    cluster_tree_backend: str = SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
    candidate_scoring_mode: str = SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    exact_cluster_tree_built: bool = False,
) -> _ModuleSelectionComputation:
    return _ModuleSelectionComputation(
        diagnostics=SignalomeModuleSelectionDiagnostics(
            strategy=strategy,
            selected_module_count=int(selected_module_count),
            requested_module_count=(
                None if requested_module_count is None else int(requested_module_count)
            ),
            threshold_used=threshold_used,
            max_clusters_evaluated=int(max_clusters_evaluated),
            candidate_scores=dict(candidate_scores),
            reason=str(reason),
            zero_variance_profile_count=int(profile_degeneracy.zero_variance_count),
            near_constant_profile_count=int(profile_degeneracy.near_constant_count),
            excluded_from_correlation_count=int(excluded_from_correlation_count),
        ),
        candidate_labels=candidate_labels,
        cluster_tree_backend=str(cluster_tree_backend),
        candidate_scoring_mode=str(candidate_scoring_mode),
        exact_cluster_tree_built=bool(exact_cluster_tree_built),
    )


def _resolve_pre_scoring_module_selection(
    *,
    requested_module_count: int | None,
    n_sites: int,
    max_clusters: int,
    profile_degeneracy: _ProfileDegeneracySummary,
    correlation_exclusion_note: str,
) -> tuple[_ModuleSelectionComputation | None, int]:
    if n_sites <= 1:
        return (
            _build_module_selection_result(
                strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
                selected_module_count=1,
                requested_module_count=requested_module_count,
                threshold_used=None,
                max_clusters_evaluated=1,
                candidate_scores={},
                reason="single phosphosite input only supports one signalome module",
                profile_degeneracy=profile_degeneracy,
                excluded_from_correlation_count=0,
                candidate_labels={},
            ),
            1,
        )

    if requested_module_count is not None:
        resolved_count = max(1, min(int(requested_module_count), n_sites))
        return (
            _build_module_selection_result(
                strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_EXPLICIT_MODULE_COUNT,
                selected_module_count=resolved_count,
                requested_module_count=int(requested_module_count),
                threshold_used=None,
                max_clusters_evaluated=min(int(max_clusters), n_sites),
                candidate_scores={},
                reason="module_count was provided explicitly by the caller",
                profile_degeneracy=profile_degeneracy,
                excluded_from_correlation_count=0,
                candidate_labels={},
            ),
            1,
        )

    resolved_max_clusters = min(int(max_clusters), n_sites)
    if resolved_max_clusters < 2:
        return (
            _build_module_selection_result(
                strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
                selected_module_count=1,
                requested_module_count=None,
                threshold_used=None,
                max_clusters_evaluated=resolved_max_clusters,
                candidate_scores={},
                reason="fewer than two cluster counts are available for evaluation",
                profile_degeneracy=profile_degeneracy,
                excluded_from_correlation_count=0,
                candidate_labels={},
            ),
            resolved_max_clusters,
        )

    if n_sites - profile_degeneracy.excluded_count <= 1:
        return (
            _build_module_selection_result(
                strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
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
                profile_degeneracy=profile_degeneracy,
                excluded_from_correlation_count=profile_degeneracy.excluded_count,
                candidate_labels={},
            ),
            resolved_max_clusters,
        )

    return None, resolved_max_clusters


def _compute_candidate_cluster_scores(
    *,
    clustering_values: np.ndarray,
    correlation_values: np.ndarray,
    candidate_range: range,
    profile_degeneracy: _ProfileDegeneracySummary,
    n_sites: int,
    scoring_mode: SignalomeClusteringScoringMode,
    cluster_tree_backend: SignalomeClusterTreeBackend,
    candidate_scoring_backend: SignalomeCandidateScoringBackend,
    max_exact_cluster_tree_sites: int | None,
    max_full_correlation_sites: int,
) -> tuple[
    dict[int, SignalomeClusterCandidateScore],
    dict[int, np.ndarray],
    str,
    str,
    bool,
]:
    """Score candidate cluster counts using full or sampled correlation paths."""

    candidate_counts = [int(cluster_count) for cluster_count in candidate_range]
    if not candidate_counts:
        return {}, {}, ""

    cluster_tree = _build_exact_cluster_tree_with_guard(
        clustering_values=clustering_values,
        n_sites=n_sites,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    )
    candidate_labels = build_cluster_labels_from_tree(
        cluster_tree=cluster_tree,
        cluster_counts=candidate_counts,
    )
    exact_cluster_tree_built = n_sites > 1

    if candidate_scoring_backend == SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL:
        if n_sites > int(max_full_correlation_sites):
            raise SignalomeScaleError(
                "Signalome full candidate-correlation scoring received "
                f"{n_sites:,} sites, which exceeds max_full_correlation_sites="
                f"{int(max_full_correlation_sites):,}. "
                "Use candidate_scoring_backend='sampled', reduce interpreted "
                "sites, or increase max_full_correlation_sites deliberately."
            )
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
        return (
            candidate_scores,
            candidate_labels,
            "",
            SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
            exact_cluster_tree_built,
        )

    candidate_scores = {}
    for cluster_count in candidate_counts:
        labels = candidate_labels[cluster_count]
        cluster_medians = [
            cluster_median_correlation_approximate(
                scoring_values=correlation_values,
                labels=labels,
                label=int(label),
                max_sites_per_cluster=MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
            )
            for label in np.unique(labels)
        ]
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
            "sampled."
        )
    else:
        approximation_note = (
            " Used sampled within-cluster correlation estimates (seeded, "
            "order-invariant sampling) to avoid materializing a full site-by-site "
            "correlation matrix."
        )
    return (
        candidate_scores,
        candidate_labels,
        approximation_note,
        SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
        exact_cluster_tree_built,
    )


def _resolve_candidate_scoring_backend(
    *,
    scoring_mode: SignalomeClusteringScoringMode,
    candidate_scoring_backend: SignalomeCandidateScoringBackend | None,
    n_sites: int,
    max_full_correlation_sites: int,
) -> SignalomeCandidateScoringBackend:
    if candidate_scoring_backend is not None:
        if (
            scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_EXACT
            and candidate_scoring_backend != SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL
        ):
            raise ValueError(
                "scoring_mode='exact' cannot be combined with "
                "candidate_scoring_backend='sampled'"
            )
        if (
            scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE
            and candidate_scoring_backend != SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED
        ):
            raise ValueError(
                "scoring_mode='approximate' cannot be combined with "
                "candidate_scoring_backend='full'"
            )
        return candidate_scoring_backend

    if scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_EXACT:
        return SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL
    if scoring_mode == SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE:
        return SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED
    if n_sites <= int(max_full_correlation_sites):
        return SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL
    return SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED


def _build_exact_cluster_tree_with_guard(
    *,
    clustering_values: np.ndarray,
    n_sites: int,
    cluster_tree_backend: SignalomeClusterTreeBackend,
    candidate_scoring_backend: SignalomeCandidateScoringBackend,
    max_exact_cluster_tree_sites: int | None,
) -> _WardClusterTree:
    if cluster_tree_backend != SIGNALOME_CLUSTER_TREE_BACKEND_EXACT:
        raise ValueError("cluster_tree_backend must be 'exact'")
    if max_exact_cluster_tree_sites is not None and n_sites > int(
        max_exact_cluster_tree_sites
    ):
        raise SignalomeScaleError(
            "Signalome exact cluster-tree construction received "
            f"{n_sites:,} sites, which exceeds max_exact_cluster_tree_sites="
            f"{int(max_exact_cluster_tree_sites):,} "
            "(cluster_tree_backend='exact'). "
            f"candidate_scoring_backend='{candidate_scoring_backend}' still "
            "requires exact cluster-tree construction in the current "
            "implementation."
        )
    return build_cluster_tree(clustering_values)


def _select_best_candidate_count(
    candidate_scores: dict[int, float],
) -> int:
    return max(candidate_scores.items(), key=lambda item: (item[1], -item[0]))[0]


def _select_threshold_candidate(
    *,
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    candidate_labels: dict[int, np.ndarray],
    max_clusters: int,
    threshold: float,
    requested_module_count: int | None,
    reason: str,
    profile_degeneracy: _ProfileDegeneracySummary,
    correlation_exclusion_note: str,
    approximation_note: str,
    cluster_tree_backend: str,
    candidate_scoring_mode: str,
    exact_cluster_tree_built: bool,
) -> _ModuleSelectionComputation | None:
    passing_candidates = filter_cluster_candidates(
        candidate_scores,
        threshold=threshold,
    )
    if not passing_candidates:
        return None
    return _build_module_selection_result(
        strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        selected_module_count=_select_best_candidate_count(passing_candidates),
        requested_module_count=requested_module_count,
        threshold_used=float(threshold),
        max_clusters_evaluated=int(max_clusters),
        candidate_scores=candidate_scores,
        reason=reason + correlation_exclusion_note + approximation_note,
        profile_degeneracy=profile_degeneracy,
        excluded_from_correlation_count=profile_degeneracy.excluded_count,
        candidate_labels=candidate_labels,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
    )


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


def filter_cluster_candidates(
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    *,
    threshold: float,
) -> dict[int, float]:
    """Return candidate counts whose minimum cluster median passes threshold."""

    return {
        cluster_count: score.mean_median_correlation
        for cluster_count, score in candidate_scores.items()
        if score.min_median_correlation >= threshold
    }


def build_cluster_tree(scoring_values: np.ndarray) -> _WardClusterTree:
    """Build a Ward agglomerative merge tree."""

    values = np.asarray(scoring_values, dtype=float)
    if values.ndim != 2:
        raise ValueError("scoring_values must be a 2D array")
    n_sites = int(values.shape[0])
    if n_sites <= 1:
        return _WardClusterTree(n_sites=n_sites, merges=())

    centroids: dict[int, np.ndarray] = {
        row: values[row].astype(float, copy=True) for row in range(n_sites)
    }
    sizes: dict[int, int] = {row: 1 for row in range(n_sites)}
    versions: dict[int, int] = {row: 0 for row in range(n_sites)}
    active: set[int] = set(range(n_sites))
    heap: list[tuple[float, int, int, int, int]] = []

    for left in range(n_sites - 1):
        for right in range(left + 1, n_sites):
            distance = _ward_distance(
                centroids[left],
                sizes[left],
                centroids[right],
                sizes[right],
            )
            heapq.heappush(
                heap,
                (
                    float(distance),
                    left,
                    right,
                    versions[left],
                    versions[right],
                ),
            )

    merges: list[tuple[int, int]] = []
    next_cluster_id = n_sites
    while len(active) > 1:
        left, right = _pop_next_valid_merge(
            heap=heap,
            active=active,
            versions=versions,
        )
        merges.append((left, right))

        left_size = sizes.pop(left)
        right_size = sizes.pop(right)
        left_centroid = centroids.pop(left)
        right_centroid = centroids.pop(right)
        active.remove(left)
        active.remove(right)

        merged_size = left_size + right_size
        merged_centroid = (
            float(left_size) * left_centroid + float(right_size) * right_centroid
        ) / float(merged_size)
        merged_id = next_cluster_id
        next_cluster_id += 1

        sizes[merged_id] = merged_size
        centroids[merged_id] = merged_centroid
        versions[merged_id] = 0
        active.add(merged_id)

        for other in active:
            if other == merged_id:
                continue
            distance = _ward_distance(
                centroids[other],
                sizes[other],
                centroids[merged_id],
                sizes[merged_id],
            )
            first = min(other, merged_id)
            second = max(other, merged_id)
            heapq.heappush(
                heap,
                (
                    float(distance),
                    first,
                    second,
                    versions[first],
                    versions[second],
                ),
            )

    return _WardClusterTree(n_sites=n_sites, merges=tuple(merges))


def build_cluster_labels_from_tree(
    *,
    cluster_tree: _WardClusterTree,
    cluster_counts: Iterable[int],
) -> dict[int, np.ndarray]:
    """Cut a Ward tree into labels for one or more cluster counts."""

    requested_counts = [int(cluster_count) for cluster_count in cluster_counts]
    if not requested_counts:
        return {}
    n_sites = int(cluster_tree.n_sites)
    unique_counts = sorted(
        {max(1, min(int(count), n_sites)) for count in requested_counts}
    )
    labels_by_count: dict[int, np.ndarray] = {}
    if n_sites == 0:
        for count in unique_counts:
            labels_by_count[count] = np.zeros(0, dtype=int)
        return {
            count: labels_by_count[max(1, min(int(count), n_sites))]
            for count in requested_counts
        }

    if n_sites == 1:
        base = np.zeros(1, dtype=int)
        for count in unique_counts:
            labels_by_count[count] = base.copy()
        return {
            count: labels_by_count[max(1, min(int(count), n_sites))]
            for count in requested_counts
        }

    current_members: dict[int, np.ndarray] = {
        site_id: np.asarray([site_id], dtype=int) for site_id in range(n_sites)
    }
    current_cluster_count = n_sites
    next_cluster_id = n_sites

    if current_cluster_count in unique_counts:
        labels_by_count[current_cluster_count] = _labels_from_members(
            members=current_members,
            n_sites=n_sites,
        )

    for left, right in cluster_tree.merges:
        left_members = current_members.pop(left)
        right_members = current_members.pop(right)
        merged = np.concatenate([left_members, right_members])
        current_members[next_cluster_id] = merged
        next_cluster_id += 1
        current_cluster_count -= 1
        if current_cluster_count in unique_counts:
            labels_by_count[current_cluster_count] = _labels_from_members(
                members=current_members,
                n_sites=n_sites,
            )
        if current_cluster_count <= min(unique_counts):
            break

    return {
        count: labels_by_count[max(1, min(int(count), n_sites))].copy()
        for count in requested_counts
    }


def fit_cluster_labels(
    scoring_values: np.ndarray,
    cluster_count: int,
    *,
    cluster_tree_backend: SignalomeClusterTreeBackend = (
        SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
    ),
    candidate_scoring_backend: SignalomeCandidateScoringBackend = (
        SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL
    ),
    max_exact_cluster_tree_sites: int | None = None,
) -> np.ndarray:
    """Fit Ward clustering and return 0-indexed labels for one count."""

    values = np.asarray(scoring_values, dtype=float)
    n_sites = int(values.shape[0])
    resolved_cluster_count = max(1, min(int(cluster_count), n_sites))
    if resolved_cluster_count == 1:
        return np.zeros(n_sites, dtype=int)
    tree = _build_exact_cluster_tree_with_guard(
        clustering_values=values,
        n_sites=n_sites,
        cluster_tree_backend=cluster_tree_backend,
        candidate_scoring_backend=candidate_scoring_backend,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    )
    return build_cluster_labels_from_tree(
        cluster_tree=tree,
        cluster_counts=[resolved_cluster_count],
    )[resolved_cluster_count].astype(int, copy=False)


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

    cluster_positions = np.flatnonzero(labels == label)
    if cluster_positions.size <= 1:
        return 0.0

    if cluster_positions.size > max_sites_per_cluster:
        cluster_positions = _sample_cluster_positions_for_approximation(
            scoring_values=scoring_values,
            cluster_positions=cluster_positions,
            sample_size=max_sites_per_cluster,
        )

    cluster_values = np.asarray(scoring_values, dtype=float)[cluster_positions]
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


def _prepare_scoring_values_for_clustering(scoring_values: np.ndarray) -> np.ndarray:
    values = np.asarray(scoring_values, dtype=float).copy()
    if values.size == 0:
        return values
    values[~np.isfinite(values)] = np.nan
    column_medians = np.nanmedian(values, axis=0)
    column_medians = np.where(np.isfinite(column_medians), column_medians, 0.0)
    row_positions, column_positions = np.where(np.isnan(values))
    if row_positions.size > 0:
        values[row_positions, column_positions] = column_medians[column_positions]
    return values


def _validate_threshold(value: float, *, field_name: str) -> None:
    if not np.isfinite(float(value)):
        raise ValueError(f"{field_name} must be finite")
    if float(value) < 0.0 or float(value) > 1.0:
        raise ValueError(f"{field_name} must be between 0 and 1")


def _ward_distance(
    left_centroid: np.ndarray,
    left_size: int,
    right_centroid: np.ndarray,
    right_size: int,
) -> float:
    delta = np.asarray(left_centroid, dtype=float) - np.asarray(
        right_centroid, dtype=float
    )
    squared_norm = float(np.dot(delta, delta))
    return (
        (float(left_size) * float(right_size))
        / float(left_size + right_size)
        * squared_norm
    )


def _pop_next_valid_merge(
    *,
    heap: list[tuple[float, int, int, int, int]],
    active: set[int],
    versions: dict[int, int],
) -> tuple[int, int]:
    while heap:
        _, left, right, left_version, right_version = heapq.heappop(heap)
        if left not in active or right not in active:
            continue
        if versions[left] != left_version or versions[right] != right_version:
            continue
        return left, right
    raise RuntimeError("failed to resolve a valid merge from ward clustering heap")


def _labels_from_members(
    *,
    members: dict[int, np.ndarray],
    n_sites: int,
) -> np.ndarray:
    labels = np.zeros(n_sites, dtype=int)
    sorted_clusters = sorted(members.items(), key=lambda item: int(item[0]))
    for label, (_, cluster_members) in enumerate(sorted_clusters):
        labels[np.asarray(cluster_members, dtype=int)] = int(label)
    return labels


__all__ = [
    "ClusterSitesResult",
    "MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER",
    "MAX_FULL_CORRELATION_SITE_COUNT",
    "NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE",
    "SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL",
    "SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED",
    "SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED",
    "SIGNALOME_CLUSTER_TREE_BACKEND_EXACT",
    "SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE",
    "SIGNALOME_CLUSTERING_SCORING_MODE_AUTO",
    "SIGNALOME_CLUSTERING_SCORING_MODE_EXACT",
    "SignalomeCandidateScoringBackend",
    "SignalomeClusterTreeBackend",
    "SignalomeClusteringScoringMode",
    "build_correlation_exclusion_note",
    "build_correlation_matrix_with_exclusions",
    "build_cluster_labels_from_tree",
    "build_cluster_tree",
    "cluster_median_correlation",
    "cluster_median_correlation_approximate",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "derive_protein_modules",
    "filter_cluster_candidates",
    "fit_cluster_labels",
    "select_module_count",
    "select_module_count_with_diagnostics",
    "summarize_profile_degeneracy",
]
