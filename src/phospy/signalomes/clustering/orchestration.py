"""Signalome clustering and module-count selection helpers."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal, cast

import numpy as np
import pandas as pd

from phospy.errors.workflows import SignalomeScaleError
from phospy.signalomes.clustering.backends.exact_python import (
    ExactWardClusterTree,
)
from phospy.signalomes.clustering.backends.exact_python import (
    build_cluster_labels_from_tree as build_exact_cluster_labels_from_tree,
)
from phospy.signalomes.clustering.backends.exact_python import (
    build_cluster_tree as build_exact_cluster_tree,
)
from phospy.signalomes.clustering.contracts import ClusterTreeEngine
from phospy.signalomes.clustering.diagnostics import (
    approximation_used_from_candidate_mode,
    build_candidate_scoring_sampling_provenance,
    build_module_selection_diagnostics,
)
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION,
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.signalomes.clustering.validation import (
    validate_cluster_count_for_site_count,
    validate_requested_module_count,
)
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
SIGNALOME_TREE_ENGINE_EXACT = "exact"
SignalomeTreeEngine = Literal["exact"]
SIGNALOME_CANDIDATE_SCORING_POLICY_FULL = "full"
SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED = "sampled"
SignalomeCandidateScoringPolicy = Literal["full", "sampled"]
SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED = "not_evaluated"
_CandidateScoringMode = SignalomeCandidateScoringPolicy | Literal["not_evaluated"]
SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD = (
    "deterministic_uniform_without_replacement"
)
SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY = (
    "order_invariant_seed_from_row_hashes_and_sample_size"
)
SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT = "explicit_module_count"
SIGNALOME_CANDIDATE_SCORING_APPLIES_TO = "candidate_module_count_evaluation_only"
SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE = "exact_cluster_tree"
SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE = "single_module_assignment"


@dataclass(frozen=True, slots=True)
class ClusterSitesResult:
    """Cluster labels and diagnostics for module-count selection."""

    site_clusters: pd.Series
    module_selection_diagnostics: SignalomeModuleSelectionDiagnostics
    tree_engine: str = SIGNALOME_TREE_ENGINE_EXACT
    candidate_scoring_mode: _CandidateScoringMode = (
        SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    )
    exact_cluster_tree_built: bool = False
    candidate_scoring_sampling: dict[str, object] | None = None
    candidate_scoring_evaluated: bool = False
    candidate_scoring_skip_reason: str | None = None
    backend_name: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    backend_version: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION
    approximation_used: bool = False
    backend_diagnostics: dict[str, object] | None = None


@dataclass(frozen=True, slots=True)
class _ModuleSelectionComputation:
    diagnostics: SignalomeModuleSelectionDiagnostics
    candidate_labels: dict[int, np.ndarray]
    tree_engine: str
    candidate_scoring_mode: _CandidateScoringMode
    exact_cluster_tree_built: bool
    candidate_scoring_evaluated: bool
    candidate_scoring_skip_reason: str | None
    candidate_scoring_sampling: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class _CandidateClusterScoreResult:
    candidate_scores: dict[int, SignalomeClusterCandidateScore]
    candidate_labels: dict[int, np.ndarray]
    approximation_note: str
    candidate_scoring_mode: _CandidateScoringMode
    exact_cluster_tree_built: bool
    candidate_scoring_evaluated: bool
    candidate_scoring_skip_reason: str | None
    candidate_scoring_sampling: dict[str, object] | None


@dataclass(frozen=True, slots=True)
class _ProfileDegeneracySummary:
    zero_variance_count: int
    near_constant_count: int
    excluded_count: int
    excluded_mask: np.ndarray


@dataclass(frozen=True, slots=True)
class _ClusterTreeOperationsAdapter:
    """Compatibility adapter for legacy cluster-tree operation hooks."""

    engine: ClusterTreeEngine

    def build_cluster_tree(self, scoring_values: np.ndarray) -> object:
        return self.engine.build_tree(scoring_values)

    def build_cluster_labels_from_tree(
        self,
        *,
        cluster_tree: object,
        cluster_counts: Iterable[int],
    ) -> dict[int, np.ndarray]:
        return self.engine.labels_for_counts(tree=cluster_tree, counts=cluster_counts)


@dataclass(frozen=True, slots=True)
class _ExactWardClusterTreeOperations:
    """Default exact-tree operations with legacy monkeypatch hook behavior."""

    def build_cluster_tree(self, scoring_values: np.ndarray) -> object:
        # Keep a runtime indirection through the exact-python module so
        # performance-contract monkeypatching of exact-python helpers remains
        # effective after orchestration extraction.
        from phospy.signalomes.clustering import exact_python as exact_clustering

        return exact_clustering._build_cluster_tree(scoring_values)

    def build_cluster_labels_from_tree(
        self,
        *,
        cluster_tree: object,
        cluster_counts: Iterable[int],
    ) -> dict[int, np.ndarray]:
        from phospy.signalomes.clustering import exact_python as exact_clustering

        return exact_clustering.build_cluster_labels_from_tree(
            cluster_tree=cast(ExactWardClusterTree, cluster_tree),
            cluster_counts=cluster_counts,
        )


ClusterTreeOperations = _ClusterTreeOperationsAdapter


_EXACT_WARD_CLUSTER_TREE_OPERATIONS = _ExactWardClusterTreeOperations()


def _resolve_cluster_tree_operations(
    cluster_tree_operations: ClusterTreeOperations | None,
) -> ClusterTreeOperations | _ExactWardClusterTreeOperations:
    if cluster_tree_operations is None:
        return _EXACT_WARD_CLUSTER_TREE_OPERATIONS
    return cluster_tree_operations


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
    tree_engine: SignalomeTreeEngine = (SIGNALOME_TREE_ENGINE_EXACT),
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> pd.Series:
    """Cluster phosphosites into site clusters."""

    return cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        cluster_tree_operations=cluster_tree_operations,
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
    tree_engine: SignalomeTreeEngine = (SIGNALOME_TREE_ENGINE_EXACT),
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> ClusterSitesResult:
    """Cluster phosphosites and capture module-selection diagnostics.

    Requested module counts are validated strictly against available clustering
    sites and are never clamped or rewritten.
    """

    scoring_values = np.asarray(scoring_matrix.to_numpy(dtype=float, copy=False))
    selection = _compute_module_selection(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        cluster_tree_operations=cluster_tree_operations,
    )
    diagnostics = selection.diagnostics
    n_sites = int(scoring_values.shape[0])
    module_count = validate_cluster_count_for_site_count(
        cluster_count=int(diagnostics.selected_module_count),
        available_clustering_site_count=n_sites,
        field_name="selected_module_count",
    )
    exact_cluster_tree_built = bool(selection.exact_cluster_tree_built)

    if module_count == 1:
        labels = np.ones(n_sites, dtype=int)
    else:
        cached_labels = selection.candidate_labels.get(module_count)
        if cached_labels is not None:
            labels = cached_labels.astype(int, copy=False) + 1
            exact_cluster_tree_built = True
        else:
            labels = (
                fit_cluster_labels(
                    scoring_values=_prepare_scoring_values_for_clustering(
                        scoring_values
                    ),
                    cluster_count=module_count,
                    tree_engine=tree_engine,
                    candidate_scoring_policy=(
                        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
                        if candidate_scoring_policy is None
                        else candidate_scoring_policy
                    ),
                    max_exact_tree_sites=max_exact_tree_sites,
                    cluster_tree_operations=cluster_tree_operations,
                )
                + 1
            )
            exact_cluster_tree_built = True

    return ClusterSitesResult(
        site_clusters=pd.Series(
            labels,
            index=scoring_matrix.index.copy(),
            dtype=int,
            name=SITE_CLUSTER_COLUMN,
        ),
        module_selection_diagnostics=diagnostics,
        tree_engine=selection.tree_engine,
        candidate_scoring_mode=selection.candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        candidate_scoring_evaluated=selection.candidate_scoring_evaluated,
        candidate_scoring_skip_reason=selection.candidate_scoring_skip_reason,
        candidate_scoring_sampling=selection.candidate_scoring_sampling,
        approximation_used=approximation_used_from_candidate_mode(
            candidate_scoring_mode=str(selection.candidate_scoring_mode),
            candidate_scoring_evaluated=bool(selection.candidate_scoring_evaluated),
        ),
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
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
) -> int:
    """Select a module count from a scoring matrix."""

    return select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
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
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
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
        max_exact_tree_sites=max_exact_tree_sites,
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
    tree_engine: SignalomeTreeEngine = (SIGNALOME_TREE_ENGINE_EXACT),
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
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
    if tree_engine != SIGNALOME_TREE_ENGINE_EXACT:
        raise ValueError("tree_engine must be 'exact'")
    if candidate_scoring_policy not in {
        None,
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    }:
        raise ValueError("candidate_scoring_policy must be one of: full, sampled")
    if max_full_candidate_scoring_sites < 1:
        raise ValueError("max_full_candidate_scoring_sites must be >= 1")
    resolved_max_exact_tree_sites = _resolve_max_exact_tree_sites(max_exact_tree_sites)

    scoring_array = np.asarray(scoring_values, dtype=float)
    if scoring_array.ndim != 2:
        raise ValueError("scoring_values must be a 2D array")
    n_sites = int(scoring_array.shape[0])
    requested_module_count = validate_requested_module_count(
        requested_module_count=requested_module_count,
        available_clustering_site_count=n_sites,
        field_name="signalome workflow request config.module_count",
    )
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
    resolved_candidate_scoring_policy = _resolve_candidate_scoring_policy(
        scoring_mode=scoring_mode,
        candidate_scoring_policy=candidate_scoring_policy,
        n_sites=n_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
    )

    clustering_values = _prepare_scoring_values_for_clustering(scoring_array)
    candidate_score_result = _compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=scoring_array,
        candidate_range=range(2, resolved_max_clusters + 1),
        profile_degeneracy=profile_degeneracy,
        n_sites=n_sites,
        scoring_mode=scoring_mode,
        tree_engine=tree_engine,
        candidate_scoring_policy=resolved_candidate_scoring_policy,
        max_exact_tree_sites=resolved_max_exact_tree_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        cluster_tree_operations=cluster_tree_operations,
    )

    primary_selection = _select_threshold_candidate(
        candidate_scores=candidate_score_result.candidate_scores,
        candidate_labels=candidate_score_result.candidate_labels,
        max_clusters=resolved_max_clusters,
        threshold=primary_threshold,
        requested_module_count=requested_module_count,
        reason=(
            "selected the highest-scoring candidate that satisfied the primary "
            "within-cluster correlation threshold"
        ),
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
        approximation_note=candidate_score_result.approximation_note,
        candidate_scoring_mode=candidate_score_result.candidate_scoring_mode,
        exact_cluster_tree_built=candidate_score_result.exact_cluster_tree_built,
        candidate_scoring_evaluated=candidate_score_result.candidate_scoring_evaluated,
        candidate_scoring_skip_reason=candidate_score_result.candidate_scoring_skip_reason,
        tree_engine=tree_engine,
        candidate_scoring_sampling=candidate_score_result.candidate_scoring_sampling,
    )
    if primary_selection is not None:
        return primary_selection

    fallback_selection = _select_threshold_candidate(
        candidate_scores=candidate_score_result.candidate_scores,
        candidate_labels=candidate_score_result.candidate_labels,
        max_clusters=resolved_max_clusters,
        threshold=fallback_threshold,
        requested_module_count=requested_module_count,
        reason=(
            "no candidate satisfied the primary threshold; selected the "
            "highest-scoring fallback candidate"
        ),
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
        approximation_note=candidate_score_result.approximation_note,
        candidate_scoring_mode=candidate_score_result.candidate_scoring_mode,
        exact_cluster_tree_built=candidate_score_result.exact_cluster_tree_built,
        candidate_scoring_evaluated=candidate_score_result.candidate_scoring_evaluated,
        candidate_scoring_skip_reason=candidate_score_result.candidate_scoring_skip_reason,
        tree_engine=tree_engine,
        candidate_scoring_sampling=candidate_score_result.candidate_scoring_sampling,
    )
    if fallback_selection is not None:
        return fallback_selection

    return _build_module_selection_result(
        strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        selected_module_count=1,
        requested_module_count=requested_module_count,
        threshold_used=None,
        max_clusters_evaluated=resolved_max_clusters,
        candidate_scores=candidate_score_result.candidate_scores,
        reason=(
            "no candidate module count satisfied the configured correlation "
            "thresholds, so the workflow fell back to one module"
        )
        + correlation_exclusion_note
        + candidate_score_result.approximation_note,
        profile_degeneracy=profile_degeneracy,
        excluded_from_correlation_count=profile_degeneracy.excluded_count,
        candidate_labels=candidate_score_result.candidate_labels,
        candidate_scoring_mode=candidate_score_result.candidate_scoring_mode,
        exact_cluster_tree_built=candidate_score_result.exact_cluster_tree_built,
        candidate_scoring_evaluated=candidate_score_result.candidate_scoring_evaluated,
        candidate_scoring_skip_reason=candidate_score_result.candidate_scoring_skip_reason,
        tree_engine=tree_engine,
        candidate_scoring_sampling=candidate_score_result.candidate_scoring_sampling,
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
    tree_engine: str = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_mode: _CandidateScoringMode = (
        SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED
    ),
    exact_cluster_tree_built: bool = False,
    candidate_scoring_evaluated: bool = False,
    candidate_scoring_skip_reason: str | None = None,
    candidate_scoring_sampling: dict[str, object] | None = None,
) -> _ModuleSelectionComputation:
    return _ModuleSelectionComputation(
        diagnostics=build_module_selection_diagnostics(
            strategy=strategy,
            selected_module_count=selected_module_count,
            requested_module_count=requested_module_count,
            threshold_used=threshold_used,
            max_clusters_evaluated=max_clusters_evaluated,
            candidate_scores=candidate_scores,
            reason=reason,
            zero_variance_profile_count=profile_degeneracy.zero_variance_count,
            near_constant_profile_count=profile_degeneracy.near_constant_count,
            excluded_from_correlation_count=excluded_from_correlation_count,
        ),
        candidate_labels=candidate_labels,
        tree_engine=str(tree_engine),
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=bool(exact_cluster_tree_built),
        candidate_scoring_evaluated=bool(candidate_scoring_evaluated),
        candidate_scoring_skip_reason=(
            None
            if candidate_scoring_skip_reason is None
            else str(candidate_scoring_skip_reason)
        ),
        candidate_scoring_sampling=candidate_scoring_sampling,
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
        resolved_count = int(requested_module_count)
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
                candidate_scoring_evaluated=False,
                candidate_scoring_skip_reason=(
                    SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT
                ),
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

    resolved_max_exact_tree_sites = _resolve_max_exact_tree_sites(max_exact_tree_sites)
    # Guard ordering policy:
    # - If full candidate-correlation scoring exceeds max_full_candidate_scoring_sites
    #   while exact-tree construction is still permitted, fail here before any
    #   exact-tree construction is attempted.
    # - If both max_full_candidate_scoring_sites and max_exact_tree_sites are
    #   exceeded, defer to the exact-tree guard below as the canonical first
    #   failure for that configuration.
    if (
        candidate_scoring_policy == SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
        and n_sites > int(max_full_candidate_scoring_sites)
        and n_sites <= int(resolved_max_exact_tree_sites)
    ):
        raise SignalomeScaleError(
            "Signalome full candidate-correlation scoring would evaluate "
            f"{n_sites:,} sites, which exceeds configured "
            f"max_full_candidate_scoring_sites={int(max_full_candidate_scoring_sites):,}. "
            "Exact cluster-tree construction has not been attempted for this "
            "request. Use candidate_scoring_policy='sampled' for candidate "
            "module-count evaluation, reduce interpreted sites, or increase "
            "max_full_candidate_scoring_sites deliberately."
        )

    cluster_tree = _build_exact_cluster_tree_with_guard(
        clustering_values=clustering_values,
        n_sites=n_sites,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=resolved_max_exact_tree_sites,
        cluster_tree_operations=cluster_tree_operations,
    )
    tree_operations = _resolve_cluster_tree_operations(cluster_tree_operations)
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


def _resolve_candidate_scoring_policy(
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


def _resolve_max_exact_tree_sites(
    max_exact_tree_sites: int | None,
) -> int:
    """Resolve exact-tree guard limit; `None` maps to the safe default limit."""

    resolved = (
        MAX_FULL_CORRELATION_SITE_COUNT
        if max_exact_tree_sites is None
        else int(max_exact_tree_sites)
    )
    if resolved < 1:
        raise ValueError("max_exact_tree_sites must be >= 1")
    return resolved


def _build_exact_cluster_tree_with_guard(
    *,
    clustering_values: np.ndarray,
    n_sites: int,
    tree_engine: SignalomeTreeEngine,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy,
    max_exact_tree_sites: int | None,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> object:
    if tree_engine != SIGNALOME_TREE_ENGINE_EXACT:
        raise ValueError("tree_engine must be 'exact'")
    resolved_max_exact_tree_sites = _resolve_max_exact_tree_sites(max_exact_tree_sites)
    if n_sites > int(resolved_max_exact_tree_sites):
        raise SignalomeScaleError(
            "Signalome exact cluster-tree construction received "
            f"{n_sites:,} sites, which exceeds max_exact_tree_sites="
            f"{int(resolved_max_exact_tree_sites):,} "
            "(tree_engine='exact'). "
            f"candidate_scoring_policy='{candidate_scoring_policy}' still "
            "requires exact cluster-tree construction in the current "
            "implementation."
        )
    tree_operations = _resolve_cluster_tree_operations(cluster_tree_operations)
    return tree_operations.build_cluster_tree(clustering_values)


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
    tree_engine: str,
    candidate_scoring_mode: _CandidateScoringMode,
    exact_cluster_tree_built: bool,
    candidate_scoring_evaluated: bool,
    candidate_scoring_skip_reason: str | None,
    candidate_scoring_sampling: dict[str, object] | None,
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
        tree_engine=tree_engine,
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        candidate_scoring_evaluated=candidate_scoring_evaluated,
        candidate_scoring_skip_reason=candidate_scoring_skip_reason,
        candidate_scoring_sampling=candidate_scoring_sampling,
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


def _build_cluster_tree(scoring_values: np.ndarray) -> ExactWardClusterTree:
    """Compatibility wrapper for pure-Python exact Ward tree construction."""

    return build_exact_cluster_tree(scoring_values)


def build_cluster_labels_from_tree(
    *,
    cluster_tree: ExactWardClusterTree,
    cluster_counts: Iterable[int],
) -> dict[int, np.ndarray]:
    """Compatibility wrapper for pure-Python exact Ward tree cutting."""

    return build_exact_cluster_labels_from_tree(
        cluster_tree=cluster_tree,
        cluster_counts=cluster_counts,
    )


def fit_cluster_labels(
    scoring_values: np.ndarray,
    cluster_count: int,
    *,
    tree_engine: SignalomeTreeEngine = (SIGNALOME_TREE_ENGINE_EXACT),
    candidate_scoring_policy: SignalomeCandidateScoringPolicy = (
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    ),
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> np.ndarray:
    """Fit Ward clustering and return 0-indexed labels for one count."""

    values = np.asarray(scoring_values, dtype=float)
    n_sites = int(values.shape[0])
    resolved_cluster_count = validate_cluster_count_for_site_count(
        cluster_count=int(cluster_count),
        available_clustering_site_count=n_sites,
        field_name="cluster_count",
    )
    if resolved_cluster_count == 1:
        return np.zeros(n_sites, dtype=int)
    resolved_max_exact_tree_sites = _resolve_max_exact_tree_sites(max_exact_tree_sites)
    tree = _build_exact_cluster_tree_with_guard(
        clustering_values=values,
        n_sites=n_sites,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=resolved_max_exact_tree_sites,
        cluster_tree_operations=cluster_tree_operations,
    )
    tree_operations = _resolve_cluster_tree_operations(cluster_tree_operations)
    return tree_operations.build_cluster_labels_from_tree(
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
) -> dict[str, object]:
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


def run_clustering_with_tree_engine(
    *,
    request: SignalomeClusteringEngineRequest,
    tree_engine: ClusterTreeEngine,
    clustering_engine: str,
    backend_version: str,
    backend_diagnostics: dict[str, object],
) -> SignalomeClusteringEngineResult:
    """Run shared orchestration with an injected tree engine implementation."""

    requested_tree_engine = request.tree_engine
    if requested_tree_engine is not None:
        resolved_requested_tree_engine = str(requested_tree_engine)
        if resolved_requested_tree_engine not in {
            str(tree_engine.name),
            SIGNALOME_TREE_ENGINE_EXACT,
            "exact_python",
            "scipy_hierarchical",
        }:
            raise ValueError(
                f"unsupported tree_engine request {resolved_requested_tree_engine!r}"
            )

    candidate_scoring_policy = (
        None
        if request.candidate_scoring_policy is None
        else str(request.candidate_scoring_policy)
    )

    clustering_result = cluster_sites_with_diagnostics(
        scoring_matrix=request.scoring_matrix,
        requested_module_count=request.requested_module_count,
        primary_threshold=request.primary_threshold,
        fallback_threshold=request.fallback_threshold,
        max_clusters=request.max_clusters,
        tree_engine=SIGNALOME_TREE_ENGINE_EXACT,
        candidate_scoring_policy=candidate_scoring_policy,  # type: ignore[arg-type]
        max_exact_tree_sites=request.max_exact_tree_sites,
        max_full_candidate_scoring_sites=request.max_full_candidate_scoring_sites,
        cluster_tree_operations=_ClusterTreeOperationsAdapter(engine=tree_engine),
    )
    protein_modules = derive_protein_modules(
        site_clusters=clustering_result.site_clusters,
        site_to_protein=request.site_to_protein,
    )
    selected_module_count = int(
        clustering_result.module_selection_diagnostics.selected_module_count
    )
    resolved_backend_diagnostics = {
        "backend_name": str(clustering_engine),
        "tree_engine": str(tree_engine.name),
        "tree_engine_version": str(tree_engine.version),
        **backend_diagnostics,
        "selected_module_count": selected_module_count,
        "input_site_count": int(request.scoring_matrix.shape[0]),
        "exact_tree_path_used": bool(clustering_result.exact_cluster_tree_built),
    }
    return SignalomeClusteringEngineResult(
        site_clusters=clustering_result.site_clusters,
        protein_modules=protein_modules,
        selected_module_count=selected_module_count,
        module_selection_diagnostics=clustering_result.module_selection_diagnostics,
        backend_name=str(clustering_engine),
        backend_version=str(backend_version),
        approximation_used=bool(clustering_result.approximation_used),
        exact_cluster_tree_built=bool(clustering_result.exact_cluster_tree_built),
        tree_engine=str(clustering_result.tree_engine),
        candidate_scoring_mode=str(clustering_result.candidate_scoring_mode),
        candidate_scoring_evaluated=bool(clustering_result.candidate_scoring_evaluated),
        candidate_scoring_skip_reason=(
            None
            if clustering_result.candidate_scoring_skip_reason is None
            else str(clustering_result.candidate_scoring_skip_reason)
        ),
        candidate_scoring_sampling=clustering_result.candidate_scoring_sampling,
        backend_diagnostics=resolved_backend_diagnostics,
        threshold_metadata={
            "primary_threshold": float(request.primary_threshold),
            "fallback_threshold": float(request.fallback_threshold),
        },
        limit_metadata={
            "max_exact_tree_sites": (
                None
                if request.max_exact_tree_sites is None
                else int(request.max_exact_tree_sites)
            ),
            "max_full_candidate_scoring_sites": int(
                request.max_full_candidate_scoring_sites
            ),
            "max_clusters": int(request.max_clusters),
        },
    )


__all__ = [
    "ClusterTreeOperations",
    "ClusterSitesResult",
    "MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER",
    "MAX_FULL_CORRELATION_SITE_COUNT",
    "NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_FULL",
    "SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED",
    "SIGNALOME_CANDIDATE_SCORING_APPLIES_TO",
    "SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED",
    "SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT",
    "SIGNALOME_TREE_ENGINE_EXACT",
    "SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE",
    "SIGNALOME_CLUSTERING_SCORING_MODE_AUTO",
    "SIGNALOME_CLUSTERING_SCORING_MODE_EXACT",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE",
    "SignalomeCandidateScoringPolicy",
    "SignalomeTreeEngine",
    "SignalomeClusteringScoringMode",
    "_CandidateClusterScoreResult",
    "_CandidateScoringMode",
    "_ProfileDegeneracySummary",
    "_build_cluster_tree",
    "_compute_candidate_cluster_scores",
    "_prepare_scoring_values_for_clustering",
    "_resolve_pre_scoring_module_selection",
    "build_correlation_exclusion_note",
    "build_correlation_matrix_with_exclusions",
    "build_cluster_labels_from_tree",
    "cluster_median_correlation",
    "cluster_median_correlation_approximate",
    "cluster_sites",
    "cluster_sites_with_diagnostics",
    "derive_protein_modules",
    "filter_cluster_candidates",
    "fit_cluster_labels",
    "run_clustering_with_tree_engine",
    "select_module_count",
    "select_module_count_with_diagnostics",
    "summarize_profile_degeneracy",
]
