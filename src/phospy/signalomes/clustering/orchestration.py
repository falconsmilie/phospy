"""Signalome clustering orchestration coordination layer."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd

from phospy.signalomes.clustering.candidate_scoring import (
    _CandidateClusterScoreResult,
    _ProfileDegeneracySummary,
)
from phospy.signalomes.clustering.candidate_scoring import (
    build_correlation_exclusion_note as _build_correlation_exclusion_note_impl,
)
from phospy.signalomes.clustering.candidate_scoring import (
    build_correlation_matrix_with_exclusions as _build_correlation_matrix_with_exclusions_impl,
)
from phospy.signalomes.clustering.candidate_scoring import (
    cluster_median_correlation as _cluster_median_correlation_impl,
)
from phospy.signalomes.clustering.candidate_scoring import (
    cluster_median_correlation_approximate as _cluster_median_correlation_approximate_impl,
)
from phospy.signalomes.clustering.candidate_scoring import (
    compute_candidate_cluster_scores as _compute_candidate_cluster_scores_impl,
)
from phospy.signalomes.clustering.candidate_scoring import (
    resolve_candidate_scoring_policy as _resolve_candidate_scoring_policy_impl,
)
from phospy.signalomes.clustering.candidate_scoring import (
    summarize_profile_degeneracy as _summarize_profile_degeneracy_impl,
)
from phospy.signalomes.clustering.candidate_selection import (
    _ModuleSelectionComputation,
)
from phospy.signalomes.clustering.candidate_selection import (
    build_module_selection_result as _build_module_selection_result_impl,
)
from phospy.signalomes.clustering.candidate_selection import (
    filter_cluster_candidates as _filter_cluster_candidates_impl,
)
from phospy.signalomes.clustering.candidate_selection import (
    resolve_pre_scoring_module_selection as _resolve_pre_scoring_module_selection_impl,
)
from phospy.signalomes.clustering.candidate_selection import (
    select_best_candidate_count as _select_best_candidate_count_impl,
)
from phospy.signalomes.clustering.candidate_selection import (
    select_threshold_candidate as _select_threshold_candidate_impl,
)
from phospy.signalomes.clustering.contracts import ClusterTreeEngine
from phospy.signalomes.clustering.diagnostics import (
    approximation_used_from_candidate_mode,
)
from phospy.signalomes.clustering.models import (
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION,
    SignalomeClusteringEngineRequest,
    SignalomeClusteringEngineResult,
)
from phospy.signalomes.clustering.policies import (
    MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
    MAX_FULL_CORRELATION_SITE_COUNT,
    NEAR_CONSTANT_PROFILE_VARIANCE_TOLERANCE,
    SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
    SIGNALOME_CANDIDATE_SCORING_MODE_NOT_EVALUATED,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY,
    SIGNALOME_CANDIDATE_SCORING_SKIP_REASON_EXPLICIT_MODULE_COUNT,
    SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE,
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_CLUSTERING_SCORING_MODE_EXACT,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
    _CandidateScoringMode,
)
from phospy.signalomes.clustering.protein_modules import derive_protein_modules
from phospy.signalomes.clustering.scale_guards import (
    resolve_max_exact_tree_sites as _resolve_max_exact_tree_sites_impl,
)
from phospy.signalomes.clustering.tree_building import (
    ClusterTreeOperations,
    ClusterTreeOperationsAdapter,
)
from phospy.signalomes.clustering.tree_building import (
    build_cluster_labels_from_tree as _build_cluster_labels_from_tree_impl,
)
from phospy.signalomes.clustering.tree_building import (
    build_cluster_tree as _build_cluster_tree_impl,
)
from phospy.signalomes.clustering.tree_building import (
    build_exact_cluster_tree_with_guard as _build_exact_cluster_tree_with_guard_impl,
)
from phospy.signalomes.clustering.tree_building import (
    fit_cluster_labels as _fit_cluster_labels_impl,
)
from phospy.signalomes.clustering.tree_building import (
    prepare_scoring_values_for_clustering as _prepare_scoring_values_for_clustering_impl,
)
from phospy.signalomes.clustering.tree_building import (
    resolve_cluster_tree_operations as _resolve_cluster_tree_operations_impl,
)
from phospy.signalomes.clustering.validation import (
    validate_cluster_count_for_site_count,
    validate_requested_module_count,
)
from phospy.signalomes.constants import SITE_CLUSTER_COLUMN
from phospy.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SignalomeClusterCandidateScore,
    SignalomeModuleSelectionDiagnostics,
)


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


def cluster_sites(
    *,
    scoring_matrix: pd.DataFrame,
    requested_module_count: int | None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
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
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None = None,
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    max_full_candidate_scoring_sites: int = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
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
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
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
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
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


def _compute_module_selection(
    *,
    scoring_values: np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
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
        cluster_tree_operations=ClusterTreeOperationsAdapter(engine=tree_engine),
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
        "tree_generation_mode": "full_exact_tree_construction",
        "tree_generation_is_approximate": False,
        "tree_generation_scope": "module_count_selection_and_final_assignment",
        "candidate_scoring_scope": SIGNALOME_CANDIDATE_SCORING_APPLIES_TO,
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


# Compatibility wrappers retained for legacy imports/monkeypatch hooks.
def _build_module_selection_result(**kwargs: object) -> _ModuleSelectionComputation:
    impl = cast(
        Callable[..., _ModuleSelectionComputation], _build_module_selection_result_impl
    )
    return impl(**kwargs)


def _resolve_pre_scoring_module_selection(**kwargs: object):
    impl = cast(
        Callable[..., tuple[_ModuleSelectionComputation | None, int]],
        _resolve_pre_scoring_module_selection_impl,
    )
    return impl(**kwargs)


def _compute_candidate_cluster_scores(**kwargs: object) -> _CandidateClusterScoreResult:
    impl = cast(
        Callable[..., _CandidateClusterScoreResult],
        _compute_candidate_cluster_scores_impl,
    )
    return impl(**kwargs)


def _resolve_candidate_scoring_policy(
    **kwargs: object,
) -> SignalomeCandidateScoringPolicy:
    impl = cast(
        Callable[..., SignalomeCandidateScoringPolicy],
        _resolve_candidate_scoring_policy_impl,
    )
    return impl(**kwargs)


def _resolve_max_exact_tree_sites(max_exact_tree_sites: int | None) -> int:
    return _resolve_max_exact_tree_sites_impl(max_exact_tree_sites)


def _build_exact_cluster_tree_with_guard(**kwargs: object) -> object:
    impl = cast(Callable[..., object], _build_exact_cluster_tree_with_guard_impl)
    return impl(**kwargs)


def _resolve_cluster_tree_operations(
    cluster_tree_operations: ClusterTreeOperations | None,
) -> object:
    return _resolve_cluster_tree_operations_impl(cluster_tree_operations)


def _select_best_candidate_count(candidate_scores: dict[int, float]) -> int:
    return _select_best_candidate_count_impl(candidate_scores)


def _select_threshold_candidate(**kwargs: object):
    impl = cast(
        Callable[..., _ModuleSelectionComputation | None],
        _select_threshold_candidate_impl,
    )
    return impl(**kwargs)


def summarize_profile_degeneracy(
    scoring_values: np.ndarray,
) -> _ProfileDegeneracySummary:
    return _summarize_profile_degeneracy_impl(scoring_values)


def build_correlation_exclusion_note(summary: _ProfileDegeneracySummary) -> str:
    return _build_correlation_exclusion_note_impl(summary)


def build_correlation_matrix_with_exclusions(
    scoring_values: np.ndarray,
    *,
    excluded_mask: np.ndarray | None = None,
) -> np.ndarray:
    return _build_correlation_matrix_with_exclusions_impl(
        scoring_values=scoring_values,
        excluded_mask=excluded_mask,
    )


def filter_cluster_candidates(
    candidate_scores: dict[int, SignalomeClusterCandidateScore],
    *,
    threshold: float,
) -> dict[int, float]:
    return _filter_cluster_candidates_impl(candidate_scores, threshold=threshold)


def _build_cluster_tree(scoring_values: np.ndarray):
    return _build_cluster_tree_impl(scoring_values)


def build_cluster_labels_from_tree(
    *,
    cluster_tree,
    cluster_counts: Iterable[int],
) -> dict[int, np.ndarray]:
    impl = cast(
        Callable[..., dict[int, np.ndarray]], _build_cluster_labels_from_tree_impl
    )
    return impl(
        cluster_tree=cluster_tree,
        cluster_counts=cluster_counts,
    )


def fit_cluster_labels(
    scoring_values: np.ndarray,
    cluster_count: int,
    *,
    tree_engine: SignalomeTreeEngine = SIGNALOME_TREE_ENGINE_EXACT,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy = (
        SIGNALOME_CANDIDATE_SCORING_POLICY_FULL
    ),
    max_exact_tree_sites: int | None = MAX_FULL_CORRELATION_SITE_COUNT,
    cluster_tree_operations: ClusterTreeOperations | None = None,
) -> np.ndarray:
    return _fit_cluster_labels_impl(
        scoring_values=scoring_values,
        cluster_count=cluster_count,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        cluster_tree_operations=cluster_tree_operations,
    )


def cluster_median_correlation(
    site_correlations: np.ndarray,
    labels: np.ndarray,
    label: int,
) -> float:
    return _cluster_median_correlation_impl(site_correlations, labels, label)


def cluster_median_correlation_approximate(
    *,
    scoring_values: np.ndarray,
    labels: np.ndarray,
    label: int,
    max_sites_per_cluster: int,
) -> float:
    return _cluster_median_correlation_approximate_impl(
        scoring_values=scoring_values,
        labels=labels,
        label=label,
        max_sites_per_cluster=max_sites_per_cluster,
    )


def _prepare_scoring_values_for_clustering(scoring_values: np.ndarray) -> np.ndarray:
    return _prepare_scoring_values_for_clustering_impl(scoring_values)


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
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD",
    "SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY",
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
    "_ModuleSelectionComputation",
    "_ProfileDegeneracySummary",
    "_build_cluster_tree",
    "_build_module_selection_result",
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
