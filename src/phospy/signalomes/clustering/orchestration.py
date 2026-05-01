"""Signalome clustering orchestration coordination layer."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.signalomes.clustering.candidate_scoring import (
    _CandidateClusterScoreResult,
    _ProfileDegeneracySummary,
    build_correlation_exclusion_note,
    build_correlation_matrix_with_exclusions,
    cluster_median_correlation,
    cluster_median_correlation_approximate,
    compute_candidate_cluster_scores,
    resolve_candidate_scoring_policy,
    summarize_profile_degeneracy,
)
from phospy.signalomes.clustering.candidate_selection import (
    _ModuleSelectionComputation,
    build_module_selection_result,
    filter_cluster_candidates,
    resolve_pre_scoring_module_selection,
    select_best_candidate_count,
    select_threshold_candidate,
)
from phospy.signalomes.clustering.contracts import ClusterTreeEngine
from phospy.signalomes.clustering.diagnostic_schemas import (
    SignalomeBackendDiagnostics,
    SignalomeCandidateScoringSamplingDiagnostics,
    SignalomeClusteringLimitMetadata,
    SignalomeClusteringThresholdMetadata,
    SignalomeTreeEngineDiagnostics,
    build_backend_diagnostics,
    validate_backend_diagnostics,
    validate_candidate_scoring_sampling_diagnostics,
)
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
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS,
    SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES,
    SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE,
    SIGNALOME_CLUSTERING_SCORING_MODE_AUTO,
    SIGNALOME_CLUSTERING_SCORING_MODE_EXACT,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE,
    SIGNALOME_TREE_ENGINE_EXACT,
    SignalomeCandidateScoringPolicy,
    SignalomeClusteringMissingValuePolicy,
    SignalomeClusteringScoringMode,
    SignalomeTreeEngine,
    _CandidateScoringMode,
)
from phospy.signalomes.clustering.protein_modules import derive_protein_modules
from phospy.signalomes.clustering.scale_guards import resolve_max_exact_tree_sites
from phospy.signalomes.clustering.tree_building import (
    ClusterTreeOperations,
    ClusterTreeOperationsAdapter,
    build_cluster_labels_from_tree,
    build_cluster_tree,
    fit_cluster_labels,
    prepare_scoring_values_for_clustering,
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
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None = (
        None
    )
    candidate_scoring_evaluated: bool = False
    candidate_scoring_skip_reason: str | None = None
    backend_name: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON
    backend_version: str = SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON_VERSION
    approximation_used: bool = False
    backend_diagnostics: SignalomeBackendDiagnostics | None = None

    def __post_init__(self) -> None:
        if self.candidate_scoring_sampling is not None:
            validate_candidate_scoring_sampling_diagnostics(
                self.candidate_scoring_sampling,
                field_name="cluster_sites_result.candidate_scoring_sampling",
            )
        if self.backend_diagnostics is not None:
            validate_backend_diagnostics(
                self.backend_diagnostics,
                field_name="cluster_sites_result.backend_diagnostics",
            )


@dataclass(frozen=True, slots=True)
class ScorePreconditioner:
    """Prepare matrix values and profile diagnostics for clustering/scoring."""

    def for_clustering(self, scoring_values: np.ndarray) -> np.ndarray:
        return prepare_scoring_values_for_clustering(scoring_values)

    def profile_degeneracy(
        self, scoring_values: np.ndarray
    ) -> _ProfileDegeneracySummary:
        return summarize_profile_degeneracy(scoring_values)

    def exclusion_note(self, summary: _ProfileDegeneracySummary) -> str:
        return build_correlation_exclusion_note(summary)


@dataclass(frozen=True, slots=True)
class ModuleScorer:
    """Resolve candidate-scoring policy and score candidate module counts."""

    def resolve_policy(
        self,
        *,
        scoring_mode: SignalomeClusteringScoringMode,
        candidate_scoring_policy: SignalomeCandidateScoringPolicy | None,
        n_sites: int,
        max_full_candidate_scoring_sites: int,
    ) -> SignalomeCandidateScoringPolicy:
        return resolve_candidate_scoring_policy(
            scoring_mode=scoring_mode,
            candidate_scoring_policy=candidate_scoring_policy,
            n_sites=n_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        )

    def score_candidates(
        self,
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
        cluster_tree_operations: ClusterTreeOperations | None,
    ) -> _CandidateClusterScoreResult:
        return compute_candidate_cluster_scores(
            clustering_values=clustering_values,
            correlation_values=correlation_values,
            candidate_range=candidate_range,
            profile_degeneracy=profile_degeneracy,
            n_sites=n_sites,
            scoring_mode=scoring_mode,
            tree_engine=tree_engine,
            candidate_scoring_policy=candidate_scoring_policy,
            max_exact_tree_sites=max_exact_tree_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
            cluster_tree_operations=cluster_tree_operations,
        )


@dataclass(frozen=True, slots=True)
class ModuleSelector:
    """Select the signalome module count and return diagnostics payloads."""

    preconditioner: ScorePreconditioner
    scorer: ModuleScorer

    def select(
        self,
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

        resolved_max_exact_tree_sites = resolve_max_exact_tree_sites(
            max_exact_tree_sites
        )
        scoring_array = np.asarray(scoring_values, dtype=float)
        if scoring_array.ndim != 2:
            raise ValueError("scoring_values must be a 2D array")

        n_sites = int(scoring_array.shape[0])
        requested_module_count = validate_requested_module_count(
            requested_module_count=requested_module_count,
            available_clustering_site_count=n_sites,
            field_name="signalome workflow request config.clustering.module_count",
        )
        profile_degeneracy = self.preconditioner.profile_degeneracy(scoring_array)
        correlation_exclusion_note = self.preconditioner.exclusion_note(
            profile_degeneracy
        )

        early_selection, resolved_max_clusters = resolve_pre_scoring_module_selection(
            requested_module_count=requested_module_count,
            n_sites=n_sites,
            max_clusters=max_clusters,
            profile_degeneracy=profile_degeneracy,
            correlation_exclusion_note=correlation_exclusion_note,
        )
        if early_selection is not None:
            return early_selection

        resolved_candidate_scoring_policy = self.scorer.resolve_policy(
            scoring_mode=scoring_mode,
            candidate_scoring_policy=candidate_scoring_policy,
            n_sites=n_sites,
            max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        )
        clustering_values = self.preconditioner.for_clustering(scoring_array)
        candidate_score_result = self.scorer.score_candidates(
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

        primary_selection = select_threshold_candidate(
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

        fallback_selection = select_threshold_candidate(
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

        return build_module_selection_result(
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


@dataclass(frozen=True, slots=True)
class ClusterTreeBuilder:
    """Build final cluster labels for the selected module count."""

    preconditioner: ScorePreconditioner

    def final_labels(
        self,
        *,
        scoring_values: np.ndarray,
        module_count: int,
        cached_labels: np.ndarray | None,
        tree_engine: SignalomeTreeEngine,
        candidate_scoring_policy: SignalomeCandidateScoringPolicy | None,
        max_exact_tree_sites: int | None,
        cluster_tree_operations: ClusterTreeOperations | None,
    ) -> tuple[np.ndarray, bool]:
        n_sites = int(scoring_values.shape[0])
        if module_count == 1:
            return np.ones(n_sites, dtype=int), False
        if cached_labels is not None:
            return cached_labels.astype(int, copy=False) + 1, True

        labels = fit_cluster_labels(
            scoring_values=self.preconditioner.for_clustering(scoring_values),
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
        return labels + 1, True


@dataclass(frozen=True, slots=True)
class ProteinMapper:
    """Map site-level cluster assignments to protein modules."""

    def map(
        self,
        *,
        site_clusters: pd.Series,
        site_to_protein: pd.Series,
    ) -> pd.Series:
        return derive_protein_modules(
            site_clusters=site_clusters,
            site_to_protein=site_to_protein,
        )


@dataclass(frozen=True, slots=True)
class SignalomeDiagnosticsBuilder:
    """Build backend diagnostics and typed metadata for engine results."""

    def backend_diagnostics(
        self,
        *,
        clustering_engine: str,
        tree_engine: ClusterTreeEngine,
        tree_engine_diagnostics: SignalomeTreeEngineDiagnostics,
        selected_module_count: int,
        input_site_count: int,
        exact_tree_path_used: bool,
    ) -> SignalomeBackendDiagnostics:
        return build_backend_diagnostics(
            backend_name=str(clustering_engine),
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


_SCORE_PRECONDITIONER = ScorePreconditioner()
_MODULE_SCORER = ModuleScorer()
_MODULE_SELECTOR = ModuleSelector(
    preconditioner=_SCORE_PRECONDITIONER,
    scorer=_MODULE_SCORER,
)
_TREE_BUILDER = ClusterTreeBuilder(preconditioner=_SCORE_PRECONDITIONER)
_PROTEIN_MAPPER = ProteinMapper()
_DIAGNOSTICS_BUILDER = SignalomeDiagnosticsBuilder()


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
    selection = _MODULE_SELECTOR.select(
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

    module_count = validate_cluster_count_for_site_count(
        cluster_count=int(selection.diagnostics.selected_module_count),
        available_clustering_site_count=int(scoring_values.shape[0]),
        field_name="selected_module_count",
    )
    labels, built_for_final_assignment = _TREE_BUILDER.final_labels(
        scoring_values=scoring_values,
        module_count=module_count,
        cached_labels=selection.candidate_labels.get(module_count),
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        cluster_tree_operations=cluster_tree_operations,
    )
    exact_cluster_tree_built = bool(
        built_for_final_assignment or selection.exact_cluster_tree_built
    )

    return ClusterSitesResult(
        site_clusters=pd.Series(
            labels,
            index=scoring_matrix.index.copy(),
            dtype=int,
            name=SITE_CLUSTER_COLUMN,
        ),
        module_selection_diagnostics=selection.diagnostics,
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
    return _MODULE_SELECTOR.select(
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
    return _MODULE_SELECTOR.select(
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


def _resolve_pre_scoring_module_selection(
    *,
    requested_module_count: int | None,
    n_sites: int,
    max_clusters: int,
    profile_degeneracy: _ProfileDegeneracySummary,
    correlation_exclusion_note: str,
) -> tuple[_ModuleSelectionComputation | None, int]:
    return resolve_pre_scoring_module_selection(
        requested_module_count=requested_module_count,
        n_sites=n_sites,
        max_clusters=max_clusters,
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
    )


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
    return compute_candidate_cluster_scores(
        clustering_values=clustering_values,
        correlation_values=correlation_values,
        candidate_range=candidate_range,
        profile_degeneracy=profile_degeneracy,
        n_sites=n_sites,
        scoring_mode=scoring_mode,
        tree_engine=tree_engine,
        candidate_scoring_policy=candidate_scoring_policy,
        max_exact_tree_sites=max_exact_tree_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
        cluster_tree_operations=cluster_tree_operations,
    )


def _resolve_candidate_scoring_policy(
    *,
    scoring_mode: SignalomeClusteringScoringMode,
    candidate_scoring_policy: SignalomeCandidateScoringPolicy | None,
    n_sites: int,
    max_full_candidate_scoring_sites: int,
) -> SignalomeCandidateScoringPolicy:
    return resolve_candidate_scoring_policy(
        scoring_mode=scoring_mode,
        candidate_scoring_policy=candidate_scoring_policy,
        n_sites=n_sites,
        max_full_candidate_scoring_sites=max_full_candidate_scoring_sites,
    )


def _select_best_candidate_count(candidate_scores: dict[int, float]) -> int:
    return select_best_candidate_count(candidate_scores)


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
    candidate_scoring_sampling: SignalomeCandidateScoringSamplingDiagnostics | None,
) -> _ModuleSelectionComputation | None:
    return select_threshold_candidate(
        candidate_scores=candidate_scores,
        candidate_labels=candidate_labels,
        max_clusters=max_clusters,
        threshold=threshold,
        requested_module_count=requested_module_count,
        reason=reason,
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=correlation_exclusion_note,
        approximation_note=approximation_note,
        tree_engine=tree_engine,
        candidate_scoring_mode=candidate_scoring_mode,
        exact_cluster_tree_built=exact_cluster_tree_built,
        candidate_scoring_evaluated=candidate_scoring_evaluated,
        candidate_scoring_skip_reason=candidate_scoring_skip_reason,
        candidate_scoring_sampling=candidate_scoring_sampling,
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
    backend_diagnostics: SignalomeTreeEngineDiagnostics,
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
    protein_modules = _PROTEIN_MAPPER.map(
        site_clusters=clustering_result.site_clusters,
        site_to_protein=request.site_to_protein,
    )
    selected_module_count = int(
        clustering_result.module_selection_diagnostics.selected_module_count
    )
    resolved_backend_diagnostics = _DIAGNOSTICS_BUILDER.backend_diagnostics(
        clustering_engine=clustering_engine,
        tree_engine=tree_engine,
        tree_engine_diagnostics=backend_diagnostics,
        selected_module_count=selected_module_count,
        input_site_count=int(request.scoring_matrix.shape[0]),
        exact_tree_path_used=bool(clustering_result.exact_cluster_tree_built),
    )
    threshold_metadata = _DIAGNOSTICS_BUILDER.threshold_metadata(
        primary_threshold=request.primary_threshold,
        fallback_threshold=request.fallback_threshold,
    )
    limit_metadata = _DIAGNOSTICS_BUILDER.limit_metadata(
        max_exact_tree_sites=request.max_exact_tree_sites,
        max_full_candidate_scoring_sites=request.max_full_candidate_scoring_sites,
        max_clusters=request.max_clusters,
    )
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
        threshold_metadata=threshold_metadata,
        limit_metadata=limit_metadata,
    )


def _build_cluster_tree(scoring_values: np.ndarray):
    return build_cluster_tree(scoring_values)


def _prepare_scoring_values_for_clustering(scoring_values: np.ndarray) -> np.ndarray:
    return prepare_scoring_values_for_clustering(scoring_values)


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
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_APPLIES_TO",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_COLUMN_MEDIAN_IMPUTATION_WITH_ZERO_FOR_ALL_MISSING_COLUMNS",
    "SIGNALOME_CLUSTERING_MISSING_VALUE_POLICY_IMPUTED_VALUES_EXPOSED_IN_OUTPUT_TABLES",
    "SIGNALOME_CLUSTERING_SCORING_MODE_APPROXIMATE",
    "SIGNALOME_CLUSTERING_SCORING_MODE_AUTO",
    "SIGNALOME_CLUSTERING_SCORING_MODE_EXACT",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE",
    "SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_SINGLE_MODULE",
    "SignalomeCandidateScoringPolicy",
    "SignalomeClusteringMissingValuePolicy",
    "SignalomeTreeEngine",
    "SignalomeClusteringScoringMode",
    "_CandidateClusterScoreResult",
    "_CandidateScoringMode",
    "_ModuleSelectionComputation",
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
