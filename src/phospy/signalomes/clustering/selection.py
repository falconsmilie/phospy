"""Module-count selection facade and helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.signalomes.models import SignalomeModuleSelectionDiagnostics


def _exact_module():
    # Local import avoids import cycles: exact backend uses this module too.
    from phospy.signalomes.clustering import exact_python as _exact

    return _exact


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: str = "auto",
    max_exact_cluster_tree_sites: int | None = 2000,
) -> int:
    """Select a module count from a scoring matrix."""

    return select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    ).selected_module_count


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: str = "auto",
    max_exact_cluster_tree_sites: int | None = 2000,
) -> SignalomeModuleSelectionDiagnostics:
    """Select a module count and return diagnostics."""

    exact = _exact_module()
    return exact.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_cluster_tree_sites=max_exact_cluster_tree_sites,
    )


def _compute_module_selection(**kwargs: object):
    """Compatibility forwarding hook for internal selection computation."""

    return _exact_module()._compute_module_selection(**kwargs)


def _resolve_pre_scoring_module_selection(**kwargs: object):
    """Compatibility forwarding hook for internal pre-scoring resolution."""

    return _exact_module()._resolve_pre_scoring_module_selection(**kwargs)


def _compute_candidate_cluster_scores(**kwargs: object):
    """Compatibility forwarding hook for candidate score computation."""

    return _exact_module()._compute_candidate_cluster_scores(**kwargs)


def _resolve_candidate_scoring_backend(**kwargs: object):
    """Compatibility forwarding hook for candidate scoring backend resolution."""

    return _exact_module()._resolve_candidate_scoring_backend(**kwargs)


def _select_best_candidate_count(*args: object, **kwargs: object):
    return _exact_module()._select_best_candidate_count(*args, **kwargs)


def _select_threshold_candidate(**kwargs: object):
    return _exact_module()._select_threshold_candidate(**kwargs)


def filter_cluster_candidates(*args: object, **kwargs: object):
    return _exact_module().filter_cluster_candidates(*args, **kwargs)


__all__ = [
    "filter_cluster_candidates",
    "select_module_count",
    "select_module_count_with_diagnostics",
]
