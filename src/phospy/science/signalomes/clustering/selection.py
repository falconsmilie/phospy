"""Module-count selection facade and helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.science.signalomes.clustering.candidate_selection import (
    filter_cluster_candidates,
)
from phospy.science.signalomes.clustering.orchestration import (
    select_module_count as _select_module_count,
)
from phospy.science.signalomes.clustering.orchestration import (
    select_module_count_with_diagnostics as _select_module_count_with_diagnostics,
)
from phospy.science.signalomes.clustering.policies import (
    SignalomeClusteringScoringMode,
)
from phospy.science.signalomes.models import (
    SignalomeModuleSelectionDiagnostics,
)


def select_module_count(
    scoring_values: pd.DataFrame | np.ndarray,
    *,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_tree_sites: int | None = 2000,
) -> int:
    """Select a module count from a scoring matrix."""

    return _select_module_count(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
    )


def select_module_count_with_diagnostics(
    *,
    scoring_values: pd.DataFrame | np.ndarray,
    requested_module_count: int | None = None,
    primary_threshold: float = 0.5,
    fallback_threshold: float = 0.1,
    max_clusters: int = 10,
    scoring_mode: SignalomeClusteringScoringMode = "auto",
    max_exact_tree_sites: int | None = 2000,
) -> SignalomeModuleSelectionDiagnostics:
    """Select a module count and return diagnostics."""

    return _select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=requested_module_count,
        primary_threshold=primary_threshold,
        fallback_threshold=fallback_threshold,
        max_clusters=max_clusters,
        scoring_mode=scoring_mode,
        max_exact_tree_sites=max_exact_tree_sites,
    )


__all__ = [
    "filter_cluster_candidates",
    "select_module_count",
    "select_module_count_with_diagnostics",
]
