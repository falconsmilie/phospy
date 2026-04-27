"""Shared helper functions for signalome workflow components."""

from __future__ import annotations

import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.signalomes.clustering import ClusterSitesResult
from phospy.workflows.signalome.component_models import SignalomeExecutionMetadata
from phospy.workflows.signalome.constants import (
    SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
)


def raise_boundary_error(
    *,
    seam: str,
    next_action: str,
    **details: object,
) -> None:
    raise WorkflowBoundaryError(
        seam=seam,
        next_action=next_action,
        details=details,
        message_prefix=SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
    )


def summarize_support(kinase_substrates: dict[str, tuple[str, ...]]) -> dict[str, int]:
    supported_site_ids: set[str] = set()
    supported_kinases = 0
    for substrates in kinase_substrates.values():
        resolved_sites = tuple(str(site_id) for site_id in substrates)
        if not resolved_sites:
            continue
        supported_kinases += 1
        supported_site_ids.update(resolved_sites)
    return {
        "supported_sites": int(len(supported_site_ids)),
        "supported_kinases": int(supported_kinases),
    }


def score_variance_kinases(downstream_score_matrix: pd.DataFrame) -> int:
    if downstream_score_matrix.empty:
        return 0
    variances = downstream_score_matrix.astype(float).var(axis=0, ddof=0)
    return int((variances > 0.0).sum())


def prediction_shape_details(
    execution_metadata: SignalomeExecutionMetadata,
) -> dict[str, int]:
    return {
        "prediction_sites": execution_metadata.prediction_sites,
        "prediction_kinases": execution_metadata.prediction_kinases,
    }


def support_details(support_counts: dict[str, int]) -> dict[str, int]:
    return {
        "supported_sites": int(support_counts["supported_sites"]),
        "supported_kinases": int(support_counts["supported_kinases"]),
    }


def module_selection_details(
    clustering_result: ClusterSitesResult,
) -> dict[str, int | str]:
    diagnostics = clustering_result.module_selection_diagnostics
    return {
        "selected_module_count": int(diagnostics.selected_module_count),
        "requested_module_count": requested_module_count_label(
            diagnostics.requested_module_count
        ),
    }


def requested_module_count_label(requested_module_count: int | None) -> int | str:
    if requested_module_count is None:
        return "auto"
    return int(requested_module_count)


__all__ = [
    "module_selection_details",
    "prediction_shape_details",
    "raise_boundary_error",
    "requested_module_count_label",
    "score_variance_kinases",
    "summarize_support",
    "support_details",
]
