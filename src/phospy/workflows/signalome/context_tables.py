"""Signalome context table construction."""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.science.signalomes.clustering import ClusterSitesResult
from phospy.science.signalomes.context import (
    build_protein_site_context_table,
    build_site_membership_table,
)
from phospy.workflows.signalome.component_helpers import (
    prediction_shape_details,
    raise_boundary_error,
    support_details,
)
from phospy.workflows.signalome.component_models import (
    SignalomeContextTableBuildResult,
    SignalomeExecutionMetadata,
    SignalomeSupportSummary,
)
from phospy.workflows.signalome.constants import SIGNALOME_EXECUTOR_CONTEXT_TABLES_SEAM
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


class SignalomeContextTableBuilder:
    """Build site and protein context tables for signalome execution."""

    def __init__(
        self,
        *,
        build_site_membership: Callable[
            ..., pd.DataFrame
        ] = build_site_membership_table,
        build_protein_context: Callable[..., pd.DataFrame] = (
            build_protein_site_context_table
        ),
    ) -> None:
        self._build_site_membership = build_site_membership
        self._build_protein_context = build_protein_context

    def run(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        module_assignments: pd.DataFrame,
        support_summary: SignalomeSupportSummary,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeContextTableBuildResult:
        try:
            site_membership = self._build_site_membership(
                module_assignments=module_assignments,
                site_clusters=clustering_result.site_clusters,
                site_metadata=request.dataset._borrow_site_metadata_frame(),
                prediction_matrix=request.prediction_matrix,
                kinase_substrates=support_summary.kinase_substrates,
                substrate_support_cutoff=config.substrate_support_cutoff,
                assignment_policy=config.assignment_policy,
            )
            protein_site_context = self._build_protein_context(
                site_membership=site_membership
            )
            return SignalomeContextTableBuildResult(
                site_membership=site_membership,
                protein_site_context=protein_site_context,
            )
        except (WorkflowStageError, ValueError, TypeError, KeyError) as exc:
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_CONTEXT_TABLES_SEAM,
                next_action=(
                    "ensure module assignments, site clusters, site metadata, "
                    "prediction scores, and supported substrate mappings are "
                    "mutually consistent before building signalome context tables"
                ),
                **prediction_shape_details(execution_metadata),
                module_assignment_rows=int(module_assignments.shape[0]),
                module_assignment_columns=int(module_assignments.shape[1]),
                site_cluster_count=int(clustering_result.site_clusters.nunique()),
                **support_details(support_summary.support_counts),
                substrate_support_cutoff=config.substrate_support_cutoff,
                assignment_policy=config.assignment_policy,
                stage_error=str(exc),
            )


__all__ = ["SignalomeContextTableBuilder"]
