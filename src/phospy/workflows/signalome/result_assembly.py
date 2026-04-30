"""Expanded-result construction and public result assembly for signalome workflow."""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

from phospy.api.results import SignalomeWorkflowResult
from phospy.errors.workflows import WorkflowStageError
from phospy.provenance.models import RunProvenance
from phospy.signalomes.clustering import ClusterSitesResult
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
    SignalomeNetworkCorrelationDiagnostics,
)
from phospy.signalomes.science import build_expanded_signalome_table
from phospy.workflows.signalome.component_helpers import (
    prediction_shape_details,
    raise_boundary_error,
    support_details,
)
from phospy.workflows.signalome.component_models import (
    SignalomeExecutionMetadata,
    SignalomeSupportSummary,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_EXECUTOR_EXPANDED_SIGNALOME_SEAM,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


class SignalomeResultAssembler:
    """Build expanded signalome output and assemble the public workflow result."""

    _EXPANDED_SIGNALOME_SEAM = SIGNALOME_EXECUTOR_EXPANDED_SIGNALOME_SEAM

    def __init__(
        self,
        *,
        build_expanded: Callable[..., pd.DataFrame] = build_expanded_signalome_table,
    ) -> None:
        self._build_expanded = build_expanded

    def build_expanded_signalome(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        module_assignments: pd.DataFrame,
        signalome_modules: pd.DataFrame,
        network_edges: pd.DataFrame,
        support_summary: SignalomeSupportSummary,
        module_count: int,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> pd.DataFrame:
        try:
            return self._build_expanded(
                module_assignments=module_assignments,
                signalome_modules=signalome_modules,
                kinase_network_edges=network_edges,
                kinase_substrates=support_summary.kinase_substrates,
                assignment_policy=config.assignment_policy,
            )
        except WorkflowStageError as exc:
            raise_boundary_error(
                seam=self._EXPANDED_SIGNALOME_SEAM,
                next_action=(
                    "ensure module assignments, module signal, and network topology "
                    "are mutually consistent for expanded signalome output"
                ),
                assignment_policy=config.assignment_policy,
                module_count=module_count,
                **support_details(support_summary.support_counts),
                **prediction_shape_details(execution_metadata),
                stage_error=str(exc),
            )

    @staticmethod
    def assemble_result(
        *,
        request: ResolvedSignalomeWorkflowRequest,
        clustering_result: ClusterSitesResult,
        module_assignments: pd.DataFrame,
        signalome_modules: pd.DataFrame,
        network_edges: pd.DataFrame,
        network_nodes: pd.DataFrame,
        candidate_correlations: pd.DataFrame,
        network_correlation_diagnostics: SignalomeNetworkCorrelationDiagnostics,
        expanded_signalome: pd.DataFrame,
        site_membership: pd.DataFrame,
        protein_site_context: pd.DataFrame,
        provenance: RunProvenance | None,
    ) -> SignalomeWorkflowResult:
        return SignalomeWorkflowResult._from_owned(
            dataset=request.dataset,
            kinase_result=request.kinase_result,
            module_assignments=SignalomeAssignments._from_owned(
                table=module_assignments
            ),
            signalome_modules=SignalomeModules._from_owned(table=signalome_modules),
            kinase_network=KinaseNetwork._from_owned(
                edges=network_edges,
                nodes=network_nodes,
                candidate_correlations=candidate_correlations,
                correlation_diagnostics=network_correlation_diagnostics,
            ),
            module_selection_diagnostics=clustering_result.module_selection_diagnostics,
            score_preconditioning_diagnostics=request.score_preconditioning_diagnostics,
            alignment_diagnostics=request.alignment_diagnostics,
            expanded_signalome=expanded_signalome,
            site_membership=site_membership,
            protein_site_context=protein_site_context,
            provenance=provenance,
        )


__all__ = ["SignalomeResultAssembler"]
