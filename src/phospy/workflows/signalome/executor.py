"""Internal executor for the signalome workflow."""

from __future__ import annotations

from phospy.api.results import SignalomeWorkflowResult
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from phospy.workflows.signalome.contracts import ResolvedSignalomeWorkflowRequest
from phospy.workflows.signalome.science import (
    build_kinase_network,
    build_module_assignments,
    build_signalome_module_table,
    select_kinase_substrates,
)


class SignalomeWorkflowExecutor:
    """Run signalome stage logic and assemble `SignalomeWorkflowResult`."""

    def run(self, request: ResolvedSignalomeWorkflowRequest) -> SignalomeWorkflowResult:
        module_assignments = build_module_assignments(
            prediction_matrix=request.prediction_matrix,
            site_to_protein=request.site_to_protein,
        )
        kinase_substrates = select_kinase_substrates(
            prediction_matrix=request.prediction_matrix,
            cutoff=float(request.config.signalome_cutoff),
        )
        signalome_modules = build_signalome_module_table(
            module_assignments=module_assignments,
            kinase_substrates=kinase_substrates,
            kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
        )
        network_edges, network_nodes = build_kinase_network(
            score_matrix=request.score_matrix,
            kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
            kinase_substrates=kinase_substrates,
            threshold=float(request.config.signalome_cutoff),
        )
        return SignalomeWorkflowResult(
            dataset=request.dataset,
            kinase_result=request.kinase_result,
            module_assignments=SignalomeAssignments(table=module_assignments),
            signalome_modules=SignalomeModules(table=signalome_modules),
            kinase_network=KinaseNetwork(
                edges=network_edges,
                nodes=network_nodes,
            ),
            expanded_signalome=None,
        )
