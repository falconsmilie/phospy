"""Internal executor for the signalome workflow."""

from __future__ import annotations

from phospy.api.results import SignalomeWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.signalomes.clustering import (
    run_signalome_clustering_engine,
)
from phospy.signalomes.context import (
    build_protein_site_context_table,
    build_site_membership_table,
)
from phospy.signalomes.science import (
    build_expanded_signalome_table,
    build_kinase_network_with_diagnostics,
    build_module_assignments,
    build_signalome_module_table,
    select_kinase_substrates,
)
from phospy.workflows.signalome.clustering_runner import SignalomeClusteringRunner
from phospy.workflows.signalome.constants import (
    SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
)
from phospy.workflows.signalome.context_tables import SignalomeContextTableBuilder
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeWorkflowRequest,
)
from phospy.workflows.signalome.module_tables import SignalomeModuleTableBuilder
from phospy.workflows.signalome.network_builder import SignalomeNetworkBuilder
from phospy.workflows.signalome.provenance import SignalomeProvenanceBuilder
from phospy.workflows.signalome.result_assembly import SignalomeResultAssembler


class SignalomeWorkflowExecutor:
    """Run signalome stage logic and assemble `SignalomeWorkflowResult`."""

    def __init__(
        self,
        *,
        clustering_runner: SignalomeClusteringRunner | None = None,
        module_table_builder: SignalomeModuleTableBuilder | None = None,
        network_builder: SignalomeNetworkBuilder | None = None,
        context_table_builder: SignalomeContextTableBuilder | None = None,
        provenance_builder: SignalomeProvenanceBuilder | None = None,
        result_assembler: SignalomeResultAssembler | None = None,
    ) -> None:
        # Dependency wiring is intentionally done here so tests can monkeypatch
        # executor-module callables and still intercept default component behavior.
        self._clustering_runner = clustering_runner or SignalomeClusteringRunner(
            run_backend_clustering=run_signalome_clustering_engine,
        )
        self._module_table_builder = (
            module_table_builder
            or SignalomeModuleTableBuilder(
                build_assignments=build_module_assignments,
                select_substrates=select_kinase_substrates,
                build_modules=build_signalome_module_table,
            )
        )
        self._network_builder = network_builder or SignalomeNetworkBuilder(
            build_network=build_kinase_network_with_diagnostics,
        )
        self._context_table_builder = (
            context_table_builder
            or SignalomeContextTableBuilder(
                build_site_membership=build_site_membership_table,
                build_protein_context=build_protein_site_context_table,
            )
        )
        self._provenance_builder = provenance_builder or SignalomeProvenanceBuilder()
        self._result_assembler = result_assembler or SignalomeResultAssembler(
            build_expanded=build_expanded_signalome_table
        )

    def run(self, request: ResolvedSignalomeWorkflowRequest) -> SignalomeWorkflowResult:
        config = request.execution_config
        execution_metadata = SignalomeClusteringRunner.collect_execution_metadata(
            request
        )
        clustering_stage = self._clustering_runner.run(
            request=request,
            config=config,
            execution_metadata=execution_metadata,
        )
        scale_guard_decision = SignalomeClusteringRunner.summarize_scale_guard(
            config=config,
            site_count=execution_metadata.downstream_score_sites,
            site_to_protein=request.site_to_protein,
            downstream_score_kinases=execution_metadata.downstream_score_kinases,
            clustering_result=clustering_stage.clustering_result,
        )
        module_stage = self._module_table_builder.run(
            request=request,
            config=config,
            clustering_result=clustering_stage.clustering_result,
            protein_modules=clustering_stage.protein_modules,
            execution_metadata=execution_metadata,
        )
        network_stage = self._network_builder.run(
            request=request,
            config=config,
            clustering_result=clustering_stage.clustering_result,
            support_summary=module_stage.support_summary,
            execution_metadata=execution_metadata,
        )
        expanded_signalome = self._result_assembler.build_expanded_signalome(
            request=request,
            config=config,
            module_assignments=module_stage.module_assignments,
            signalome_modules=module_stage.signalome_modules,
            network_edges=network_stage.edges,
            support_summary=module_stage.support_summary,
            module_count=module_stage.module_count,
            execution_metadata=execution_metadata,
        )
        context_stage = self._context_table_builder.run(
            request=request,
            config=config,
            clustering_result=clustering_stage.clustering_result,
            module_assignments=module_stage.module_assignments,
            support_summary=module_stage.support_summary,
            execution_metadata=execution_metadata,
        )
        provenance = self._provenance_builder.build(
            request=request,
            config=config,
            clustering_result=clustering_stage.clustering_result,
            module_assignments=module_stage.module_assignments,
            signalome_modules=module_stage.signalome_modules,
            network_edges=network_stage.edges,
            network_nodes=network_stage.nodes,
            candidate_correlations=network_stage.candidate_correlations,
            network_correlation_diagnostics=network_stage.correlation_diagnostics,
            expanded_signalome=expanded_signalome,
            site_membership=context_stage.site_membership,
            protein_site_context=context_stage.protein_site_context,
            scale_guard_decision=scale_guard_decision,
        )
        return self._result_assembler.assemble_result(
            request=request,
            clustering_result=clustering_stage.clustering_result,
            module_assignments=module_stage.module_assignments,
            signalome_modules=module_stage.signalome_modules,
            network_edges=network_stage.edges,
            network_nodes=network_stage.nodes,
            candidate_correlations=network_stage.candidate_correlations,
            network_correlation_diagnostics=network_stage.correlation_diagnostics,
            expanded_signalome=expanded_signalome,
            site_membership=context_stage.site_membership,
            protein_site_context=context_stage.protein_site_context,
            provenance=provenance,
        )

    @staticmethod
    def _raise_boundary_error(
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
