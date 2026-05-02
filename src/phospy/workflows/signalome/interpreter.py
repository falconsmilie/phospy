"""Internal interpreter for signalome workflow requests."""

from __future__ import annotations

from phospy.api.requests import SignalomeWorkflowRequest
from phospy.prediction.scoring import select_downstream_score_matrix
from phospy.scientific_policies import resolve_score_preconditioning_policy
from phospy.signalomes.clustering.policies import (
    resolve_candidate_scoring_policy_definition,
)
from phospy.workflows.signalome.alignment_diagnostics import (
    SignalomeAlignmentDiagnosticsBuilder,
)
from phospy.workflows.signalome.boundary_errors import (
    raise_signalome_boundary_error,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)
from phospy.workflows.signalome.matrix_alignment import SignalomeMatrixAligner
from phospy.workflows.signalome.protein_resolution import SignalomeProteinResolver
from phospy.workflows.signalome.score_matrix_selection import (
    SignalomeScoreMatrixSelector,
)
from phospy.workflows.signalome.score_preconditioning import (
    SignalomeScorePreconditioner,
)


class SignalomeWorkflowInterpreter:
    """Resolve a validated signalome request into execution-ready inputs."""

    def __init__(
        self,
        *,
        score_matrix_selector: SignalomeScoreMatrixSelector | None = None,
        matrix_aligner: SignalomeMatrixAligner | None = None,
        score_preconditioner: SignalomeScorePreconditioner | None = None,
        protein_resolver: SignalomeProteinResolver | None = None,
        alignment_diagnostics_builder: (
            SignalomeAlignmentDiagnosticsBuilder | None
        ) = None,
    ) -> None:
        self._score_matrix_selector = (
            score_matrix_selector
            or SignalomeScoreMatrixSelector(
                select_matrix=select_downstream_score_matrix
            )
        )
        self._matrix_aligner = matrix_aligner or SignalomeMatrixAligner()
        self._score_preconditioner = (
            score_preconditioner or SignalomeScorePreconditioner()
        )
        self._protein_resolver = protein_resolver or SignalomeProteinResolver()
        self._alignment_diagnostics_builder = (
            alignment_diagnostics_builder or SignalomeAlignmentDiagnosticsBuilder()
        )

    def run(
        self, request: SignalomeWorkflowRequest
    ) -> ResolvedSignalomeWorkflowRequest:
        score_selection = self._score_matrix_selector.run(
            request.kinase_result.scoring_result
        )
        aligned_matrices = self._matrix_aligner.run(
            dataset_sites=request.kinase_result.dataset.phospho.index,
            prediction_matrix=request.kinase_result.prediction_result.pred_mat,
            downstream_score_matrix=score_selection.downstream_score_matrix,
            downstream_score_source=score_selection.downstream_score_source,
        )
        execution_config = self._resolve_execution_config(request)
        preconditioning_result = self._score_preconditioner.run(
            score_matrix=aligned_matrices.aligned_downstream_score_matrix,
            policy=execution_config.score_preconditioning_policy,
        )
        if preconditioning_result.downstream_score_matrix.empty:
            raise_signalome_boundary_error(
                seam=SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
                next_action=(
                    "ensure upstream downstream scores retain at least one "
                    "interpretable site after score preconditioning"
                ),
                aligned_score_sites=int(
                    aligned_matrices.aligned_downstream_score_matrix.shape[0]
                ),
                aligned_score_kinases=int(
                    aligned_matrices.aligned_downstream_score_matrix.shape[1]
                ),
                dropped_all_missing_row_count=(
                    preconditioning_result.diagnostics.dropped_all_missing_row_count
                ),
                retained_row_count=preconditioning_result.diagnostics.retained_row_count,
                score_preconditioning_policy=execution_config.score_preconditioning_policy,
            )
        retained_site_index = preconditioning_result.downstream_score_matrix.index
        retained_prediction_matrix = aligned_matrices.aligned_prediction_matrix.loc[
            retained_site_index
        ]
        site_to_protein = self._protein_resolver.run(
            dataset=request.kinase_result.dataset,
            site_index=retained_site_index,
            removed_by_score_preconditioning_count=int(
                aligned_matrices.aligned_site_index.size - retained_site_index.size
            ),
        )
        alignment_diagnostics = self._alignment_diagnostics_builder.run(
            dataset_sites=aligned_matrices.dataset_site_index,
            prediction_sites=aligned_matrices.resolved_prediction_matrix.index,
            score_sites=aligned_matrices.resolved_downstream_score_matrix.index,
            shared_sites=aligned_matrices.aligned_site_index,
            retained_sites=retained_site_index,
            prediction_kinases=aligned_matrices.resolved_prediction_matrix.columns,
            score_kinases=aligned_matrices.resolved_downstream_score_matrix.columns,
            shared_kinases=aligned_matrices.aligned_kinase_index,
            interpreted_protein_sites=aligned_matrices.aligned_site_index,
            retained_protein_sites=retained_site_index,
        )
        return ResolvedSignalomeWorkflowRequest(
            dataset=request.kinase_result.dataset,
            kinase_result=request.kinase_result,
            execution_config=execution_config,
            downstream_score_matrix=preconditioning_result.downstream_score_matrix,
            downstream_score_source=score_selection.downstream_score_source,
            prediction_matrix=retained_prediction_matrix,
            site_to_protein=site_to_protein,
            score_preconditioning_diagnostics=preconditioning_result.diagnostics,
            alignment_diagnostics=alignment_diagnostics,
            downstream_score_selection_policy=(
                score_selection.downstream_score_selection_policy
            ),
        )

    @staticmethod
    def _resolve_execution_config(
        request: SignalomeWorkflowRequest,
    ) -> ResolvedSignalomeExecutionConfig:
        return ResolvedSignalomeExecutionConfig(
            substrate_support_cutoff=float(
                request.config.scientific.substrate_support_cutoff
            ),
            network_correlation_threshold=float(
                request.config.output.network_correlation_threshold
            ),
            network_policy=request.config.output.network_policy,
            assignment_policy=request.config.scientific.assignment_policy,
            score_preconditioning_policy=(
                request.config.validation.score_preconditioning_policy
            ),
            module_selection_primary_threshold=float(
                request.config.clustering.module_selection_primary_correlation_threshold
            ),
            module_selection_fallback_threshold=float(
                request.config.clustering.module_selection_fallback_correlation_threshold
            ),
            module_selection_max_clusters=int(
                request.config.clustering.module_selection_max_clusters
            ),
            candidate_scoring_policy=request.config.clustering.candidate_scoring_policy,
            max_exact_tree_sites=int(request.config.performance.max_exact_tree_sites),
            max_full_candidate_scoring_sites=int(
                request.config.performance.max_full_candidate_scoring_sites
            ),
            requested_module_count=(
                None
                if request.config.clustering.module_count is None
                else int(request.config.clustering.module_count)
            ),
            clustering_engine=str(request.config.clustering.clustering_engine),
            candidate_scoring_policy_definition=(
                resolve_candidate_scoring_policy_definition(
                    candidate_scoring_policy=(
                        request.config.clustering.candidate_scoring_policy
                    )
                )
            ),
            score_preconditioning_policy_definition=(
                resolve_score_preconditioning_policy(
                    policy=str(request.config.validation.score_preconditioning_policy)
                )
            ),
        )


__all__ = ["SignalomeWorkflowInterpreter"]
