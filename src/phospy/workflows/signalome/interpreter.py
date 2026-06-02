"""Internal interpreter for signalome workflow requests."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.requests import SignalomeWorkflowRequest
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.prediction.scoring import select_downstream_score_matrix
from phospy.science.signalomes.clustering.policies import (
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
from phospy.workflows.signalome.scientific_policies import (
    resolve_score_preconditioning_policy,
)
from phospy.workflows.signalome.score_matrix_selection import (
    SignalomeScoreMatrixSelector,
)
from phospy.workflows.signalome.score_preconditioning import (
    SignalomeScorePreconditioner,
)

_SITE_KEY_COLUMN = "site_key"


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
        dataset = request.kinase_result.dataset
        dataset_phospho = dataset._borrow_phospho_frame()
        dataset_site_metadata = dataset._borrow_site_metadata_frame()
        dataset_sites = pd.Index(
            dataset_phospho.index.astype(str),
            name=dataset_phospho.index.name,
        )
        prediction_result = request.kinase_result.prediction_result
        prediction_matrix_input = prediction_result._borrow_pred_mat_frame()
        score_selection = self._score_matrix_selector.run(
            request.kinase_result.scoring_result
        )
        aligned_matrices = self._matrix_aligner.run(
            dataset_sites=dataset_sites,
            prediction_matrix=prediction_matrix_input,
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
        (
            aligned_site_key_index,
            retained_site_key_index,
        ) = _resolve_site_key_indexes(
            site_metadata=dataset_site_metadata,
            aligned_site_index=aligned_matrices.aligned_site_index,
            retained_site_index=retained_site_index,
        )
        prediction_output_index = pd.Index(
            retained_site_key_index.tolist(),
            name=(
                prediction_matrix_input.index.name
                if prediction_matrix_input.index.name is not None
                else retained_site_key_index.name
            ),
        )
        downstream_output_index = pd.Index(
            retained_site_key_index.tolist(),
            name=(
                score_selection.downstream_score_matrix.index.name
                if score_selection.downstream_score_matrix.index.name is not None
                else retained_site_key_index.name
            ),
        )
        retained_prediction_matrix = retained_prediction_matrix.copy(deep=False)
        retained_prediction_matrix.index = prediction_output_index
        downstream_score_matrix = preconditioning_result.downstream_score_matrix.copy(
            deep=False
        )
        downstream_score_matrix.index = downstream_output_index
        site_to_protein = self._protein_resolver.run(
            dataset=dataset,
            site_index=retained_site_index,
            removed_by_score_preconditioning_count=int(
                aligned_matrices.aligned_site_index.size - retained_site_index.size
            ),
        )
        site_to_protein.index = retained_site_key_index.copy()
        dataset_sites_for_diagnostics = _map_site_index_to_site_keys(
            site_metadata=dataset_site_metadata,
            site_index=aligned_matrices.dataset_site_index,
        )
        prediction_sites_for_diagnostics = _map_site_index_to_site_keys(
            site_metadata=dataset_site_metadata,
            site_index=aligned_matrices.resolved_prediction_matrix.index,
            allow_unmapped=True,
        )
        score_sites_for_diagnostics = _map_site_index_to_site_keys(
            site_metadata=dataset_site_metadata,
            site_index=aligned_matrices.resolved_downstream_score_matrix.index,
            allow_unmapped=True,
        )
        shared_sites_for_diagnostics = _map_site_index_to_site_keys(
            site_metadata=dataset_site_metadata,
            site_index=aligned_matrices.aligned_site_index,
        )
        retained_sites_for_diagnostics = _map_site_index_to_site_keys(
            site_metadata=dataset_site_metadata,
            site_index=retained_site_index,
        )
        alignment_diagnostics = self._alignment_diagnostics_builder.run(
            dataset_sites=dataset_sites_for_diagnostics,
            prediction_sites=prediction_sites_for_diagnostics,
            score_sites=score_sites_for_diagnostics,
            shared_sites=shared_sites_for_diagnostics,
            retained_sites=retained_sites_for_diagnostics,
            prediction_kinases=aligned_matrices.resolved_prediction_matrix.columns,
            score_kinases=aligned_matrices.resolved_downstream_score_matrix.columns,
            shared_kinases=aligned_matrices.aligned_kinase_index,
            interpreted_protein_sites=aligned_site_key_index,
            retained_protein_sites=retained_site_key_index,
        )
        return ResolvedSignalomeWorkflowRequest(
            dataset=dataset,
            kinase_result=request.kinase_result,
            execution_config=execution_config,
            downstream_score_matrix=downstream_score_matrix,
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


def _resolve_site_key_indexes(
    *,
    site_metadata: pd.DataFrame,
    aligned_site_index: pd.Index,
    retained_site_index: pd.Index,
) -> tuple[pd.Index, pd.Index]:
    if _SITE_KEY_COLUMN not in site_metadata.columns:
        raise WorkflowBoundaryError(
            "signalome interpreter requires dataset.site_metadata.site_key"
        )
    aligned_metadata = site_metadata.reindex(aligned_site_index)
    aligned_site_keys = (
        aligned_metadata.loc[:, _SITE_KEY_COLUMN].fillna("").astype(str).str.strip()
    )
    if (aligned_site_keys == "").any():
        raise WorkflowBoundaryError(
            "signalome interpreter cannot resolve site_key for aligned sites from "
            "dataset.site_metadata.site_key"
        )
    if aligned_site_keys.duplicated().any():
        raise WorkflowBoundaryError(
            "signalome interpreter requires unique site_key values across aligned sites"
        )
    if aligned_site_keys.tolist() != aligned_site_index.astype(str).tolist():
        raise WorkflowBoundaryError(
            "signalome interpreter requires dataset.site_metadata.site_key to "
            "match aligned site indexes"
        )
    aligned_site_key_index = pd.Index(
        aligned_site_keys.tolist(),
        name=_resolved_index_name_for_site_keys(
            resolved_site_keys=aligned_site_keys.tolist(),
            source_site_index=aligned_site_index,
        ),
    )
    retained_mask = aligned_site_index.isin(retained_site_index)
    retained_site_key_index = pd.Index(
        aligned_site_keys.loc[retained_mask].tolist(),
        name=_resolved_index_name_for_site_keys(
            resolved_site_keys=aligned_site_keys.loc[retained_mask].tolist(),
            source_site_index=retained_site_index,
        ),
    )
    return aligned_site_key_index, retained_site_key_index


def _map_site_index_to_site_keys(
    *,
    site_metadata: pd.DataFrame,
    site_index: pd.Index,
    allow_unmapped: bool = False,
) -> pd.Index:
    source_index = pd.Index(site_index.astype(str).tolist(), name=site_index.name)
    if _SITE_KEY_COLUMN not in site_metadata.columns:
        raise WorkflowBoundaryError(
            "signalome interpreter requires dataset.site_metadata.site_key"
        )
    metadata = site_metadata.copy(deep=False)
    metadata.index = pd.Index(
        metadata.index.astype(str).tolist(), name=metadata.index.name
    )
    site_keys = metadata.loc[:, _SITE_KEY_COLUMN].fillna("").astype(str).str.strip()
    mapping = dict(zip(metadata.index.tolist(), site_keys.tolist(), strict=False))
    mapped: list[str] = []
    for value in source_index.astype(str).tolist():
        resolved = mapping.get(value, "")
        if resolved == "":
            if allow_unmapped:
                mapped.append(value)
                continue
            raise WorkflowBoundaryError(
                "signalome interpreter cannot resolve site_key for diagnostic site "
                f"index value {value!r}"
            )
        mapped.append(resolved)
    return pd.Index(
        mapped,
        name=_resolved_index_name_for_site_keys(
            resolved_site_keys=mapped,
            source_site_index=source_index,
        ),
    )


def _resolved_index_name_for_site_keys(
    *,
    resolved_site_keys: list[str],
    source_site_index: pd.Index,
) -> str:
    source_values = source_site_index.astype(str).tolist()
    if len(resolved_site_keys) == len(source_values) and all(
        resolved == source
        for resolved, source in zip(resolved_site_keys, source_values, strict=True)
    ):
        source_name = source_site_index.name
        if source_name is not None:
            return str(source_name)
    return _SITE_KEY_COLUMN


__all__ = ["SignalomeWorkflowInterpreter"]
