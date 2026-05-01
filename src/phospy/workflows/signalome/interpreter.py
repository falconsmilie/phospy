"""Internal interpreter for signalome workflow requests."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy.api.configs import (
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    SignalomeScorePreconditioningPolicy,
)
from phospy.api.requests import SignalomeWorkflowRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.prediction.scoring import select_downstream_score_matrix
from phospy.signalomes.constants import KINASE_COLUMN, PROTEIN_COLUMN, SITE_ID_COLUMN
from phospy.signalomes.models import (
    SignalomeAlignmentDiagnostics,
    SignalomeAlignmentInputDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_INTERPRETER_KINASE_OVERLAP_SEAM,
    SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM,
    SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
    SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM,
    SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA,
    SIGNALOME_WORKFLOW_BOUNDARY_MESSAGE_PREFIX,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


class SignalomeWorkflowInterpreter:
    """Resolve a validated signalome request into execution-ready inputs."""

    _SITE_ID_COLUMN = SITE_ID_COLUMN
    _KINASE_COLUMN = KINASE_COLUMN
    _PROTEIN_COLUMN = PROTEIN_COLUMN
    _SITE_ALIGNMENT_SEAM = SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM
    _KINASE_OVERLAP_SEAM = SIGNALOME_INTERPRETER_KINASE_OVERLAP_SEAM
    _PROTEIN_MAPPING_SEAM = SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM
    _SCORE_PRECONDITIONING_SEAM = SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM
    _REASON_MISSING_FROM_DATASET = "missing_from_dataset"
    _REASON_MISSING_FROM_PREDICTION_SCORES = "missing_from_prediction_scores"
    _REASON_MISSING_FROM_DOWNSTREAM_SCORES = "missing_from_downstream_scores"
    _REASON_MISSING_FROM_KINASE_SUPPORT = "missing_kinase_support"
    _REASON_MISSING_PROTEIN_IDENTIFIER = "missing_protein_identifier"
    _REASON_REMOVED_BY_SCORE_PRECONDITIONING = "removed_by_score_preconditioning"
    _REASON_REMOVED_BY_VALIDATION_POLICY = "removed_by_validation_policy"

    def run(
        self, request: SignalomeWorkflowRequest
    ) -> ResolvedSignalomeWorkflowRequest:
        scoring_result = request.kinase_result.scoring_result
        downstream_score_matrix, downstream_score_source = (
            select_downstream_score_matrix(
                profile_scores=scoring_result.profile_scores,
                rank_weighted_fusion_scores=scoring_result.rank_weighted_fusion_scores,
            )
        )
        resolved_downstream_score_matrix = self._as_aligned_numeric_frame(
            downstream_score_matrix,
            index_name=self._SITE_ID_COLUMN,
            columns_name=self._KINASE_COLUMN,
        )
        resolved_prediction_matrix = self._as_aligned_numeric_frame(
            request.kinase_result.prediction_result.pred_mat,
            index_name=self._SITE_ID_COLUMN,
            columns_name=self._KINASE_COLUMN,
        )
        dataset_site_index = pd.Index(
            request.kinase_result.dataset.phospho.index.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        aligned_site_index = self._resolve_shared_site_index(
            dataset_sites=dataset_site_index,
            prediction_sites=resolved_prediction_matrix.index,
            score_sites=resolved_downstream_score_matrix.index,
        )
        aligned_kinase_index = self._resolve_shared_kinase_index(
            prediction_kinases=resolved_prediction_matrix.columns,
            score_kinases=resolved_downstream_score_matrix.columns,
        )
        aligned_prediction_matrix = resolved_prediction_matrix.loc[
            aligned_site_index, aligned_kinase_index
        ]
        aligned_downstream_score_matrix = resolved_downstream_score_matrix.loc[
            aligned_site_index, aligned_kinase_index
        ]
        execution_config = self._resolve_execution_config(request)
        (
            preconditioned_downstream_score_matrix,
            score_preconditioning_diagnostics,
        ) = self._precondition_downstream_score_matrix(
            aligned_downstream_score_matrix,
            policy=execution_config.score_preconditioning_policy,
        )
        if preconditioned_downstream_score_matrix.empty:
            self._raise_boundary_error(
                seam=self._SCORE_PRECONDITIONING_SEAM,
                next_action=(
                    "ensure upstream downstream scores retain at least one "
                    "interpretable site after score preconditioning"
                ),
                aligned_score_sites=int(aligned_downstream_score_matrix.shape[0]),
                aligned_score_kinases=int(aligned_downstream_score_matrix.shape[1]),
                dropped_all_missing_row_count=(
                    score_preconditioning_diagnostics.dropped_all_missing_row_count
                ),
                retained_row_count=score_preconditioning_diagnostics.retained_row_count,
                score_preconditioning_policy=(
                    execution_config.score_preconditioning_policy
                ),
            )
        retained_site_index = preconditioned_downstream_score_matrix.index
        retained_prediction_matrix = aligned_prediction_matrix.loc[retained_site_index]
        site_to_protein = self._resolve_site_to_protein(
            dataset=request.kinase_result.dataset,
            site_index=retained_site_index,
            removed_by_score_preconditioning_count=int(
                aligned_site_index.size - retained_site_index.size
            ),
        )
        alignment_diagnostics = self._build_alignment_diagnostics(
            dataset_sites=dataset_site_index,
            prediction_sites=resolved_prediction_matrix.index,
            score_sites=resolved_downstream_score_matrix.index,
            shared_sites=aligned_site_index,
            retained_sites=retained_site_index,
            prediction_kinases=resolved_prediction_matrix.columns,
            score_kinases=resolved_downstream_score_matrix.columns,
            shared_kinases=aligned_kinase_index,
            interpreted_protein_sites=aligned_site_index,
            retained_protein_sites=retained_site_index,
        )
        return ResolvedSignalomeWorkflowRequest(
            dataset=request.kinase_result.dataset,
            kinase_result=request.kinase_result,
            execution_config=execution_config,
            downstream_score_matrix=preconditioned_downstream_score_matrix,
            downstream_score_source=downstream_score_source,
            prediction_matrix=retained_prediction_matrix,
            site_to_protein=site_to_protein,
            score_preconditioning_diagnostics=score_preconditioning_diagnostics,
            alignment_diagnostics=alignment_diagnostics,
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
        )

    @staticmethod
    def _as_aligned_numeric_frame(
        frame: pd.DataFrame,
        *,
        index_name: str,
        columns_name: str,
    ) -> pd.DataFrame:
        resolved = frame.astype(float)
        resolved.index = pd.Index(resolved.index.astype(str), name=index_name)
        resolved.columns = pd.Index(resolved.columns.astype(str), name=columns_name)
        return resolved

    def _resolve_shared_site_index(
        self,
        *,
        dataset_sites: pd.Index,
        prediction_sites: pd.Index,
        score_sites: pd.Index,
    ) -> pd.Index:
        dataset_site_index = pd.Index(
            dataset_sites.astype(str), name=self._SITE_ID_COLUMN
        )
        prediction_site_index = pd.Index(
            prediction_sites.astype(str), name=self._SITE_ID_COLUMN
        )
        score_site_index = pd.Index(score_sites.astype(str), name=self._SITE_ID_COLUMN)
        prediction_site_set = set(prediction_site_index.tolist())
        score_site_set = set(score_site_index.tolist())
        shared_sites = [
            site_id
            for site_id in dataset_site_index
            if site_id in prediction_site_set and site_id in score_site_set
        ]
        if not shared_sites:
            self._raise_boundary_error(
                seam=self._SITE_ALIGNMENT_SEAM,
                next_action=(
                    "ensure prediction and scoring outputs share phosphosite IDs with "
                    "kinase_result.dataset.phospho.index"
                ),
                dataset_sites=int(dataset_site_index.size),
                prediction_sites=int(prediction_site_index.size),
                score_sites=int(score_site_index.size),
                shared_sites=0,
            )
        return pd.Index(shared_sites, name=self._SITE_ID_COLUMN)

    def _resolve_shared_kinase_index(
        self,
        *,
        prediction_kinases: pd.Index,
        score_kinases: pd.Index,
    ) -> pd.Index:
        prediction_kinase_index = pd.Index(
            prediction_kinases.astype(str), name=self._KINASE_COLUMN
        )
        score_kinase_index = pd.Index(
            score_kinases.astype(str), name=self._KINASE_COLUMN
        )
        score_kinase_set = set(score_kinase_index.tolist())
        shared_kinases = [
            kinase for kinase in prediction_kinase_index if kinase in score_kinase_set
        ]
        if not shared_kinases:
            self._raise_boundary_error(
                seam=self._KINASE_OVERLAP_SEAM,
                next_action=(
                    "rerun kinase workflow so scoring_result and "
                    "prediction_result are generated from the same kinase lane"
                ),
                prediction_kinases=int(prediction_kinase_index.size),
                score_kinases=int(score_kinase_index.size),
                shared_kinases=0,
            )
        return pd.Index(shared_kinases, name=self._KINASE_COLUMN)

    def _resolve_site_to_protein(
        self,
        *,
        dataset: AnalysisReadyPhosphoDataset,
        site_index: pd.Index,
        removed_by_score_preconditioning_count: int,
    ) -> pd.Series:
        metadata = dataset.site_metadata
        if self._PROTEIN_COLUMN not in metadata.columns:
            self._raise_boundary_error(
                seam=self._PROTEIN_MAPPING_SEAM,
                next_action=(
                    "populate dataset.site_metadata.protein_id with explicit protein "
                    "identifiers for retained signalome sites after score "
                    "preconditioning"
                ),
                protein_resolution_source=SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA,
                interpreted_sites=int(site_index.size),
                resolved_protein_sites=0,
                unresolved_protein_sites=int(site_index.size),
                removed_by_score_preconditioning=(
                    int(removed_by_score_preconditioning_count)
                ),
                retained_with_missing_protein_id=int(site_index.size),
                retained_and_valid=0,
            )
        resolved = (
            metadata.reindex(site_index)
            .loc[:, self._PROTEIN_COLUMN]
            .fillna("")
            .astype(str)
            .str.strip()
        )
        resolution_source = SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA
        unresolved_mask = resolved.astype(str).str.strip() == ""
        if unresolved_mask.any():
            resolved_sites = int((~unresolved_mask).sum())
            unresolved_sites = resolved.loc[unresolved_mask].index.astype(str).tolist()
            self._raise_boundary_error(
                seam=self._PROTEIN_MAPPING_SEAM,
                next_action=(
                    "populate dataset.site_metadata.protein_id with explicit protein "
                    "identifiers for retained signalome sites after score "
                    "preconditioning"
                ),
                protein_resolution_source=resolution_source,
                interpreted_sites=int(site_index.size),
                resolved_protein_sites=resolved_sites,
                unresolved_protein_sites=int(unresolved_mask.sum()),
                retained_with_missing_protein_id=int(unresolved_mask.sum()),
                retained_and_valid=resolved_sites,
                removed_by_score_preconditioning=(
                    int(removed_by_score_preconditioning_count)
                ),
                missing_protein_id_sites=unresolved_sites,
            )
        resolved.index = site_index.copy()
        resolved.name = self._PROTEIN_COLUMN
        return resolved.astype(str)

    def _precondition_downstream_score_matrix(
        self,
        score_matrix: pd.DataFrame,
        *,
        policy: SignalomeScorePreconditioningPolicy,
    ) -> tuple[pd.DataFrame, SignalomeScorePreconditioningDiagnostics]:
        """Prepare aligned upstream scores for score-driven signalome stages.

        All-missing rows are treated as unsupported evidence and handled by the
        configured policy.
        Partially missing rows are retained for pairwise-complete correlation.
        """
        if score_matrix.empty:
            return (
                score_matrix,
                self._score_preconditioning_diagnostics(
                    input_row_count=0,
                    dropped_all_missing_row_count=0,
                    retained_row_count=0,
                    policy=policy,
                ),
            )
        score_values = score_matrix.to_numpy(dtype=float, copy=False)
        infinite_mask = np.isinf(score_values)
        if infinite_mask.any():
            self._raise_boundary_error(
                seam=self._SCORE_PRECONDITIONING_SEAM,
                next_action=(
                    "rerun kinase workflow and ensure scoring outputs contain "
                    "finite values only"
                ),
                aligned_score_sites=int(score_matrix.shape[0]),
                aligned_score_kinases=int(score_matrix.shape[1]),
                infinite_score_entries=int(infinite_mask.sum()),
            )
        supported_rows = score_matrix.notna().any(axis=1)
        supported_row_mask = supported_rows.to_numpy(dtype=bool, copy=False)
        input_row_count = int(score_matrix.shape[0])
        retained_row_count = int(supported_row_mask.sum())
        dropped_all_missing_row_count = int(input_row_count - retained_row_count)
        diagnostics = self._score_preconditioning_diagnostics(
            input_row_count=input_row_count,
            dropped_all_missing_row_count=dropped_all_missing_row_count,
            retained_row_count=retained_row_count,
            policy=policy,
        )
        if policy == SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP:
            if dropped_all_missing_row_count > 0:
                self._raise_boundary_error(
                    seam=self._SCORE_PRECONDITIONING_SEAM,
                    next_action=(
                        "set "
                        "config.validation.score_preconditioning_policy="
                        "'allow_and_report' "
                        "to proceed with explicit row dropping, or ensure upstream "
                        "downstream scores contain non-missing support for every "
                        "interpreted site"
                    ),
                    aligned_score_sites=input_row_count,
                    aligned_score_kinases=int(score_matrix.shape[1]),
                    dropped_all_missing_row_count=dropped_all_missing_row_count,
                    retained_row_count=retained_row_count,
                    score_preconditioning_policy=policy,
                )
            return score_matrix, diagnostics
        if policy != SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT:
            self._raise_boundary_error(
                seam=self._SCORE_PRECONDITIONING_SEAM,
                next_action=(
                    "use a supported score preconditioning policy from "
                    "SignalomeConfig.validation.score_preconditioning_policy"
                ),
                score_preconditioning_policy=policy,
            )
        if supported_row_mask.all():
            return score_matrix, diagnostics
        return score_matrix.iloc[supported_row_mask, :], diagnostics

    @staticmethod
    def _score_preconditioning_diagnostics(
        *,
        input_row_count: int,
        dropped_all_missing_row_count: int,
        retained_row_count: int,
        policy: SignalomeScorePreconditioningPolicy,
    ) -> SignalomeScorePreconditioningDiagnostics:
        return SignalomeScorePreconditioningDiagnostics(
            input_row_count=int(input_row_count),
            dropped_all_missing_row_count=int(dropped_all_missing_row_count),
            retained_row_count=int(retained_row_count),
            policy=policy,
        )

    def _build_alignment_diagnostics(
        self,
        *,
        dataset_sites: pd.Index,
        prediction_sites: pd.Index,
        score_sites: pd.Index,
        shared_sites: pd.Index,
        retained_sites: pd.Index,
        prediction_kinases: pd.Index,
        score_kinases: pd.Index,
        shared_kinases: pd.Index,
        interpreted_protein_sites: pd.Index,
        retained_protein_sites: pd.Index,
    ) -> SignalomeAlignmentDiagnostics:
        dataset_site_index = pd.Index(
            dataset_sites.astype(str), name=self._SITE_ID_COLUMN
        )
        prediction_site_index = pd.Index(
            prediction_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        score_site_index = pd.Index(score_sites.astype(str), name=self._SITE_ID_COLUMN)
        shared_site_index = pd.Index(
            shared_sites.astype(str), name=self._SITE_ID_COLUMN
        )
        retained_site_index = pd.Index(
            retained_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        dataset_sites_diag = self._build_site_alignment_input_diagnostics(
            provided=dataset_site_index,
            retained=retained_site_index,
            missing_from_first=(
                prediction_site_index,
                self._REASON_MISSING_FROM_PREDICTION_SCORES,
            ),
            missing_from_second=(
                score_site_index,
                self._REASON_MISSING_FROM_DOWNSTREAM_SCORES,
            ),
        )
        prediction_sites_diag = self._build_site_alignment_input_diagnostics(
            provided=prediction_site_index,
            retained=retained_site_index,
            missing_from_first=(
                dataset_site_index,
                self._REASON_MISSING_FROM_DATASET,
            ),
            missing_from_second=(
                score_site_index,
                self._REASON_MISSING_FROM_DOWNSTREAM_SCORES,
            ),
        )
        downstream_sites_diag = self._build_site_alignment_input_diagnostics(
            provided=score_site_index,
            retained=retained_site_index,
            missing_from_first=(
                dataset_site_index,
                self._REASON_MISSING_FROM_DATASET,
            ),
            missing_from_second=(
                prediction_site_index,
                self._REASON_MISSING_FROM_PREDICTION_SCORES,
            ),
        )

        prediction_kinase_index = pd.Index(
            prediction_kinases.astype(str),
            name=self._KINASE_COLUMN,
        )
        score_kinase_index = pd.Index(
            score_kinases.astype(str), name=self._KINASE_COLUMN
        )
        shared_kinase_index = pd.Index(
            shared_kinases.astype(str), name=self._KINASE_COLUMN
        )
        provided_kinases = prediction_kinase_index.append(
            score_kinase_index.difference(prediction_kinase_index)
        )
        kinase_reasons = {
            self._REASON_MISSING_FROM_PREDICTION_SCORES: 0,
            self._REASON_MISSING_FROM_DOWNSTREAM_SCORES: 0,
            self._REASON_MISSING_FROM_KINASE_SUPPORT: 0,
        }
        retained_kinase_set = set(shared_kinase_index.tolist())
        prediction_kinase_set = set(prediction_kinase_index.tolist())
        score_kinase_set = set(score_kinase_index.tolist())
        for kinase in provided_kinases:
            kinase_id = str(kinase)
            if kinase_id in retained_kinase_set:
                continue
            if kinase_id not in prediction_kinase_set:
                kinase_reasons[self._REASON_MISSING_FROM_PREDICTION_SCORES] += 1
            elif kinase_id not in score_kinase_set:
                kinase_reasons[self._REASON_MISSING_FROM_DOWNSTREAM_SCORES] += 1
            else:
                kinase_reasons[self._REASON_MISSING_FROM_KINASE_SUPPORT] += 1
        kinases_diag = SignalomeAlignmentInputDiagnostics(
            provided_count=int(provided_kinases.size),
            retained_count=int(shared_kinase_index.size),
            dropped_count=int(provided_kinases.size - shared_kinase_index.size),
            dropped_reasons=kinase_reasons,
        )

        interpreted_protein_index = pd.Index(
            interpreted_protein_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        retained_protein_index = pd.Index(
            retained_protein_sites.astype(str),
            name=self._SITE_ID_COLUMN,
        )
        protein_dropped = int(
            interpreted_protein_index.size - retained_protein_index.size
        )
        protein_diag = SignalomeAlignmentInputDiagnostics(
            provided_count=int(interpreted_protein_index.size),
            retained_count=int(retained_protein_index.size),
            dropped_count=protein_dropped,
            dropped_reasons={
                self._REASON_REMOVED_BY_SCORE_PRECONDITIONING: protein_dropped,
                self._REASON_MISSING_PROTEIN_IDENTIFIER: 0,
                self._REASON_REMOVED_BY_VALIDATION_POLICY: 0,
            },
        )
        _ = shared_site_index
        return SignalomeAlignmentDiagnostics(
            dataset_sites=dataset_sites_diag,
            prediction_score_sites=prediction_sites_diag,
            downstream_score_sites=downstream_sites_diag,
            kinases=kinases_diag,
            protein_identifiers=protein_diag,
        )

    def _build_site_alignment_input_diagnostics(
        self,
        *,
        provided: pd.Index,
        retained: pd.Index,
        missing_from_first: tuple[pd.Index, str],
        missing_from_second: tuple[pd.Index, str],
    ) -> SignalomeAlignmentInputDiagnostics:
        provided_index = pd.Index(provided.astype(str), name=self._SITE_ID_COLUMN)
        retained_index = pd.Index(retained.astype(str), name=self._SITE_ID_COLUMN)
        first_index, first_reason = missing_from_first
        second_index, second_reason = missing_from_second
        first_site_set = set(pd.Index(first_index.astype(str)).tolist())
        second_site_set = set(pd.Index(second_index.astype(str)).tolist())
        retained_site_set = set(retained_index.tolist())
        reasons = {
            first_reason: 0,
            second_reason: 0,
            self._REASON_REMOVED_BY_SCORE_PRECONDITIONING: 0,
            self._REASON_REMOVED_BY_VALIDATION_POLICY: 0,
        }
        for site_id in provided_index:
            site = str(site_id)
            if site in retained_site_set:
                continue
            if site not in first_site_set:
                reasons[first_reason] += 1
            elif site not in second_site_set:
                reasons[second_reason] += 1
            else:
                reasons[self._REASON_REMOVED_BY_SCORE_PRECONDITIONING] += 1
        return SignalomeAlignmentInputDiagnostics(
            provided_count=int(provided_index.size),
            retained_count=int(retained_index.size),
            dropped_count=int(provided_index.size - retained_index.size),
            dropped_reasons=reasons,
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
