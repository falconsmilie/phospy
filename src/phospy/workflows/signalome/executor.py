"""Internal executor for the signalome workflow."""

from __future__ import annotations

import pandas as pd

from phospy.api.results import SignalomeWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError, WorkflowStageError
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

    _MODULE_ID_COLUMN = "module_id"
    _NETWORK_SEAM = "signalome.executor.network"
    _KINASE_SUPPORT_SEAM = "signalome.executor.kinase_support"
    _MODULE_CONSTRUCTION_SEAM = "signalome.executor.module_construction"

    def run(self, request: ResolvedSignalomeWorkflowRequest) -> SignalomeWorkflowResult:
        substrate_support_cutoff = float(request.config.substrate_support_cutoff)
        network_correlation_threshold = float(
            request.config.network_correlation_threshold
        )
        try:
            module_assignments = build_module_assignments(
                prediction_matrix=request.prediction_matrix,
                site_to_protein=request.site_to_protein,
            )
        except WorkflowStageError as exc:
            self._raise_boundary_error(
                seam=self._MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure interpreted prediction inputs provide unique site IDs and "
                    "resolvable site-to-protein assignments"
                ),
                prediction_sites=int(request.prediction_matrix.shape[0]),
                prediction_kinases=int(request.prediction_matrix.shape[1]),
                substrate_support_cutoff=substrate_support_cutoff,
                network_correlation_threshold=network_correlation_threshold,
                stage_error=str(exc),
            )
        kinase_substrates = select_kinase_substrates(
            prediction_matrix=request.prediction_matrix,
            cutoff=substrate_support_cutoff,
        )
        support_counts = self._summarize_support(kinase_substrates)
        if support_counts["supported_kinases"] == 0:
            self._raise_boundary_error(
                seam=self._KINASE_SUPPORT_SEAM,
                next_action=(
                    "lower substrate_support_cutoff or ensure prediction scores "
                    "include kinase-site support above the configured threshold"
                ),
                prediction_sites=int(request.prediction_matrix.shape[0]),
                prediction_kinases=int(request.prediction_matrix.shape[1]),
                supported_sites=support_counts["supported_sites"],
                supported_kinases=support_counts["supported_kinases"],
                substrate_support_cutoff=substrate_support_cutoff,
            )

        try:
            signalome_modules = build_signalome_module_table(
                module_assignments=module_assignments,
                kinase_substrates=kinase_substrates,
                kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
            )
        except WorkflowStageError as exc:
            self._raise_boundary_error(
                seam=self._MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure supported substrates map to interpreted proteins so module "
                    "aggregation can be computed"
                ),
                supported_kinases=support_counts["supported_kinases"],
                supported_sites=support_counts["supported_sites"],
                prediction_sites=int(request.prediction_matrix.shape[0]),
                prediction_kinases=int(request.prediction_matrix.shape[1]),
                substrate_support_cutoff=substrate_support_cutoff,
                network_correlation_threshold=network_correlation_threshold,
                stage_error=str(exc),
            )
        module_count = int(
            module_assignments.loc[
                module_assignments.loc[:, self._MODULE_ID_COLUMN].astype("int64") > 0,
                self._MODULE_ID_COLUMN,
            ].nunique(dropna=True)
        )
        module_rows_with_support = (
            int((signalome_modules.sum(axis=1) > 0.0).sum())
            if not signalome_modules.empty
            else 0
        )
        if (
            signalome_modules.empty
            or module_rows_with_support == 0
            or (module_count < 2 and support_counts["supported_kinases"] < 2)
        ):
            self._raise_boundary_error(
                seam=self._MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "increase interpreted signal diversity and ensure at least two "
                    "supported kinases contribute non-zero module signal"
                ),
                module_count=module_count,
                module_rows_with_support=module_rows_with_support,
                supported_kinases=support_counts["supported_kinases"],
                supported_sites=support_counts["supported_sites"],
                prediction_sites=int(request.prediction_matrix.shape[0]),
                prediction_kinases=int(request.prediction_matrix.shape[1]),
                substrate_support_cutoff=substrate_support_cutoff,
                network_correlation_threshold=network_correlation_threshold,
            )

        score_variance_kinases = self._score_variance_kinases(
            request.downstream_score_matrix
        )
        try:
            network_edges, network_nodes = build_kinase_network(
                downstream_score_matrix=request.downstream_score_matrix,
                kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
                kinase_substrates=kinase_substrates,
                threshold=network_correlation_threshold,
            )
        except WorkflowStageError as exc:
            self._raise_boundary_error(
                seam=self._NETWORK_SEAM,
                next_action=(
                    "ensure interpreted score and prediction matrices share the same "
                    "kinase set and contain variable score signal"
                ),
                shared_kinases=int(request.prediction_matrix.shape[1]),
                supported_kinases=support_counts["supported_kinases"],
                downstream_score_sites=int(request.downstream_score_matrix.shape[0]),
                downstream_score_source=request.downstream_score_source,
                score_variance_kinases=score_variance_kinases,
                network_correlation_threshold=network_correlation_threshold,
                stage_error=str(exc),
            )
        if network_edges.empty:
            self._raise_boundary_error(
                seam=self._NETWORK_SEAM,
                next_action=(
                    "lower network_correlation_threshold or provide more variable "
                    "score profiles so kinase correlations can be estimated"
                ),
                shared_kinases=int(request.prediction_matrix.shape[1]),
                supported_kinases=support_counts["supported_kinases"],
                downstream_score_sites=int(request.downstream_score_matrix.shape[0]),
                downstream_score_source=request.downstream_score_source,
                score_variance_kinases=score_variance_kinases,
                network_correlation_threshold=network_correlation_threshold,
            )

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
            ),
            expanded_signalome=None,
        )

    @staticmethod
    def _summarize_support(
        kinase_substrates: dict[str, tuple[str, ...]],
    ) -> dict[str, int]:
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

    @staticmethod
    def _score_variance_kinases(downstream_score_matrix: pd.DataFrame) -> int:
        if downstream_score_matrix.empty:
            return 0
        variances = downstream_score_matrix.astype(float).var(axis=0, ddof=0)
        return int((variances > 0.0).sum())

    @staticmethod
    def _raise_boundary_error(
        *,
        seam: str,
        next_action: str,
        **details: int | float | str,
    ) -> None:
        details_text = ", ".join(f"{key}={value}" for key, value in details.items())
        raise WorkflowBoundaryError(
            "signalome workflow boundary validation failed at "
            f"seam={seam}; {details_text}; next_action={next_action}"
        )
