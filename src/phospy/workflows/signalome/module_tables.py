"""Module assignment and module summary table builders."""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.science.signalomes.clustering import ClusterSitesResult
from phospy.science.signalomes.constants import MODULE_ID_COLUMN
from phospy.science.signalomes.science import (
    build_module_assignments,
    build_signalome_module_table,
    select_kinase_substrates,
)
from phospy.workflows.signalome.component_helpers import (
    module_selection_details,
    prediction_shape_details,
    raise_boundary_error,
    summarize_support,
    support_details,
)
from phospy.workflows.signalome.component_models import (
    SignalomeExecutionMetadata,
    SignalomeModuleTableBuildResult,
    SignalomeSupportSummary,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_EXECUTOR_KINASE_SUPPORT_SEAM,
    SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)


class SignalomeModuleTableBuilder:
    """Build module-assignment and module summary tables for signalome execution."""

    _MODULE_ID_COLUMN = MODULE_ID_COLUMN

    def __init__(
        self,
        *,
        build_assignments: Callable[..., pd.DataFrame] = build_module_assignments,
        select_substrates: Callable[..., dict[str, tuple[str, ...]]] = (
            select_kinase_substrates
        ),
        build_modules: Callable[..., pd.DataFrame] = build_signalome_module_table,
    ) -> None:
        self._build_assignments = build_assignments
        self._select_substrates = select_substrates
        self._build_modules = build_modules

    def run(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        protein_modules: pd.Series,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeModuleTableBuildResult:
        module_assignments = self._build_module_assignments(
            request=request,
            config=config,
            clustering_result=clustering_result,
            protein_modules=protein_modules,
            execution_metadata=execution_metadata,
        )
        support_summary = self._select_supported_substrates(
            request=request,
            config=config,
            execution_metadata=execution_metadata,
        )
        support_counts = support_summary.support_counts
        try:
            signalome_modules = self._build_modules(
                module_assignments=module_assignments,
                kinase_substrates=support_summary.kinase_substrates,
                kinase_order=request.prediction_matrix.columns.astype(str).tolist(),
                assignment_policy=config.assignment_policy,
            )
        except WorkflowStageError as exc:
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure supported substrates map to interpreted proteins so module "
                    "aggregation can be computed"
                ),
                **support_details(support_counts),
                **prediction_shape_details(execution_metadata),
                substrate_support_cutoff=config.substrate_support_cutoff,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
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
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "increase interpreted signal diversity and ensure at least two "
                    "supported kinases contribute non-zero module signal"
                ),
                module_count=module_count,
                module_rows_with_support=module_rows_with_support,
                **support_details(support_counts),
                **prediction_shape_details(execution_metadata),
                substrate_support_cutoff=config.substrate_support_cutoff,
                network_correlation_threshold=config.network_correlation_threshold,
                selected_module_count=module_selection_details(clustering_result)[
                    "selected_module_count"
                ],
            )
        return SignalomeModuleTableBuildResult(
            module_assignments=module_assignments,
            signalome_modules=signalome_modules,
            module_count=module_count,
            support_summary=support_summary,
        )

    def _build_module_assignments(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        clustering_result: ClusterSitesResult,
        protein_modules: pd.Series,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> pd.DataFrame:
        try:
            return self._build_assignments(
                prediction_matrix=request.prediction_matrix,
                site_to_protein=request.site_to_protein,
                site_metadata=request.dataset._borrow_site_metadata_frame(),
                protein_modules=protein_modules,
            )
        except WorkflowStageError as exc:
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
                next_action=(
                    "ensure interpreted prediction inputs provide unique site_key rows and "
                    "resolvable site-to-protein assignments"
                ),
                **prediction_shape_details(execution_metadata),
                substrate_support_cutoff=config.substrate_support_cutoff,
                network_correlation_threshold=config.network_correlation_threshold,
                network_policy=config.network_policy,
                **module_selection_details(clustering_result),
                stage_error=str(exc),
            )

    def _select_supported_substrates(
        self,
        *,
        request: ResolvedSignalomeWorkflowRequest,
        config: ResolvedSignalomeExecutionConfig,
        execution_metadata: SignalomeExecutionMetadata,
    ) -> SignalomeSupportSummary:
        kinase_substrates = self._select_substrates(
            prediction_matrix=request.prediction_matrix,
            cutoff=config.substrate_support_cutoff,
        )
        support_counts = summarize_support(kinase_substrates)
        if support_counts["supported_kinases"] == 0:
            raise_boundary_error(
                seam=SIGNALOME_EXECUTOR_KINASE_SUPPORT_SEAM,
                next_action=(
                    "lower substrate_support_cutoff or ensure prediction scores "
                    "include kinase-site support above the configured threshold"
                ),
                **prediction_shape_details(execution_metadata),
                **support_details(support_counts),
                substrate_support_cutoff=config.substrate_support_cutoff,
            )
        return SignalomeSupportSummary(
            kinase_substrates=kinase_substrates,
            support_counts=support_counts,
        )


__all__ = ["SignalomeModuleTableBuilder"]
