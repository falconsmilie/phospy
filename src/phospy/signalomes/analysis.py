from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..errors import InputCompatibilityError
from ..validation.requests.signalome import (
    SignalomeInputs,
    SignalomeRequest,
)
from .assignments import (
    build_expanded_signalomes,
    build_kinase_module_relationship_table,
    build_protein_assignment_table,
    build_signalome_module_table,
    build_site_assignments,
    derive_protein_modules,
    select_kinase_substrates,
)
from .clustering import (
    SignalomeModuleSelectionDiagnostics,
    SignalomeModuleSelectionPolicy,
    cluster_sites_with_diagnostics,
)
from .results import (
    ExpandedSignalome,
    SignalomeAssignments,
    SignalomeKinaseNetwork,
    SignalomeModules,
    SignalomeResult,
)

__all__ = ["SignalomeRunner", "build_signalome_result", "execute_signalome_inputs"]


@dataclass(slots=True)
class SignalomeRunner:
    """Construct a signalome result from trusted aligned signalome inputs."""

    @dataclass(frozen=True, slots=True)
    class _ClusteringStage:
        site_clusters: pd.Series
        protein_modules: pd.Series
        module_selection_diagnostics: SignalomeModuleSelectionDiagnostics

    @dataclass(frozen=True, slots=True)
    class _AssignmentsStage:
        site_assignments: pd.DataFrame
        kinase_substrates: dict[str, tuple[str, ...]]

    @dataclass(frozen=True, slots=True)
    class _NetworkStage:
        kinase_network: pd.DataFrame
        kinase_correlation_matrix: pd.DataFrame
        network: SignalomeKinaseNetwork

    @dataclass(frozen=True, slots=True)
    class _ModulesStage:
        module_table: pd.DataFrame
        protein_assignments: pd.DataFrame
        kinase_module_relationships: pd.DataFrame

    def execute(self, inputs: SignalomeInputs) -> SignalomeResult:
        clustering_stage = self._cluster_sites(inputs)
        assignments_stage = self._assign_sites(inputs, clustering_stage)
        network_stage = self._build_network(inputs, assignments_stage)
        modules_stage = self._build_modules(assignments_stage)
        expanded_signalomes = self._expand_signalomes(
            inputs=inputs,
            assignments_stage=assignments_stage,
            network_stage=network_stage,
            modules_stage=modules_stage,
        )

        return SignalomeResult(
            scoring_matrix=inputs.scoring_matrix,
            pred_mat=inputs.pred_mat,
            expression_matrix=inputs.expression_matrix,
            modules=SignalomeModules(
                module_table=modules_stage.module_table,
                kinase_module_relationships=modules_stage.kinase_module_relationships,
            ),
            assignments=SignalomeAssignments(
                site_assignments=assignments_stage.site_assignments,
                protein_assignments=modules_stage.protein_assignments,
            ),
            network=network_stage.network,
            kinase_substrate_map=assignments_stage.kinase_substrates,
            expanded_signalomes=expanded_signalomes,
            module_selection_diagnostics=(
                clustering_stage.module_selection_diagnostics
            ),
        )

    def _cluster_sites(self, inputs: SignalomeInputs) -> _ClusteringStage:
        clustering_result = cluster_sites_with_diagnostics(
            scoring_matrix=inputs.scoring_matrix,
            requested_module_count=inputs.module_count,
            policy=inputs.module_selection_policy,
        )
        protein_modules = derive_protein_modules(
            site_clusters=clustering_result.site_clusters,
            site_to_protein=inputs.site_to_protein,
        )

        return self._ClusteringStage(
            site_clusters=clustering_result.site_clusters,
            protein_modules=protein_modules,
            module_selection_diagnostics=clustering_result.module_selection_diagnostics,
        )

    def _assign_sites(
        self,
        inputs: SignalomeInputs,
        clustering_stage: _ClusteringStage,
    ) -> _AssignmentsStage:
        site_assignments = build_site_assignments(
            pred_mat=inputs.pred_mat,
            protein_modules=clustering_stage.protein_modules,
            site_to_protein=inputs.site_to_protein,
        )
        kinase_substrates = select_kinase_substrates(
            pred_mat=inputs.pred_mat,
            cutoff=inputs.signalome_cutoff,
        )

        return self._AssignmentsStage(
            site_assignments=site_assignments,
            kinase_substrates=kinase_substrates,
        )

    def _build_network(
        self,
        inputs: SignalomeInputs,
        assignments_stage: _AssignmentsStage,
    ) -> _NetworkStage:
        kinase_network, kinase_correlation_matrix = build_kinase_network(
            scoring_matrix=inputs.scoring_matrix,
            threshold=inputs.kinase_network_threshold,
        )
        network = build_kinase_network_view(
            kinase_network=kinase_network,
            kinase_correlation_matrix=kinase_correlation_matrix,
            kinase_substrates=assignments_stage.kinase_substrates,
        )

        return self._NetworkStage(
            kinase_network=kinase_network,
            kinase_correlation_matrix=kinase_correlation_matrix,
            network=network,
        )

    def _build_modules(self, assignments_stage: _AssignmentsStage) -> _ModulesStage:
        module_table = build_signalome_module_table(
            site_assignments=assignments_stage.site_assignments,
            kinase_substrates=assignments_stage.kinase_substrates,
        )
        protein_assignments = build_protein_assignment_table(
            site_assignments=assignments_stage.site_assignments,
        )
        kinase_module_relationships = build_kinase_module_relationship_table(
            module_table=module_table,
        )

        return self._ModulesStage(
            module_table=module_table,
            protein_assignments=protein_assignments,
            kinase_module_relationships=kinase_module_relationships,
        )

    def _expand_signalomes(
        self,
        *,
        inputs: SignalomeInputs,
        assignments_stage: _AssignmentsStage,
        network_stage: _NetworkStage,
        modules_stage: _ModulesStage,
    ) -> dict[str, ExpandedSignalome]:
        return build_expanded_signalomes(
            kinases_of_interest=inputs.kinases_of_interest,
            kinase_network=network_stage.network.neighbor_map,
            kinase_substrates=assignments_stage.kinase_substrates,
            signalome_modules=modules_stage.module_table,
            site_assignments=assignments_stage.site_assignments,
            expression_matrix=inputs.expression_matrix,
            min_kinase_module_share_percent=inputs.min_kinase_module_share_percent,
        )


def execute_signalome_inputs(inputs: SignalomeInputs) -> SignalomeResult:
    """Build a signalome result from trusted signalome inputs."""

    return SignalomeRunner().execute(inputs)


def build_signalome_result(
    *,
    scoring_matrix: pd.DataFrame,
    pred_mat: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    kinases_of_interest: Sequence[str],
    site_to_protein: Mapping[str, str] | pd.Series | None = None,
    kinase_network_threshold: float = 0.9,
    signalome_cutoff: float = 0.5,
    module_count: int | None = None,
    min_kinase_module_share_percent: float = 1.0,
    module_selection_policy: SignalomeModuleSelectionPolicy | None = None,
) -> SignalomeResult:
    """Build a structured signalome result from validated aligned inputs."""

    request = SignalomeRequest.validate_request(
        kinases_of_interest=kinases_of_interest,
        site_to_protein=(None if site_to_protein is None else dict(site_to_protein)),
        kinase_network_threshold=kinase_network_threshold,
        signalome_cutoff=signalome_cutoff,
        module_count=module_count,
        min_kinase_module_share_percent=min_kinase_module_share_percent,
        module_selection_policy=module_selection_policy,
    )
    inputs = SignalomeInputs.from_trusted_inputs(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=request.kinases_of_interest,
        site_to_protein=request.site_to_protein,
        kinase_network_threshold=request.kinase_network_threshold,
        signalome_cutoff=request.signalome_cutoff,
        module_count=request.module_count,
        min_kinase_module_share_percent=request.min_kinase_module_share_percent,
        module_selection_policy=request.module_selection_policy,
        scoring_context="scoring_matrix",
        pred_mat_context="pred_mat",
        expression_context="expression_matrix",
    )

    return execute_signalome_inputs(inputs)


def build_kinase_network(
    *,
    scoring_matrix: pd.DataFrame,
    threshold: float = 0.9,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a thresholded kinase-kinase co-correlation adjacency matrix."""

    correlation_matrix = scoring_matrix.corr(method="pearson")
    thresholded_network = correlation_matrix.where(
        correlation_matrix >= threshold,
        0.0,
    ).copy()

    for kinase in thresholded_network.index:
        thresholded_network.loc[kinase, kinase] = 0.0

    return thresholded_network, correlation_matrix


def build_kinase_network_view(
    *,
    kinase_network: pd.DataFrame,
    kinase_correlation_matrix: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
) -> SignalomeKinaseNetwork:
    """Wrap derived network tables and neighbour lookup in a structured result view."""

    kinase_index = pd.Index(kinase_correlation_matrix.index.astype(str), name="kinase")
    aligned_network = kinase_network.loc[kinase_index, kinase_index]
    aligned_correlation_matrix = kinase_correlation_matrix.loc[
        kinase_index, kinase_index
    ]
    kinase_names = kinase_index.to_numpy(dtype=object, copy=False)

    network_values = aligned_network.to_numpy(dtype=float, copy=False)
    positive_edge_mask = network_values > 0.0
    diagonal_positions = np.arange(len(kinase_names), dtype=int)
    positive_edge_mask[diagonal_positions, diagonal_positions] = False

    neighbor_map = {
        str(kinase): tuple(str(neighbor) for neighbor in kinase_names[row_mask])
        for kinase, row_mask in zip(kinase_names, positive_edge_mask, strict=True)
    }

    node_table = pd.DataFrame(
        {
            "degree": positive_edge_mask.sum(axis=1, dtype=int),
            "n_substrates": np.asarray(
                [
                    len(tuple(kinase_substrates.get(str(kinase), ())))
                    for kinase in kinase_names
                ],
                dtype=int,
            ),
        },
        index=kinase_index,
    ).astype({"degree": int, "n_substrates": int})

    upper_triangle_source, upper_triangle_target = np.triu_indices(
        len(kinase_names),
        k=1,
    )
    edge_mask = positive_edge_mask[upper_triangle_source, upper_triangle_target]
    edge_table = pd.DataFrame(
        {
            "source_kinase": kinase_names[upper_triangle_source[edge_mask]],
            "target_kinase": kinase_names[upper_triangle_target[edge_mask]],
            "correlation": aligned_correlation_matrix.to_numpy(dtype=float, copy=False)[
                upper_triangle_source[edge_mask],
                upper_triangle_target[edge_mask],
            ],
        }
    )
    if edge_table.empty:
        edge_table = pd.DataFrame(
            columns=["source_kinase", "target_kinase", "correlation"]
        )
    edge_table = edge_table.astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
        }
    )
    edge_table = edge_table.sort_values(
        ["source_kinase", "target_kinase"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)

    return SignalomeKinaseNetwork(
        correlation_matrix=kinase_correlation_matrix,
        node_table=node_table,
        edge_table=edge_table,
        neighbor_map=neighbor_map,
    )


def subset_signalome_network(
    kinase_network: pd.DataFrame,
    kinases_of_interest: Sequence[str],
) -> pd.DataFrame:
    """Return the symmetric kinase network subset around selected kinases."""

    selected = [
        kinase for kinase in kinases_of_interest if kinase in kinase_network.index
    ]
    return kinase_network.loc[selected, selected]


def build_signalome_support_matrix(
    site_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    kinases_of_interest: Sequence[str],
) -> pd.DataFrame:
    """Build a kinase-by-protein support matrix from aligned assignments."""

    support_rows: dict[str, list[int]] = {}
    site_index = pd.Index(site_assignments.index.astype(str), name="site_id")
    if site_index.has_duplicates:
        msg = "site_assignments index must contain unique site IDs"
        raise InputCompatibilityError(msg)
    site_values = site_index.to_numpy(dtype=object, copy=False)
    site_id_set = set(site_values.tolist())

    for kinase in kinases_of_interest:
        substrates = kinase_substrates.get(kinase, ())
        aligned_site_set = {
            str(site_id) for site_id in substrates if str(site_id) in site_id_set
        }
        if len(aligned_site_set) == 0:
            msg = (
                f"No aligned phosphosites were found for kinase '{kinase}' in the "
                "signalome assignment table"
            )
            raise InputCompatibilityError(msg)

        support_rows[kinase] = (
            np.isin(site_values, list(aligned_site_set)).astype(int).tolist()
        )

    return pd.DataFrame.from_dict(support_rows, orient="index", columns=site_index)
