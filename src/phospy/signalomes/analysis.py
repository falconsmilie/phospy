from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..errors import InputCompatibilityError
from ..validation.requests import (
    SignalomeRequest,
    ValidatedSignalomeRequest,
    _build_validated_signalome_request,
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
from .clustering import cluster_sites
from .results import (
    SignalomeAssignments,
    SignalomeKinaseNetwork,
    SignalomeModules,
    SignalomeResult,
)

__all__ = ["SignalomePlan", "SignalomeRunner", "build_signalome_result"]


@dataclass(frozen=True, slots=True)
class SignalomePlan:
    """Trusted aligned inputs and options for signalome construction."""

    scoring_matrix: pd.DataFrame
    pred_mat: pd.DataFrame
    expression_matrix: pd.DataFrame
    site_to_protein: pd.Series
    kinases_of_interest: tuple[str, ...]
    kinase_network_threshold: float
    signalome_cutoff: float
    module_count: int | None
    min_kinase_module_share_percent: float


class SignalomeRunner:
    """Construct a signalome result from a validated signalome plan."""

    def execute(self, plan: SignalomePlan) -> SignalomeResult:
        scoring_matrix = plan.scoring_matrix
        pred_mat = plan.pred_mat
        expression_matrix = plan.expression_matrix

        site_clusters = cluster_sites(
            scoring_matrix=scoring_matrix,
            requested_module_count=plan.module_count,
        )
        protein_modules = derive_protein_modules(
            site_clusters=site_clusters,
            site_to_protein=plan.site_to_protein,
        )
        site_assignments = build_site_assignments(
            pred_mat=pred_mat,
            protein_modules=protein_modules,
            site_to_protein=plan.site_to_protein,
        )
        selected_kinase_substrates = select_kinase_substrates(
            pred_mat=pred_mat,
            cutoff=plan.signalome_cutoff,
        )
        kinase_network, kinase_correlation_matrix = build_kinase_network(
            scoring_matrix=scoring_matrix,
            threshold=plan.kinase_network_threshold,
        )
        signalome_modules = build_signalome_module_table(
            site_assignments=site_assignments,
            kinase_substrates=selected_kinase_substrates,
        )
        protein_assignments = build_protein_assignment_table(
            site_assignments=site_assignments,
        )
        kinase_module_relationships = build_kinase_module_relationship_table(
            module_table=signalome_modules,
        )
        network = build_kinase_network_view(
            kinase_network=kinase_network,
            kinase_correlation_matrix=kinase_correlation_matrix,
            kinase_substrates=selected_kinase_substrates,
        )
        expanded_signalomes = build_expanded_signalomes(
            kinases_of_interest=plan.kinases_of_interest,
            kinase_network=kinase_network,
            kinase_substrates=selected_kinase_substrates,
            signalome_modules=signalome_modules,
            site_assignments=site_assignments,
            expression_matrix=expression_matrix,
            min_kinase_module_share_percent=plan.min_kinase_module_share_percent,
        )

        return SignalomeResult(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            modules=SignalomeModules(
                module_table=signalome_modules,
                kinase_module_relationships=kinase_module_relationships,
            ),
            assignments=SignalomeAssignments(
                site_assignments=site_assignments,
                protein_assignments=protein_assignments,
            ),
            network=network,
            kinase_substrate_map=selected_kinase_substrates,
            expanded_signalomes=expanded_signalomes,
        )


def build_signalome_plan(
    request: ValidatedSignalomeRequest,
) -> SignalomePlan:
    """Convert a validated request into an executable signalome plan."""

    return SignalomePlan(
        scoring_matrix=request.scoring_matrix,
        pred_mat=request.pred_mat,
        expression_matrix=request.expression_matrix,
        site_to_protein=request.site_to_protein,
        kinases_of_interest=request.request.kinases_of_interest,
        kinase_network_threshold=request.request.kinase_network_threshold,
        signalome_cutoff=request.request.signalome_cutoff,
        module_count=request.request.module_count,
        min_kinase_module_share_percent=(
            request.request.min_kinase_module_share_percent
        ),
    )


def execute_validated_signalome_request(
    request: ValidatedSignalomeRequest,
) -> SignalomeResult:
    """Build a signalome result from a trusted validated request."""

    return SignalomeRunner().execute(build_signalome_plan(request))


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
) -> SignalomeResult:
    """Build a structured signalome result from validated aligned inputs."""

    request = SignalomeRequest.validate_request(
        kinases_of_interest=kinases_of_interest,
        site_to_protein=(None if site_to_protein is None else dict(site_to_protein)),
        kinase_network_threshold=kinase_network_threshold,
        signalome_cutoff=signalome_cutoff,
        module_count=module_count,
        min_kinase_module_share_percent=min_kinase_module_share_percent,
    )
    validated = _build_validated_signalome_request(
        request=request,
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        scoring_context="scoring_matrix",
        pred_mat_context="pred_mat",
        expression_context="expression_matrix",
    )

    return execute_validated_signalome_request(validated)


def build_kinase_network(
    *,
    scoring_matrix: pd.DataFrame,
    threshold: float,
) -> tuple[dict[str, tuple[str, ...]], pd.DataFrame]:
    """Build the kinase neighbor map and correlation matrix."""

    kinase_correlation_matrix = scoring_matrix.corr().fillna(0.0)
    if kinase_correlation_matrix.empty:
        msg = "scoring_matrix must contain at least one kinase column"
        raise InputCompatibilityError(msg)

    kinase_names = kinase_correlation_matrix.columns.astype(str)
    kinase_correlation_values = kinase_correlation_matrix.to_numpy(
        dtype=float, copy=True
    )
    np.fill_diagonal(kinase_correlation_values, 0.0)
    kinase_correlation_matrix = pd.DataFrame(
        kinase_correlation_values,
        index=kinase_names.copy(),
        columns=kinase_names.copy(),
    )

    adjacency_mask = kinase_correlation_values > threshold
    kinase_name_values = kinase_names.to_numpy(dtype=object, copy=False)
    kinase_network = {
        str(kinase_name): tuple(kinase_name_values[adjacency_mask[position]].tolist())
        for position, kinase_name in enumerate(kinase_name_values)
    }

    return kinase_network, kinase_correlation_matrix


def build_kinase_network_view(
    *,
    kinase_network: Mapping[str, Sequence[str]],
    kinase_correlation_matrix: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
) -> SignalomeKinaseNetwork:
    """Build the network-centric signalome view from derived network parts."""

    neighbor_sets = {
        str(kinase): frozenset(str(neighbor) for neighbor in neighbors)
        for kinase, neighbors in kinase_network.items()
    }
    node_rows = [
        {
            "kinase": kinase,
            "degree": len(neighbor_sets[kinase]),
            "n_substrates": len(tuple(kinase_substrates.get(kinase, ()))),
        }
        for kinase in sorted(neighbor_sets)
    ]
    node_table = pd.DataFrame.from_records(node_rows).set_index("kinase")
    node_table = node_table.astype({"degree": int, "n_substrates": int})
    node_table.index.name = "kinase"

    columns = kinase_correlation_matrix.columns.astype(str)
    column_positions = {
        str(column): position for position, column in enumerate(columns)
    }
    correlation_values = kinase_correlation_matrix.to_numpy(dtype=float, copy=False)
    edge_rows = [
        {
            "source_kinase": source,
            "target_kinase": target,
            "correlation": float(
                correlation_values[
                    column_positions[source],
                    column_positions[target],
                ]
            ),
        }
        for source, neighbors in neighbor_sets.items()
        for target in sorted(neighbors)
        if source in column_positions
        and target in column_positions
        and column_positions[source] < column_positions[target]
    ]

    edge_table = pd.DataFrame.from_records(edge_rows)
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
        neighbor_map={str(key): tuple(value) for key, value in kinase_network.items()},
    )
