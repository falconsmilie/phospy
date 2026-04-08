from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .signalome_assignments import (
    build_expanded_signalomes,
    build_kinase_module_relationship_table,
    build_protein_assignment_table,
    build_signalome_module_table,
    build_site_assignments,
    derive_protein_modules,
    resolve_site_to_protein,
    select_kinase_substrates,
)
from .signalome_clustering import cluster_sites
from .signalome_models import (
    SignalomeAssignments,
    SignalomeKinaseNetwork,
    SignalomeModules,
    SignalomeResult,
)
from .validation.errors import InputCompatibilityError

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
    """Build a structured signalome result from trusted aligned inputs."""

    plan = SignalomePlan(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        site_to_protein=resolve_site_to_protein(
            site_ids=scoring_matrix.index.astype(str).tolist(),
            site_to_protein=site_to_protein,
        ),
        kinases_of_interest=tuple(kinases_of_interest),
        kinase_network_threshold=kinase_network_threshold,
        signalome_cutoff=signalome_cutoff,
        module_count=module_count,
        min_kinase_module_share_percent=min_kinase_module_share_percent,
    )
    return SignalomeRunner().execute(plan)


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

    kinase_correlation_values = kinase_correlation_matrix.to_numpy(copy=True)
    np.fill_diagonal(kinase_correlation_values, 0.0)
    kinase_correlation_matrix = pd.DataFrame(
        kinase_correlation_values,
        index=kinase_correlation_matrix.index.copy(),
        columns=kinase_correlation_matrix.columns.copy(),
    )

    kinase_network: dict[str, tuple[str, ...]] = {}
    for kinase in kinase_correlation_matrix.columns.astype(str):
        correlated = kinase_correlation_matrix.index[
            kinase_correlation_matrix.loc[:, kinase] > threshold
        ].astype(str)
        kinase_network[kinase] = tuple(correlated.tolist())

    return kinase_network, kinase_correlation_matrix


def build_kinase_network_view(
    *,
    kinase_network: Mapping[str, Sequence[str]],
    kinase_correlation_matrix: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
) -> SignalomeKinaseNetwork:
    """Build the network-centric signalome view from derived network parts."""

    node_rows = [
        {
            "kinase": str(kinase),
            "degree": len(tuple(neighbors)),
            "n_substrates": len(tuple(kinase_substrates.get(str(kinase), ()))),
        }
        for kinase, neighbors in sorted(kinase_network.items())
    ]
    node_table = pd.DataFrame.from_records(node_rows).set_index("kinase")
    node_table = node_table.astype({"degree": int, "n_substrates": int})
    node_table.index.name = "kinase"

    edge_rows: list[dict[str, object]] = []
    columns = kinase_correlation_matrix.columns.astype(str)
    for left_position, source in enumerate(columns):
        for target in columns[left_position + 1 :]:
            target_name = str(target)
            if target_name not in set(kinase_network.get(source, ())):
                continue
            correlation = float(kinase_correlation_matrix.loc[source, target_name])
            edge_rows.append(
                {
                    "source_kinase": source,
                    "target_kinase": target_name,
                    "correlation": correlation,
                }
            )

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
