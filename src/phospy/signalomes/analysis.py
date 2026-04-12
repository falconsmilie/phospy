from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

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
from .clustering import cluster_sites
from .results import (
    SignalomeAssignments,
    SignalomeKinaseNetwork,
    SignalomeModules,
    SignalomeResult,
)

__all__ = ["SignalomeRunner", "build_signalome_result", "execute_signalome_inputs"]


@dataclass(slots=True)
class SignalomeRunner:
    """Construct a signalome result from trusted aligned signalome inputs."""

    def execute(self, inputs: SignalomeInputs) -> SignalomeResult:
        scoring_matrix = inputs.scoring_matrix
        pred_mat = inputs.pred_mat
        expression_matrix = inputs.expression_matrix
        site_clusters = cluster_sites(
            scoring_matrix=scoring_matrix,
            requested_module_count=inputs.module_count,
        )
        protein_modules = derive_protein_modules(
            site_clusters=site_clusters,
            site_to_protein=inputs.site_to_protein,
        )
        site_assignments = build_site_assignments(
            pred_mat=pred_mat,
            protein_modules=protein_modules,
            site_to_protein=inputs.site_to_protein,
        )
        selected_kinase_substrates = select_kinase_substrates(
            pred_mat=pred_mat,
            cutoff=inputs.signalome_cutoff,
        )
        kinase_network, kinase_correlation_matrix = build_kinase_network(
            scoring_matrix=scoring_matrix,
            threshold=inputs.kinase_network_threshold,
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
            kinases_of_interest=inputs.kinases_of_interest,
            kinase_network=kinase_network,
            kinase_substrates=selected_kinase_substrates,
            signalome_modules=signalome_modules,
            site_assignments=site_assignments,
            expression_matrix=expression_matrix,
            min_kinase_module_share_percent=inputs.min_kinase_module_share_percent,
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

    kinase_order = [str(kinase) for kinase in kinase_correlation_matrix.index]
    node_rows: list[dict[str, object]] = []
    neighbor_map: dict[str, tuple[str, ...]] = {}

    for kinase in kinase_order:
        neighbors = tuple(
            str(neighbor)
            for neighbor, value in kinase_network.loc[kinase].items()
            if float(value) > 0.0
        )
        neighbor_map[kinase] = neighbors
        node_rows.append(
            {
                "kinase": kinase,
                "degree": len(neighbors),
                "n_substrates": len(tuple(kinase_substrates.get(kinase, ()))),
            }
        )

    node_table = pd.DataFrame.from_records(node_rows).set_index("kinase")
    node_table.index.name = "kinase"
    node_table = node_table.astype({"degree": int, "n_substrates": int})

    edge_rows: list[dict[str, object]] = []
    for source_position, source_kinase in enumerate(kinase_order):
        for target_kinase in kinase_order[source_position + 1 :]:
            correlation = float(kinase_network.loc[source_kinase, target_kinase])
            if correlation <= 0.0:
                continue
            edge_rows.append(
                {
                    "source_kinase": source_kinase,
                    "target_kinase": target_kinase,
                    "correlation": float(
                        kinase_correlation_matrix.loc[source_kinase, target_kinase]
                    ),
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
    site_index = site_assignments.index

    for kinase in kinases_of_interest:
        substrates = kinase_substrates.get(kinase, ())
        aligned_sites = [site for site in substrates if site in site_index]
        if len(aligned_sites) == 0:
            msg = (
                f"No aligned phosphosites were found for kinase '{kinase}' in the "
                "signalome assignment table"
            )
            raise InputCompatibilityError(msg)

        support_rows[kinase] = [
            1 if site in aligned_sites else 0 for site in site_index
        ]

    return pd.DataFrame.from_dict(support_rows, orient="index", columns=site_index)
