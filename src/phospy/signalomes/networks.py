from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from .results import SignalomeResult

__all__ = [
    "SignalomeNetworkData",
    "SignalomeNetworkEdge",
    "SignalomeNetworkNode",
    "build_signalome_network_data",
]


@dataclass(frozen=True, slots=True)
class SignalomeNetworkNode:
    """Explicit node model for derived signalome kinase-network outputs."""

    kinase: str
    degree: int
    n_substrates: int
    module_count: int
    total_share_percent: float
    is_kinase_of_interest: bool


@dataclass(frozen=True, slots=True)
class SignalomeNetworkEdge:
    """Explicit edge model for derived signalome kinase-network outputs."""

    source_kinase: str
    target_kinase: str
    correlation: float
    shared_module_count: int
    shared_modules: tuple[int, ...]
    source_is_kinase_of_interest: bool
    target_is_kinase_of_interest: bool


@dataclass(frozen=True, slots=True)
class SignalomeNetworkData:
    """Graph-friendly kinase-network data derived from a signalome result.

    This model is intentionally render-free. It exposes deterministic node,
    edge, and adjacency outputs that graph or plotting layers can consume.
    """

    adjacency_matrix: pd.DataFrame
    node_table: pd.DataFrame
    edge_table: pd.DataFrame
    node_models: tuple[SignalomeNetworkNode, ...]
    edge_models: tuple[SignalomeNetworkEdge, ...]

    def adjacency(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical kinase adjacency matrix."""

        if copy:
            return self.adjacency_matrix.copy(deep=True)
        return self.adjacency_matrix

    def nodes(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical kinase node table."""

        if copy:
            return self.node_table.copy(deep=True)
        return self.node_table

    def edges(self, *, copy: bool = True) -> pd.DataFrame:
        """Return the canonical kinase edge list."""

        if copy:
            return self.edge_table.copy(deep=True)
        return self.edge_table

    def to_frames(self, *, copy: bool = True) -> dict[str, pd.DataFrame]:
        """Return the named graph-friendly network tables."""

        return {
            "signalome_network_nodes": self.nodes(copy=copy),
            "signalome_network_edges": self.edges(copy=copy),
            "signalome_network_adjacency": self.adjacency(copy=copy),
        }

    def to_csv(self, directory: str | Path) -> dict[str, Path]:
        """Write the canonical network tables to CSV files."""

        target_dir = Path(directory)
        target_dir.mkdir(parents=True, exist_ok=True)

        written_paths: dict[str, Path] = {}
        for name, frame in self.to_frames(copy=True).items():
            path = target_dir / f"{name}.csv"
            write_index = not isinstance(frame.index, pd.RangeIndex)
            frame.to_csv(
                path,
                encoding="utf-8",
                float_format="%.17g",
                lineterminator="\n",
                index=write_index,
            )
            written_paths[name] = path
        return written_paths


def build_signalome_network_data(
    signalome_result: SignalomeResult,
) -> SignalomeNetworkData:
    """Build graph-friendly kinase-network data from a canonical signalome result."""

    adjacency_matrix = signalome_result.network.adjacency(copy=True)
    node_table, node_models = _build_node_outputs(signalome_result)
    edge_table, edge_models = _build_edge_outputs(signalome_result)
    return SignalomeNetworkData(
        adjacency_matrix=adjacency_matrix,
        node_table=node_table,
        edge_table=edge_table,
        node_models=node_models,
        edge_models=edge_models,
    )


def _build_node_outputs(
    signalome_result: SignalomeResult,
) -> tuple[pd.DataFrame, tuple[SignalomeNetworkNode, ...]]:
    network_nodes = signalome_result.network.nodes(copy=True)
    relationships = signalome_result.modules.to_relationship_table(copy=True)
    kinases_of_interest = set(signalome_result.kinases_of_interest)
    kinase_order = [
        str(kinase) for kinase in signalome_result.signalome_modules.columns
    ]

    module_counts = (
        relationships.groupby("kinase").size().astype(int)
        if not relationships.empty
        else pd.Series(dtype=int)
    )
    total_share_percent = (
        relationships.groupby("kinase")["share_percent"].sum().astype(float)
        if not relationships.empty
        else pd.Series(dtype=float)
    )

    rows: list[dict[str, object]] = []
    models: list[SignalomeNetworkNode] = []
    for kinase in kinase_order:
        node = SignalomeNetworkNode(
            kinase=kinase,
            degree=int(network_nodes.loc[kinase, "degree"]),
            n_substrates=int(network_nodes.loc[kinase, "n_substrates"]),
            module_count=int(module_counts.get(kinase, 0)),
            total_share_percent=float(total_share_percent.get(kinase, 0.0)),
            is_kinase_of_interest=kinase in kinases_of_interest,
        )
        models.append(node)
        rows.append(
            {
                "kinase": node.kinase,
                "degree": node.degree,
                "n_substrates": node.n_substrates,
                "module_count": node.module_count,
                "total_share_percent": node.total_share_percent,
                "is_kinase_of_interest": node.is_kinase_of_interest,
            }
        )

    node_table = pd.DataFrame.from_records(rows).set_index("kinase")
    node_table.index.name = "kinase"
    node_table = node_table.astype(
        {
            "degree": int,
            "n_substrates": int,
            "module_count": int,
            "total_share_percent": float,
            "is_kinase_of_interest": bool,
        }
    )
    return node_table, tuple(models)


def _build_edge_outputs(
    signalome_result: SignalomeResult,
) -> tuple[pd.DataFrame, tuple[SignalomeNetworkEdge, ...]]:
    network_edges = signalome_result.network.edges(copy=True)
    relationships = signalome_result.modules.to_relationship_table(copy=True)
    kinases_of_interest = set(signalome_result.kinases_of_interest)

    if relationships.empty:
        relationship_map: dict[str, tuple[int, ...]] = {}
    else:
        relationship_map = {
            str(kinase): tuple(sorted(group["module_id"].astype(int).tolist()))
            for kinase, group in relationships.groupby("kinase", sort=True)
        }

    models: list[SignalomeNetworkEdge] = []
    rows: list[dict[str, object]] = []
    for row in network_edges.itertuples(index=False):
        source_kinase = str(row.source_kinase)
        target_kinase = str(row.target_kinase)
        source_modules = set(relationship_map.get(source_kinase, ()))
        target_modules = set(relationship_map.get(target_kinase, ()))
        shared_modules = tuple(sorted(source_modules.intersection(target_modules)))
        edge = SignalomeNetworkEdge(
            source_kinase=source_kinase,
            target_kinase=target_kinase,
            correlation=float(row.correlation),
            shared_module_count=len(shared_modules),
            shared_modules=shared_modules,
            source_is_kinase_of_interest=source_kinase in kinases_of_interest,
            target_is_kinase_of_interest=target_kinase in kinases_of_interest,
        )
        models.append(edge)
        rows.append(
            {
                "source_kinase": edge.source_kinase,
                "target_kinase": edge.target_kinase,
                "correlation": edge.correlation,
                "shared_module_count": edge.shared_module_count,
                "shared_modules": "["
                + ",".join(str(module_id) for module_id in edge.shared_modules)
                + "]",
                "source_is_kinase_of_interest": edge.source_is_kinase_of_interest,
                "target_is_kinase_of_interest": edge.target_is_kinase_of_interest,
            }
        )

    edge_table = pd.DataFrame.from_records(rows)
    if edge_table.empty:
        edge_table = pd.DataFrame(
            columns=[
                "source_kinase",
                "target_kinase",
                "correlation",
                "shared_module_count",
                "shared_modules",
                "source_is_kinase_of_interest",
                "target_is_kinase_of_interest",
            ]
        )
    edge_table = edge_table.astype(
        {
            "source_kinase": str,
            "target_kinase": str,
            "correlation": float,
            "shared_module_count": int,
            "shared_modules": str,
            "source_is_kinase_of_interest": bool,
            "target_is_kinase_of_interest": bool,
        }
    )
    edge_table = edge_table.sort_values(
        ["source_kinase", "target_kinase"],
        ascending=[True, True],
        kind="stable",
    ).reset_index(drop=True)
    return edge_table, tuple(models)
