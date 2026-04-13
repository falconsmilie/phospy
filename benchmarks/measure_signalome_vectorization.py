#!/usr/bin/env python3
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _time_call(repeats: int, func, *args, **kwargs) -> float:
    durations: list[float] = []
    for _ in range(repeats):
        start = time.perf_counter()
        func(*args, **kwargs)
        durations.append(time.perf_counter() - start)
    return sum(durations) / len(durations)


def _build_signalome_inputs(
    *,
    n_sites: int = 1200,
    n_kinases: int = 48,
    n_kinases_of_interest: int = 8,
    random_state: int = 17,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    dict[str, tuple[str, ...]],
    pd.DataFrame,
    pd.DataFrame,
    tuple[str, ...],
]:
    rng = np.random.default_rng(random_state)
    site_ids = pd.Index(
        [f"PROTEIN_{i};S{i};" for i in range(1, n_sites + 1)], dtype=object
    )
    kinase_names = pd.Index(
        [f"KINASE_{i}" for i in range(1, n_kinases + 1)],
        dtype=object,
    )

    raw_correlation = rng.uniform(0.0, 1.0, size=(n_kinases, n_kinases))
    correlation_values = (raw_correlation + raw_correlation.T) / 2.0
    np.fill_diagonal(correlation_values, 1.0)
    kinase_correlation_matrix = pd.DataFrame(
        correlation_values,
        index=kinase_names.copy(),
        columns=kinase_names.copy(),
    )
    kinase_network = kinase_correlation_matrix.where(
        kinase_correlation_matrix >= 0.85, 0.0
    )
    np.fill_diagonal(kinase_network.to_numpy(dtype=float, copy=False), 0.0)

    kinase_substrates = {
        str(kinase): tuple(
            site_ids[
                rng.choice(n_sites, size=max(n_sites // 12, 1), replace=False)
            ].tolist()
        )
        for kinase in kinase_names
    }
    signalome_modules = pd.DataFrame(
        rng.uniform(0.0, 100.0, size=(12, n_kinases)),
        index=pd.Index(range(1, 13), name="module_id"),
        columns=kinase_names.copy(),
    ).round(3)
    site_assignments = pd.DataFrame(
        {
            "protein_id": [f"PROTEIN_{i}" for i in range(1, n_sites + 1)],
            "module_id": rng.integers(1, 13, size=n_sites),
            "top_kinase": rng.choice(
                kinase_names.to_numpy(dtype=object, copy=False), size=n_sites
            ),
            "top_kinase_candidates": ["[]"] * n_sites,
            "top_kinase_tie_count": np.ones(n_sites, dtype=int),
            "top_kinase_is_ambiguous": np.zeros(n_sites, dtype=bool),
            "top_score": rng.uniform(0.5, 1.0, size=n_sites),
        },
        index=pd.Index(site_ids, name="site_id"),
    )
    expression_matrix = pd.DataFrame(
        rng.normal(size=(n_sites, 8)),
        index=site_ids.copy(),
        columns=[f"sample_{i}" for i in range(1, 9)],
    )
    kinases_of_interest = tuple(
        str(kinase) for kinase in kinase_names[:n_kinases_of_interest]
    )
    return (
        kinase_network,
        kinase_correlation_matrix,
        kinase_substrates,
        signalome_modules,
        site_assignments,
        expression_matrix,
        kinases_of_interest,
    )


def _baseline_build_kinase_network_view(
    *,
    kinase_network: pd.DataFrame,
    kinase_correlation_matrix: pd.DataFrame,
    kinase_substrates: dict[str, tuple[str, ...]],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, tuple[str, ...]]]:
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
    return node_table, edge_table, neighbor_map


def _baseline_build_expanded_signalomes(
    *,
    kinases_of_interest: tuple[str, ...],
    kinase_network: dict[str, tuple[str, ...]],
    kinase_substrates: dict[str, tuple[str, ...]],
    signalome_modules: pd.DataFrame,
    site_assignments: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    min_kinase_module_share_percent: float,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame]]:
    expanded: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    available_sites = pd.Index(site_assignments.index.astype(str), dtype=object)

    for kinase in kinases_of_interest:
        linked_kinases = tuple(dict.fromkeys((kinase, *kinase_network.get(kinase, ()))))
        regulated_module_ids = tuple(
            int(module_id)
            for module_id in signalome_modules.index[
                signalome_modules.loc[:, kinase] > min_kinase_module_share_percent
            ].tolist()
        )
        substrate_site_ids = tuple(
            site_id
            for linked_kinase in linked_kinases
            for site_id in kinase_substrates.get(linked_kinase, ())
        )
        substrate_site_index = available_sites.intersection(
            pd.Index(substrate_site_ids, dtype=object)
        )
        site_mask = site_assignments.loc[substrate_site_index, "module_id"].isin(
            regulated_module_ids
        )
        selected_site_ids = site_mask.index[site_mask]
        expanded[kinase] = (
            expression_matrix.loc[selected_site_ids].copy(deep=True),
            site_assignments.loc[selected_site_ids].copy(deep=True),
        )
    return expanded


def main() -> None:
    from phospy.signalomes.analysis import build_kinase_network_view
    from phospy.signalomes.assignments import build_expanded_signalomes

    (
        kinase_network,
        kinase_correlation_matrix,
        kinase_substrates,
        signalome_modules,
        site_assignments,
        expression_matrix,
        kinases_of_interest,
    ) = _build_signalome_inputs()
    network_view = build_kinase_network_view(
        kinase_network=kinase_network,
        kinase_correlation_matrix=kinase_correlation_matrix,
        kinase_substrates=kinase_substrates,
    )
    repeats = 3

    baseline_network_runtime = _time_call(
        repeats,
        _baseline_build_kinase_network_view,
        kinase_network=kinase_network,
        kinase_correlation_matrix=kinase_correlation_matrix,
        kinase_substrates=kinase_substrates,
    )
    optimized_network_runtime = _time_call(
        repeats,
        build_kinase_network_view,
        kinase_network=kinase_network,
        kinase_correlation_matrix=kinase_correlation_matrix,
        kinase_substrates=kinase_substrates,
    )
    baseline_expansion_runtime = _time_call(
        repeats,
        _baseline_build_expanded_signalomes,
        kinases_of_interest=kinases_of_interest,
        kinase_network=network_view.neighbor_map,
        kinase_substrates=kinase_substrates,
        signalome_modules=signalome_modules,
        site_assignments=site_assignments,
        expression_matrix=expression_matrix,
        min_kinase_module_share_percent=25.0,
    )
    optimized_expansion_runtime = _time_call(
        repeats,
        build_expanded_signalomes,
        kinases_of_interest=kinases_of_interest,
        kinase_network=network_view.neighbor_map,
        kinase_substrates=kinase_substrates,
        signalome_modules=signalome_modules,
        site_assignments=site_assignments,
        expression_matrix=expression_matrix,
        min_kinase_module_share_percent=25.0,
    )

    print(f"baseline_network_view_mean_seconds={baseline_network_runtime:.6f}")
    print(f"optimized_network_view_mean_seconds={optimized_network_runtime:.6f}")
    print(
        "network_view_speedup="
        f"{baseline_network_runtime / optimized_network_runtime:.3f}x"
    )
    print(f"baseline_expansion_mean_seconds={baseline_expansion_runtime:.6f}")
    print(f"optimized_expansion_mean_seconds={optimized_expansion_runtime:.6f}")
    print(
        "expansion_speedup="
        f"{baseline_expansion_runtime / optimized_expansion_runtime:.3f}x"
    )


if __name__ == "__main__":
    main()
