#!/usr/bin/env python3
"""Benchmark signalome/prediction science hot paths against legacy-style baselines.

Targets:
- `phospy.workflows.signalome.science.build_signalome_module_table`
- `phospy.workflows.signalome.science.build_expanded_signalome_table`
- `phospy.workflows.kinase.science.build_prediction_outputs`
"""

from __future__ import annotations

import json
import sys
import time
from collections.abc import Mapping, Sequence
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


def _legacy_build_signalome_module_table(
    *,
    module_assignments: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
    kinase_order: Sequence[str],
) -> pd.DataFrame:
    module_index = pd.Index(
        sorted(
            {
                int(value)
                for value in module_assignments.loc[:, "module_id"]
                if int(value) > 0
            }
        ),
        name="module_id",
    )
    kinase_index = pd.Index([str(kinase) for kinase in kinase_order], name="kinase")
    module_table = pd.DataFrame(
        0.0,
        index=module_index.copy(),
        columns=kinase_index.copy(),
    )

    protein_to_module = (
        module_assignments.loc[:, ["protein_id", "module_id"]]
        .drop_duplicates(subset=["protein_id"])
        .set_index("protein_id")
        .loc[:, "module_id"]
        .astype("int64")
    )
    protein_to_module = protein_to_module.loc[protein_to_module > 0]
    site_to_protein = module_assignments.loc[:, "protein_id"].astype(str)
    site_to_protein.index = pd.Index(
        site_to_protein.index.astype(str),
        name="site_id",
    )

    for kinase in kinase_index:
        substrate_sites = pd.Index(
            [str(site_id) for site_id in kinase_substrates.get(str(kinase), ())],
            name="site_id",
        )
        if substrate_sites.empty:
            continue
        substrate_proteins = (
            site_to_protein.reindex(substrate_sites).dropna().astype(str)
        )
        if substrate_proteins.empty:
            continue
        unique_proteins = pd.Index(sorted(set(substrate_proteins.tolist())))
        module_hits = (
            protein_to_module.reindex(unique_proteins).dropna().astype("int64")
        )
        if module_hits.empty:
            continue
        counts = module_hits.value_counts().astype(float)
        module_table.loc[counts.index.astype(int), kinase] = counts.to_numpy(
            dtype=float, copy=False
        )

    row_totals = module_table.sum(axis=1)
    non_zero_rows = row_totals > 0.0
    if non_zero_rows.any():
        module_table.loc[non_zero_rows] = (
            module_table.loc[non_zero_rows].div(row_totals.loc[non_zero_rows], axis=0)
            * 100.0
        )
    return module_table.astype(float).round(3)


def _legacy_build_prediction_outputs(
    *,
    prediction_score_matrix: pd.DataFrame,
    selected_kinases: pd.Index,
    candidate_substrates: dict[str, list[str]],
    top_k: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pred_mat = pd.DataFrame(
        np.nan,
        index=prediction_score_matrix.index.copy(),
        columns=selected_kinases.copy(),
    )
    pred_mat.index.name = prediction_score_matrix.index.name
    pred_mat.columns.name = "kinase"

    substrate_rows: list[dict[str, object]] = []
    for kinase in selected_kinases:
        candidate_sites = candidate_substrates.get(str(kinase), [])
        available_sites = [
            site for site in candidate_sites if site in prediction_score_matrix.index
        ]
        if not available_sites:
            continue
        ranked_sites = (
            prediction_score_matrix.loc[available_sites, kinase]
            .astype(float)
            .dropna()
            .sort_values(ascending=False)
            .head(top_k)
        )
        if ranked_sites.empty:
            continue
        pred_mat.loc[ranked_sites.index, kinase] = ranked_sites.values
        for rank, (site_id, score) in enumerate(ranked_sites.items(), start=1):
            substrate_rows.append(
                {
                    "kinase": str(kinase),
                    "substrate_site": site_id,
                    "score": float(score),
                    "rank": rank,
                }
            )

    substrate_list = pd.DataFrame(
        substrate_rows,
        columns=["kinase", "substrate_site", "score", "rank"],
    )
    return pred_mat, substrate_list


def _build_signalome_inputs(
    *,
    n_sites: int = 12000,
    n_kinases: int = 220,
    random_state: int = 23,
) -> tuple[pd.DataFrame, dict[str, tuple[str, ...]], list[str]]:
    rng = np.random.default_rng(random_state)
    site_ids = [f"SITE_{index + 1}" for index in range(n_sites)]
    kinase_order = [f"KINASE_{index + 1}" for index in range(n_kinases)]
    protein_ids = [
        f"P_{(index % max(n_sites // 2, 1)) + 1}" for index in range(n_sites)
    ]
    module_ids = rng.integers(0, 36, size=n_sites).astype(int)
    module_assignments = pd.DataFrame(
        {
            "protein_id": protein_ids,
            "module_id": module_ids,
            "top_kinase": rng.choice(
                np.asarray(kinase_order, dtype=object), size=n_sites
            ).astype(str),
            "top_score": rng.uniform(0.55, 0.99, size=n_sites).astype(float),
        },
        index=pd.Index(site_ids, name="site_id"),
    )
    module_assignments.loc[:, "top_kinase_weights"] = [
        ((str(top_kinase), 1.0),)
        for top_kinase in module_assignments.loc[:, "top_kinase"].tolist()
    ]
    site_array = np.asarray(site_ids, dtype=object)
    kinase_substrates: dict[str, tuple[str, ...]] = {}
    for kinase in kinase_order:
        substrate_count = int(rng.integers(160, 520))
        sampled_sites = site_array[
            rng.choice(len(site_array), size=substrate_count, replace=False)
        ]
        kinase_substrates[kinase] = tuple(str(site_id) for site_id in sampled_sites)
    return module_assignments, kinase_substrates, kinase_order


def _build_expanded_inputs(
    *,
    module_assignments: pd.DataFrame,
    kinase_order: Sequence[str],
    random_state: int = 31,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(random_state)
    module_ids = sorted(
        {
            int(module_id)
            for module_id in module_assignments.loc[:, "module_id"].tolist()
            if int(module_id) > 0
        }
    )
    signalome_values = rng.uniform(
        0.0,
        100.0,
        size=(len(module_ids), len(kinase_order)),
    )
    sparse_mask = rng.random(size=signalome_values.shape) < 0.57
    signalome_values[sparse_mask] = 0.0
    signalome_modules = pd.DataFrame(
        signalome_values,
        index=pd.Index(module_ids, name="module_id"),
        columns=pd.Index([str(kinase) for kinase in kinase_order], name="kinase"),
        dtype=float,
    )

    unique_edges: set[tuple[str, str]] = set()
    kinase_array = np.asarray([str(kinase) for kinase in kinase_order], dtype=object)
    if kinase_array.size > 1:
        for idx in range(int(kinase_array.size - 1)):
            left = str(kinase_array[idx])
            right = str(kinase_array[idx + 1])
            unique_edges.add(tuple(sorted((left, right))))
        target_edge_count = int(max(len(kinase_order) * 4, len(kinase_order) + 1))
        while len(unique_edges) < target_edge_count:
            positions = rng.choice(kinase_array.size, size=2, replace=False)
            left = str(kinase_array[int(positions[0])])
            right = str(kinase_array[int(positions[1])])
            unique_edges.add(tuple(sorted((left, right))))
    edge_rows = [
        {
            "source_kinase": source,
            "target_kinase": target,
            "correlation": float(rng.uniform(0.6, 0.99)),
        }
        for source, target in sorted(unique_edges)
    ]
    kinase_network_edges = pd.DataFrame.from_records(
        edge_rows,
        columns=["source_kinase", "target_kinase", "correlation"],
    ).astype({"source_kinase": str, "target_kinase": str, "correlation": float})
    return signalome_modules, kinase_network_edges


def _legacy_build_expanded_signalome_table(
    *,
    module_assignments: pd.DataFrame,
    signalome_modules: pd.DataFrame,
    kinase_network_edges: pd.DataFrame,
    kinase_substrates: Mapping[str, Sequence[str]],
) -> pd.DataFrame:
    site_index = pd.Index(module_assignments.index.astype(str), name="site_id")
    indexed_assignments = module_assignments.copy(deep=False)
    indexed_assignments.index = site_index

    module_id_values = indexed_assignments.loc[:, "module_id"].astype("int64")
    protein_ids = indexed_assignments.loc[:, "protein_id"].astype(str)
    top_kinases = indexed_assignments.loc[:, "top_kinase"].astype(str)
    top_scores = indexed_assignments.loc[:, "top_score"].astype(float)

    kinase_order = [str(kinase) for kinase in signalome_modules.columns.astype(str)]
    neighbor_map: dict[str, set[str]] = {kinase: set() for kinase in kinase_order}
    for row in kinase_network_edges.loc[
        :, ["source_kinase", "target_kinase"]
    ].itertuples(index=False):
        source = str(row.source_kinase)
        target = str(row.target_kinase)
        neighbor_map.setdefault(source, set()).add(target)
        neighbor_map.setdefault(target, set()).add(source)

    site_positions = {
        str(site_id): int(position)
        for position, site_id in enumerate(
            site_index.to_numpy(dtype=object, copy=False)
        )
    }
    support_by_kinase: dict[str, np.ndarray] = {}
    for kinase, substrates in kinase_substrates.items():
        weights = np.zeros(int(site_index.size), dtype=float)
        for site_id in substrates:
            position = site_positions.get(str(site_id))
            if position is None:
                continue
            weights[position] = 1.0
        support_by_kinase[str(kinase)] = weights

    site_positions_array = np.arange(site_index.size, dtype=np.int64)
    site_module_ids = module_id_values.to_numpy(dtype=np.int64, copy=False)
    site_proteins = protein_ids.to_numpy(dtype=object, copy=False)
    site_top_kinases = top_kinases.to_numpy(dtype=object, copy=False)
    site_top_scores = top_scores.to_numpy(dtype=float, copy=False)
    site_ids = site_index.to_numpy(dtype=object, copy=False)

    rows: list[dict[str, object]] = []
    for focal_kinase in kinase_order:
        linked_kinases = tuple(
            dict.fromkeys(
                (focal_kinase, *tuple(sorted(neighbor_map.get(focal_kinase, set()))))
            )
        )
        regulated_module_ids = tuple(
            int(module_id)
            for module_id, share in signalome_modules.loc[:, focal_kinase].items()
            if float(share) > 1.0
        )
        regulated_module_set = set(regulated_module_ids)

        linked_kinases_json = json.dumps(
            list(linked_kinases),
            separators=(",", ":"),
            ensure_ascii=True,
        )
        regulated_module_ids_json = json.dumps(
            list(regulated_module_ids),
            separators=(",", ":"),
            ensure_ascii=True,
        )

        matched_site_count = 0
        for position, site_id, module_id in zip(
            site_positions_array,
            site_ids,
            site_module_ids,
            strict=True,
        ):
            if int(module_id) not in regulated_module_set:
                continue
            support_kinases: list[str] = []
            support_weight = 0.0
            for linked_kinase in linked_kinases:
                kinase_support = support_by_kinase.get(linked_kinase)
                if kinase_support is None:
                    continue
                weight = float(kinase_support[int(position)])
                if weight <= 0.0:
                    continue
                support_kinases.append(linked_kinase)
                support_weight += weight
            if support_weight <= 0.0:
                continue
            matched_site_count += 1
            rows.append(
                {
                    "kinase": focal_kinase,
                    "row_kind": "site",
                    "assignment_policy": "cutoff_binary",
                    "linked_kinases": linked_kinases_json,
                    "regulated_module_ids": regulated_module_ids_json,
                    "site_id": str(site_id),
                    "site_order": int(position),
                    "protein_id": str(site_proteins[int(position)]),
                    "module_id": int(module_id),
                    "support_kinases": json.dumps(
                        support_kinases,
                        separators=(",", ":"),
                        ensure_ascii=True,
                    ),
                    "support_weight": float(support_weight),
                    "top_kinase": str(site_top_kinases[int(position)]),
                    "top_score": float(site_top_scores[int(position)]),
                }
            )
        if matched_site_count == 0:
            rows.append(
                {
                    "kinase": focal_kinase,
                    "row_kind": "summary",
                    "assignment_policy": "cutoff_binary",
                    "linked_kinases": linked_kinases_json,
                    "regulated_module_ids": regulated_module_ids_json,
                    "site_id": "",
                    "site_order": -1,
                    "protein_id": "",
                    "module_id": 0,
                    "support_kinases": "[]",
                    "support_weight": 0.0,
                    "top_kinase": "",
                    "top_score": np.nan,
                }
            )

    expanded = pd.DataFrame.from_records(rows).astype(
        {
            "kinase": str,
            "row_kind": str,
            "assignment_policy": str,
            "linked_kinases": str,
            "regulated_module_ids": str,
            "site_id": str,
            "site_order": "int64",
            "protein_id": str,
            "module_id": "int64",
            "support_kinases": str,
            "support_weight": float,
            "top_kinase": str,
            "top_score": float,
        }
    )
    return expanded.sort_values(
        ["kinase", "row_kind", "site_order", "site_id"],
        ascending=[True, True, True, True],
        kind="stable",
    ).reset_index(drop=True)


def _build_prediction_inputs(
    *,
    n_sites: int = 16000,
    n_kinases: int = 180,
    candidate_sites_per_kinase: int = 900,
    top_k: int = 50,
    random_state: int = 29,
) -> tuple[pd.DataFrame, pd.Index, dict[str, list[str]], int]:
    rng = np.random.default_rng(random_state)
    site_ids = pd.Index(
        [f"SITE_{index + 1}" for index in range(n_sites)], name="site_id"
    )
    kinase_ids = pd.Index(
        [f"KINASE_{index + 1}" for index in range(n_kinases)], name="kinase"
    )
    score_values = rng.uniform(0.0, 1.0, size=(n_sites, n_kinases))
    nan_mask = rng.random(size=(n_sites, n_kinases)) < 0.08
    score_values[nan_mask] = np.nan
    prediction_score_matrix = pd.DataFrame(
        score_values,
        index=site_ids.copy(),
        columns=kinase_ids.copy(),
    )
    site_array = site_ids.to_numpy(dtype=object, copy=False)
    candidate_substrates: dict[str, list[str]] = {}
    for kinase in kinase_ids:
        positions = rng.choice(
            n_sites,
            size=min(candidate_sites_per_kinase, n_sites),
            replace=False,
        )
        candidate_substrates[str(kinase)] = [str(site_array[pos]) for pos in positions]
    return prediction_score_matrix, kinase_ids.copy(), candidate_substrates, top_k


def main() -> None:
    from phospy.workflows.kinase.science import build_prediction_outputs
    from phospy.workflows.signalome.science import (
        build_expanded_signalome_table,
        build_signalome_module_table,
    )

    module_assignments, kinase_substrates, kinase_order = _build_signalome_inputs()
    signalome_modules, kinase_network_edges = _build_expanded_inputs(
        module_assignments=module_assignments,
        kinase_order=kinase_order,
    )
    prediction_score_matrix, selected_kinases, candidate_substrates, top_k = (
        _build_prediction_inputs()
    )

    baseline_modules = _legacy_build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates=kinase_substrates,
        kinase_order=kinase_order,
    )
    optimized_modules = build_signalome_module_table(
        module_assignments=module_assignments,
        kinase_substrates=kinase_substrates,
        kinase_order=kinase_order,
    )
    pd.testing.assert_frame_equal(
        baseline_modules, optimized_modules, check_dtype=False
    )
    baseline_expanded = _legacy_build_expanded_signalome_table(
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates=kinase_substrates,
    )
    optimized_expanded = build_expanded_signalome_table(
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates=kinase_substrates,
    )
    pd.testing.assert_frame_equal(
        baseline_expanded,
        optimized_expanded,
        check_dtype=False,
    )

    baseline_pred_mat, baseline_substrate_list = _legacy_build_prediction_outputs(
        prediction_score_matrix=prediction_score_matrix,
        selected_kinases=selected_kinases,
        candidate_substrates=candidate_substrates,
        top_k=top_k,
    )
    optimized_pred_mat, optimized_substrate_list = build_prediction_outputs(
        prediction_score_matrix=prediction_score_matrix,
        selected_kinases=selected_kinases,
        candidate_substrates=candidate_substrates,
        top_k=top_k,
    )
    pd.testing.assert_frame_equal(
        baseline_pred_mat,
        optimized_pred_mat,
        check_dtype=False,
    )
    pd.testing.assert_frame_equal(
        baseline_substrate_list.reset_index(drop=True),
        optimized_substrate_list.reset_index(drop=True),
        check_dtype=False,
    )

    repeats = 4
    signalome_legacy_seconds = _time_call(
        repeats,
        _legacy_build_signalome_module_table,
        module_assignments=module_assignments,
        kinase_substrates=kinase_substrates,
        kinase_order=kinase_order,
    )
    signalome_optimized_seconds = _time_call(
        repeats,
        build_signalome_module_table,
        module_assignments=module_assignments,
        kinase_substrates=kinase_substrates,
        kinase_order=kinase_order,
    )
    expanded_legacy_seconds = _time_call(
        repeats,
        _legacy_build_expanded_signalome_table,
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates=kinase_substrates,
    )
    expanded_optimized_seconds = _time_call(
        repeats,
        build_expanded_signalome_table,
        module_assignments=module_assignments,
        signalome_modules=signalome_modules,
        kinase_network_edges=kinase_network_edges,
        kinase_substrates=kinase_substrates,
    )
    prediction_legacy_seconds = _time_call(
        repeats,
        _legacy_build_prediction_outputs,
        prediction_score_matrix=prediction_score_matrix,
        selected_kinases=selected_kinases,
        candidate_substrates=candidate_substrates,
        top_k=top_k,
    )
    prediction_optimized_seconds = _time_call(
        repeats,
        build_prediction_outputs,
        prediction_score_matrix=prediction_score_matrix,
        selected_kinases=selected_kinases,
        candidate_substrates=candidate_substrates,
        top_k=top_k,
    )

    print(f"repeats={repeats}")
    print(f"signalome_legacy_mean_seconds={signalome_legacy_seconds:.6f}")
    print(f"signalome_optimized_mean_seconds={signalome_optimized_seconds:.6f}")
    print(
        "signalome_speedup="
        f"{(signalome_legacy_seconds / signalome_optimized_seconds):.3f}x"
    )
    print(f"expanded_legacy_mean_seconds={expanded_legacy_seconds:.6f}")
    print(f"expanded_optimized_mean_seconds={expanded_optimized_seconds:.6f}")
    print(
        "expanded_speedup="
        f"{(expanded_legacy_seconds / expanded_optimized_seconds):.3f}x"
    )
    print(f"prediction_legacy_mean_seconds={prediction_legacy_seconds:.6f}")
    print(f"prediction_optimized_mean_seconds={prediction_optimized_seconds:.6f}")
    print(
        "prediction_speedup="
        f"{(prediction_legacy_seconds / prediction_optimized_seconds):.3f}x"
    )


if __name__ == "__main__":
    main()
