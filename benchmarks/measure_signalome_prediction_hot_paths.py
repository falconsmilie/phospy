#!/usr/bin/env python3
from __future__ import annotations

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
    protein_ids = [
        f"P_{(index % max(n_sites // 2, 1)) + 1}" for index in range(n_sites)
    ]
    module_ids = rng.integers(0, 36, size=n_sites).astype(int)
    module_assignments = pd.DataFrame(
        {
            "protein_id": protein_ids,
            "module_id": module_ids,
        },
        index=pd.Index(site_ids, name="site_id"),
    )
    kinase_order = [f"KINASE_{index + 1}" for index in range(n_kinases)]
    site_array = np.asarray(site_ids, dtype=object)
    kinase_substrates: dict[str, tuple[str, ...]] = {}
    for kinase in kinase_order:
        substrate_count = int(rng.integers(160, 520))
        sampled_sites = site_array[
            rng.choice(len(site_array), size=substrate_count, replace=False)
        ]
        kinase_substrates[kinase] = tuple(str(site_id) for site_id in sampled_sites)
    return module_assignments, kinase_substrates, kinase_order


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
    from phospy.workflows.signalome.science import build_signalome_module_table

    module_assignments, kinase_substrates, kinase_order = _build_signalome_inputs()
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
    print(f"prediction_legacy_mean_seconds={prediction_legacy_seconds:.6f}")
    print(f"prediction_optimized_mean_seconds={prediction_optimized_seconds:.6f}")
    print(
        "prediction_speedup="
        f"{(prediction_legacy_seconds / prediction_optimized_seconds):.3f}x"
    )


if __name__ == "__main__":
    main()
