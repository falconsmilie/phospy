#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy.api import PredMatWorkflow, SignalomeWorkflow
from phospy.signalomes import (
    SignalomeMapData,
    SignalomeNetworkData,
    SignalomeResult,
)


def load_demo_inputs(
    data_dir: Path,
) -> tuple[pd.DataFrame, dict[str, list[str]], dict[str, str], dict[str, list[str]]]:
    phospho_matrix = pd.read_csv(
        data_dir / "predmat_phospho_matrix.csv",
        index_col=0,
    )
    phospho_matrix.index = phospho_matrix.index.map(str)
    substrate_map = json.loads(
        (data_dir / "predmat_substrate_map.json").read_text(encoding="utf-8")
    )
    site_sequences = json.loads(
        (data_dir / "predmat_site_sequences.json").read_text(encoding="utf-8")
    )
    motif_sequences = json.loads(
        (data_dir / "predmat_motif_sequences.json").read_text(encoding="utf-8")
    )
    return phospho_matrix, substrate_map, site_sequences, motif_sequences


def build_site_to_protein(site_ids: pd.Index) -> dict[str, str]:
    return {str(site_id): str(site_id).split(";", 1)[0] for site_id in site_ids}


def run_demo(
    outdir: Path, *, svm_mode: str = "default"
) -> tuple[
    SignalomeResult,
    SignalomeMapData,
    SignalomeNetworkData,
    dict[str, dict[str, Path]],
]:
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "examples" / "data"
    phospho_matrix, substrate_map, site_sequences, motif_sequences = load_demo_inputs(
        data_dir
    )
    outdir.mkdir(parents=True, exist_ok=True)

    pred_mat_result = PredMatWorkflow(flank_size=2, svm_mode=svm_mode).run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        min_substrates=2,
        min_motif_size=2,
        ensemble_size=3,
        top=4,
        score_threshold=0.75,
        inclusion=3,
        n_iterations=2,
        random_state=17,
    )
    signalome_result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A", "KINASE_B"],
        site_to_protein=build_site_to_protein(phospho_matrix.index),
        signalome_cutoff=0.5,
        # For larger real-world datasets, set module_count explicitly to skip
        # the extra automatic module-selection scoring pass.
    )
    map_data = signalome_result.to_map_data()
    network_data = signalome_result.to_network_data()

    written = {
        "signalome": signalome_result.to_csv(outdir / "signalome"),
        "map": map_data.to_csv(outdir / "signalome_map"),
        "network": network_data.to_csv(outdir / "signalome_network"),
    }
    return signalome_result, map_data, network_data, written


def main() -> None:
    with TemporaryDirectory(prefix="phospy-signalome-workflow-") as tmp_dir:
        signalome_result, map_data, network_data, written = run_demo(
            Path(tmp_dir), svm_mode="default"
        )
        print(f"Wrote signalome workflow tables to {tmp_dir} using svm_mode=default")
        print("Signalome modules")
        print(signalome_result.modules.to_frame().round(2))
        print("Map modules")
        print(map_data.modules().round(3))
        print("Network edges")
        print(network_data.edges())
        print("Written files")
        for group_name, paths in written.items():
            print(group_name)
            print("\n".join(str(path) for path in paths.values()))


if __name__ == "__main__":
    main()
