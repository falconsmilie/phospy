#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy import PredMatWorkflow
from phospy.io import load_pred_mat
from phospy.workflow import PredMatWorkflowResult


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


def run_demo(
    outdir: Path, *, svm_mode: str = "default"
) -> tuple[PredMatWorkflowResult, Path]:
    repo_root = Path(__file__).resolve().parents[1]
    data_dir = repo_root / "examples" / "data"
    phospho_matrix, substrate_map, site_sequences, motif_sequences = load_demo_inputs(
        data_dir
    )
    outdir.mkdir(parents=True, exist_ok=True)

    workflow = PredMatWorkflow(flank_size=2, svm_mode=svm_mode)
    result = workflow.run(
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

    export_path = result.pred_mat_result.to_csv(outdir / "predMat.csv")
    reloaded_pred_mat = load_pred_mat(export_path)
    pd.testing.assert_frame_equal(
        reloaded_pred_mat,
        result.pred_mat_result.data_frame,
    )
    return result, export_path


def main() -> None:
    with TemporaryDirectory(prefix="phospy-predmat-") as tmp_dir:
        result, export_path = run_demo(Path(tmp_dir), svm_mode="default")
        print(f"Wrote predMat to {export_path} using svm_mode=default")
        print("predMat")
        print(result.pred_mat_result.data_frame.round(4))


if __name__ == "__main__":
    main()
