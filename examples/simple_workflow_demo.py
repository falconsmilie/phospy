#!/usr/bin/env python3
"""Simple common-path workflow demo.

Use this lane when you have biologically shaped phospho input and want PhosPy to
handle preprocessing, analysis-ready dataset construction, bundled reference
resolution, predMat generation, and kinase activity analysis.
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy.api import (
    KinaseActivityConfig,
    PredictionRunConfig,
    SimpleKinaseWorkflow,
)
from phospy.api.workflow_results import SimpleKinaseWorkflowResult


def build_demo_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    phospho_df = pd.DataFrame(
        {
            "uid": ["u1", "u2", "u3", "u4", "u5"],
            "gene_names": ["TSC2", "GSK3B", "TBC1D1", "TBC1D1", "EIF4B"],
            "gene_p_site": [
                "TSC2_S939",
                "GSK3B_S9",
                "TBC1D1_T590",
                "TBC1D1_S231",
                "EIF4B_S422",
            ],
            "localization_prob": [0.99, 0.99, 0.99, 0.99, 0.99],
            "centralized_sequence": [
                "HRPVSVHSGSTTLQQFSQPCSRQVTPTPNSP",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                "IITRIQQHLSQQRARSSTPCQGSPEYFTRHVL",
                "ALFVRNQNVQQLHHSSTLPRSLSPPSSQSKGY",
                "SPRRTSRESQVVSQTPRRESEKQESRRRSRSL",
            ],
            "p_group1": [8.0, 7.0, 6.0, 5.0, 4.0],
            "p_group2": [8.5, 7.5, 6.5, 5.5, 4.5],
            "p_group3": [9.0, 8.0, 7.0, 6.0, 5.0],
            "p_group4": [7.5, 6.5, 5.5, 4.5, 3.5],
            "p_group5": [7.0, 6.0, 5.0, 4.0, 3.0],
            "p_group6": [6.5, 5.5, 4.5, 3.5, 2.5],
        }
    )
    total_df = pd.DataFrame(
        {
            "genes": ["TSC2", "GSK3B", "TBC1D1", "EIF4B"],
            "group1": [2.0, 2.1, 2.2, 2.3],
            "group2": [2.1, 2.2, 2.3, 2.4],
            "group3": [2.2, 2.3, 2.4, 2.5],
            "group4": [1.9, 2.0, 2.1, 2.2],
            "group5": [1.8, 1.9, 2.0, 2.1],
            "group6": [1.7, 1.8, 1.9, 2.0],
        }
    )
    return total_df, phospho_df


def run_demo(
    outdir: Path,
    *,
    use_files: bool = True,
) -> tuple[SimpleKinaseWorkflowResult, dict[str, Path]]:
    total_df, phospho_df = build_demo_inputs()
    outdir.mkdir(parents=True, exist_ok=True)

    workflow = SimpleKinaseWorkflow(flank_size=7)
    shared_run_kwargs = {
        "species": "rat",
        "prediction_config": PredictionRunConfig(
            min_substrates=1,
            min_motif_size=1,
            ensemble_size=2,
            top=3,
            inclusion=2,
            n_iterations=2,
            random_state=7,
        ),
        "activity_config": KinaseActivityConfig(
            threshold=0.1,
            min_substrates=1,
            top_n_substrates=3,
        ),
    }

    if use_files:
        total_path = outdir / "total.tsv"
        phospho_path = outdir / "phospho.tsv"
        total_df.to_csv(total_path, sep="\t", index=False)
        phospho_df.to_csv(phospho_path, sep="\t", index=False)
        result = workflow.run(
            total=total_path,
            phospho=phospho_path,
            **shared_run_kwargs,
        )
    else:
        result = workflow.run(
            total=total_df,
            phospho=phospho_df,
            **shared_run_kwargs,
        )

    written = {
        "pred_mat": result.pred_mat_result.to_csv(outdir / "predMat.csv"),
        "weighted_activity": outdir / "weighted_activity.csv",
        "ksea_scores": outdir / "ksea_scores.csv",
    }
    result.kinase_activity_result.weighted_activity.to_csv(written["weighted_activity"])
    result.kinase_activity_result.ksea_scores.to_csv(written["ksea_scores"])
    return result, written


def main() -> None:
    with TemporaryDirectory(prefix="phospy-simple-workflow-") as tmp_dir:
        result, written = run_demo(Path(tmp_dir), use_files=True)
        print("Simple workflow demo")
        print("Reference lane")
        print(
            {
                "species": result.reference_bundle.species,
                "reference": result.reference_bundle.source_metadata.reference,
            }
        )
        print()
        print("predMat")
        print(result.pred_mat_result.to_frame(copy=False).round(4))
        print()
        print("Weighted activity")
        print(result.kinase_activity_result.weighted_activity.round(4))
        print()
        print("Written files")
        print("\n".join(str(path) for path in written.values()))


if __name__ == "__main__":
    main()
