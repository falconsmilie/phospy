#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from phospy.api import SignalomeRunConfig, SignalomeWorkflow
from phospy.prediction import KinaseScoringResult, PredMatResult
from phospy.signalomes import SignalomeNetworkData, SignalomeResult


def build_demo_inputs() -> tuple[KinaseScoringResult, PredMatResult, pd.DataFrame]:
    site_ids = [
        "PROTEIN_1;S1;",
        "PROTEIN_1;S2;",
        "PROTEIN_2;S3;",
        "PROTEIN_2;S4;",
    ]
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0, 4.0, 4.0],
            "KINASE_B": [1.1, 1.1, 4.1, 4.1],
            "KINASE_C": [4.0, 4.0, 1.0, 1.0],
        },
        index=site_ids,
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.93, 0.20, 0.25],
            "KINASE_B": [0.90, 0.89, 0.94, 0.92],
            "KINASE_C": [0.10, 0.15, 0.91, 0.90],
        },
        index=site_ids,
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 1.1, 3.0, 3.1],
            "sample_2": [1.2, 1.0, 2.9, 3.0],
            "sample_3": [0.9, 1.2, 3.2, 2.8],
        },
        index=site_ids,
    )
    return (
        KinaseScoringResult(
            profile_scores=scoring_matrix,
            combined_scores=scoring_matrix.copy(deep=True),
        ),
        PredMatResult(pred_mat),
        expression_matrix,
    )


def run_demo(
    outdir: Path,
) -> tuple[SignalomeResult, SignalomeNetworkData, dict[str, Path]]:
    scoring_result, prediction_result, expression_matrix = build_demo_inputs()
    signalome_result = SignalomeWorkflow().run(
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A", "KINASE_B"],
        config=SignalomeRunConfig(
            kinase_network_threshold=0.9,
            signalome_cutoff=0.75,
            module_count=2,
        ),
    )
    network_data = signalome_result.to_network_data()
    written = network_data.to_csv(outdir / "signalome_network")
    return signalome_result, network_data, written


def main() -> None:
    with TemporaryDirectory(prefix="phospy-kinase-network-") as tmp_dir:
        signalome_result, network_data, written = run_demo(Path(tmp_dir))
        print(
            f"Wrote signalome network tables to {Path(tmp_dir) / 'signalome_network'}"
        )
        print("Signalome modules")
        print(signalome_result.signalome_modules.round(2))
        print("Network nodes")
        print(network_data.nodes())
        print("Network edges")
        print(network_data.edges())
        print("Written files")
        print("\n".join(str(path) for path in written.values()))


if __name__ == "__main__":
    main()
