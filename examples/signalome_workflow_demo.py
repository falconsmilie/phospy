#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from phospy.api import (
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
    SimpleKinaseWorkflow,
)
from phospy.signalomes import (
    SignalomeMapData,
    SignalomeNetworkData,
    SignalomeResult,
)


def run_demo(
    outdir: Path,
) -> tuple[
    SignalomeResult,
    SignalomeMapData,
    SignalomeNetworkData,
    dict[str, dict[str, Path]],
]:
    repo_root = Path(__file__).resolve().parents[1]
    outdir.mkdir(parents=True, exist_ok=True)

    simple_result = SimpleKinaseWorkflow(flank_size=7).run(
        total=repo_root / "examples" / "data" / "simple_workflow" / "total.tsv",
        phospho=repo_root / "examples" / "data" / "simple_workflow" / "phospho.tsv",
        species="rat",
        prediction_config=PredictionRunConfig(
            min_substrates=1,
            min_motif_size=1,
            ensemble_size=2,
            top=3,
            inclusion=2,
            n_iterations=2,
            random_state=7,
        ),
    )
    pred_mat_result = simple_result.pred_mat_result
    prediction_result = simple_result.prediction_result
    scoring_result = simple_result.scoring_result

    kinases_of_interest = list(pred_mat_result.kinase_names[:2])
    signalome_result = SignalomeWorkflow().run_from_analysis_ready(
        dataset=simple_result.analysis_ready_dataset,
        scoring_result=scoring_result,
        prediction_result=prediction_result,
        kinases_of_interest=kinases_of_interest,
        metadata_fallback_policy="metadata",
        allow_gene_symbol_fallback=True,
        config=SignalomeRunConfig(signalome_cutoff=0.5),
    )
    map_data = signalome_result.to_map_data()
    network_data = signalome_result.to_network_data()

    written = {
        "signalome": signalome_result.to_csv(outdir / "signalome"),
        "map": map_data.to_csv(outdir / "signalome_map"),
        "network": network_data.to_csv(outdir / "signalome_network"),
    }

    simple_result.close()
    return signalome_result, map_data, network_data, written


def main() -> None:
    with TemporaryDirectory(prefix="phospy-signalome-workflow-") as tmp_dir:
        signalome_result, map_data, network_data, written = run_demo(Path(tmp_dir))
        print(f"Wrote signalome workflow tables to {tmp_dir}")
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
