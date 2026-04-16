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
from phospy.signalomes import SignalomeMapData, SignalomeResult


def run_demo(
    outdir: Path,
) -> tuple[SignalomeResult, SignalomeMapData, dict[str, Path]]:
    repo_root = Path(__file__).resolve().parents[1]

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
    kinases_of_interest = list(
        simple_result.prediction_result.pred_mat_result.kinase_names[:2]
    )
    signalome_result = SignalomeWorkflow().run_from_analysis_ready(
        dataset=simple_result.analysis_ready_dataset,
        scoring_result=simple_result.scoring_result,
        prediction_result=simple_result.prediction_result,
        kinases_of_interest=kinases_of_interest,
        metadata_fallback_policy="metadata",
        allow_gene_symbol_fallback=True,
        config=SignalomeRunConfig(signalome_cutoff=0.5),
    )
    map_data = signalome_result.to_map_data()
    written = map_data.to_csv(outdir / "signalome_map")
    simple_result.close()
    return signalome_result, map_data, written


def main() -> None:
    with TemporaryDirectory(prefix="phospy-signalome-map-") as tmp_dir:
        signalome_result, map_data, written = run_demo(Path(tmp_dir))
        print(f"Wrote signalome map tables to {Path(tmp_dir) / 'signalome_map'}")
        print("Signalome modules")
        print(signalome_result.modules.to_frame().round(2))
        print("Map modules")
        print(map_data.modules().round(3))
        print("Written files")
        print("\n".join(str(path) for path in written.values()))


if __name__ == "__main__":
    main()
