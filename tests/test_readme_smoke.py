from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd

from phospy import KinaseActivityAnalyzer, PhosphoDataset, PhosRPipeline
from phospy.constants import (
    CORE_OUTPUT_ARTIFACT_BASENAMES,
    KINASE_OUTPUT_FILENAMES,
)
from phospy.io import load_pred_mat

EXAMPLE_OUTPUT_FILES = {
    *(f"{basename}.csv" for basename in CORE_OUTPUT_ARTIFACT_BASENAMES),
    *KINASE_OUTPUT_FILENAMES,
}


def test_readme_example_analyzer_runs_end_to_end(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parents[1]

    dataset = PhosphoDataset.from_files(
        repo_root / "examples" / "data" / "total.tsv",
        repo_root / "examples" / "data" / "phospho.tsv",
        phospho_encoding="utf-16le",
    )
    core = dataset.preprocessing.run(max_unmatched_fraction=0.1)

    analyzer = KinaseActivityAnalyzer()
    result = analyzer.run(
        pred_mat=analyzer.load_pred_mat(
            repo_root / "examples" / "data" / "predMat.csv"
        ),
        phospho_matrix=core.site_matrix.matrix,
        threshold=0.6,
        min_substrates=1,
        top_n_substrates=1,
    )
    analyzer.write_outputs(result, outdir=tmp_path)

    assert result.target_counts.to_dict() == {"PRKACA": 3, "BTK": 2}
    assert EXAMPLE_OUTPUT_FILES - {
        "df_phospho_corrected.csv",
        "df_phospho_filtered.csv",
        "df_total_filtered.csv",
        "df_total_unique.csv",
        "mat_phospho_corrected.csv",
        "phosr_input.csv",
        "site_sequences.csv",
    } <= {path.name for path in tmp_path.iterdir()}


def test_readme_example_pipeline_runs_end_to_end(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    outdir = tmp_path / "output"

    pipeline = PhosRPipeline.from_files(
        total_path=repo_root / "examples" / "data" / "total.tsv",
        phospho_path=repo_root / "examples" / "data" / "phospho.tsv",
        pred_mat_path=repo_root / "examples" / "data" / "predMat.csv",
        phospho_encoding="utf-16le",
        max_unmatched_fraction=0.1,
        kinase_activity_threshold=0.6,
        kinase_activity_min_substrates=1,
        kinase_activity_top_n_substrates=1,
    )
    outputs = pipeline.run(outdir=outdir)

    assert outputs.core.site_matrix.matrix.index.tolist() == ["BTK;Y551;"]
    assert outputs.kinase_activity is not None
    assert EXAMPLE_OUTPUT_FILES.issubset({path.name for path in outdir.iterdir()})

    expected_matrix = pd.read_csv(
        repo_root / "examples" / "output" / "mat_phospho_corrected.csv",
        index_col=0,
    )
    pd.testing.assert_frame_equal(
        outputs.core.site_matrix.matrix,
        expected_matrix,
        check_index_type=False,
        check_column_type=False,
    )

    expected_target_counts = pd.read_csv(
        repo_root / "examples" / "output" / "kinase_target_counts.csv"
    )
    actual_target_counts = outputs.kinase_activity.target_counts.rename_axis(
        "kinase"
    ).reset_index(name="n_targets")
    pd.testing.assert_frame_equal(actual_target_counts, expected_target_counts)


def _load_example_module(path: Path):
    spec = spec_from_file_location(path.stem, path)
    assert spec is not None
    assert spec.loader is not None
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pred_mat_workflow_demo_runs_end_to_end(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    module = _load_example_module(repo_root / "examples" / "predmat_workflow_demo.py")

    result, export_path = module.run_demo(tmp_path)

    assert export_path.name == "predMat.csv"
    assert export_path.exists()

    reloaded = load_pred_mat(export_path)
    pd.testing.assert_frame_equal(reloaded, result.pred_mat_result.data_frame)
    assert list(result.pred_mat_result.data_frame.columns) == [
        "KINASE_A",
        "KINASE_B",
    ]


def test_signalome_workflow_demo_runs_end_to_end(tmp_path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    module = _load_example_module(repo_root / "examples" / "signalome_workflow_demo.py")

    signalome_result, map_data, network_data, written = module.run_demo(tmp_path)

    assert sorted(written) == ["map", "network", "signalome"]
    assert signalome_result.modules.to_frame().shape == (2, 2)
    assert map_data.modules().shape[0] == 2
    assert list(network_data.nodes().index) == ["KINASE_A", "KINASE_B"]

    signalome_modules_path = written["signalome"]["signalome_modules"]
    signalome_map_modules_path = written["map"]["signalome_map_modules"]
    signalome_network_nodes_path = written["network"]["signalome_network_nodes"]

    assert signalome_modules_path.exists()
    assert signalome_map_modules_path.exists()
    assert signalome_network_nodes_path.exists()

    reloaded_signalome_modules = pd.read_csv(signalome_modules_path, index_col=0)
    reloaded_map_modules = pd.read_csv(signalome_map_modules_path, index_col=0)
    reloaded_network_nodes = pd.read_csv(signalome_network_nodes_path, index_col=0)

    reloaded_signalome_modules.index.name = (
        signalome_result.modules.to_frame().index.name
    )
    reloaded_signalome_modules.columns.name = (
        signalome_result.modules.to_frame().columns.name
    )
    reloaded_map_modules.index.name = map_data.modules().index.name
    reloaded_network_nodes.index.name = network_data.nodes().index.name

    pd.testing.assert_frame_equal(
        reloaded_signalome_modules.astype(
            signalome_result.modules.to_frame().dtypes.to_dict()
        ),
        signalome_result.modules.to_frame(),
    )
    pd.testing.assert_frame_equal(
        reloaded_map_modules.astype(map_data.modules().dtypes.to_dict()),
        map_data.modules(),
    )
    pd.testing.assert_frame_equal(
        reloaded_network_nodes.astype(network_data.nodes().dtypes.to_dict()),
        network_data.nodes(),
    )
