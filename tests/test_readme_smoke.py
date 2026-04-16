from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd
import pandas.testing as pdt

from phospy.activities import KinaseActivityAnalyzer
from phospy.api import PredictionRunConfig, SimpleKinaseWorkflow
from phospy.datasets import PhosphoDataset
from phospy.internal.constants import KINASE_OUTPUT_FILENAMES
from phospy.preprocessing import CorePreprocessingConfig

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_WORKFLOW_REFERENCE_DIR = (
    ROOT / "tests" / "fixtures" / "public_workflow_reference"
)


def _read_indexed_fixture(name: str) -> pd.DataFrame:
    return pd.read_csv(PUBLIC_WORKFLOW_REFERENCE_DIR / name, index_col=0)


def test_readme_example_analyzer_runs_end_to_end(tmp_path) -> None:
    dataset = PhosphoDataset.from_files(
        ROOT / "examples" / "data" / "total.tsv",
        ROOT / "examples" / "data" / "phospho.tsv",
        phospho_encoding="utf-16le",
    )
    core = dataset.preprocessing.run(
        config=CorePreprocessingConfig(max_unmatched_fraction=0.1)
    )

    analyzer = KinaseActivityAnalyzer()
    result = analyzer.run(
        pred_mat=analyzer.load_pred_mat(ROOT / "examples" / "data" / "predMat.csv"),
        phospho_matrix=core.site_matrix.matrix,
        threshold=0.6,
        min_substrates=1,
        top_n_substrates=1,
    )
    analyzer.write_outputs(result, outdir=tmp_path)

    assert result.target_counts.to_dict() == {"PRKACA": 3, "BTK": 2}
    assert set(KINASE_OUTPUT_FILENAMES) <= {path.name for path in tmp_path.iterdir()}


def test_readme_example_simple_workflow_runs_end_to_end() -> None:
    with SimpleKinaseWorkflow(flank_size=7).run(
        total=ROOT / "examples" / "data" / "simple_workflow" / "total.tsv",
        phospho=ROOT / "examples" / "data" / "simple_workflow" / "phospho.tsv",
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
    ) as result:
        assert result.reference_bundle.source_metadata.reference == "l6_native"
        assert not result.pred_mat_result.to_frame(copy=False).empty
        assert isinstance(result.prediction_result.substrate_list, dict)
        assert not result.scoring_result.profile_scores.empty
        assert result.scoring_result.combined_scores is not None
        assert not result.kinase_activity_result.weighted_activity.empty
        assert not result.kinase_activity_result.ksea_scores.empty
        assert result.pred_mat_result is result.prediction_result.pred_mat_result
        assert result.profile_scores is result.scoring_result.profile_scores
        assert result.combined_scores is result.scoring_result.combined_scores
        assert result.weights is result.scoring_result.weights
        assert result.substrate_list is result.prediction_result.substrate_list
        assert not hasattr(result, "pred_mat")

        actual_pred_mat = result.pred_mat_result.to_frame(copy=False).copy(deep=True)
        actual_weighted_activity = result.kinase_activity_result.weighted_activity.copy(
            deep=True
        )
        actual_ksea_scores = result.kinase_activity_result.ksea_scores.copy(deep=True)
        actual_ksea_counts = result.kinase_activity_result.ksea_counts.to_frame(
            name="n_substrates"
        )
        actual_target_counts = result.kinase_activity_result.target_counts.to_frame(
            name="n_targets"
        )

    expected_pred_mat = _read_indexed_fixture("simple_workflow_predmat.csv")
    expected_weighted_activity = _read_indexed_fixture(
        "simple_workflow_weighted_activity.csv"
    )
    expected_ksea_scores = _read_indexed_fixture("simple_workflow_ksea_scores.csv")
    expected_ksea_counts = _read_indexed_fixture("simple_workflow_ksea_counts.csv")
    expected_target_counts = _read_indexed_fixture("simple_workflow_target_counts.csv")
    expected_pred_mat.index = expected_pred_mat.index.astype(
        actual_pred_mat.index.dtype
    )

    pdt.assert_frame_equal(actual_pred_mat, expected_pred_mat)
    pdt.assert_frame_equal(actual_weighted_activity, expected_weighted_activity)
    pdt.assert_frame_equal(actual_ksea_scores, expected_ksea_scores)
    pdt.assert_frame_equal(actual_ksea_counts, expected_ksea_counts)
    pdt.assert_frame_equal(actual_target_counts, expected_target_counts)


def _load_example_module(path: Path):
    spec = spec_from_file_location(path.stem, path)
    assert spec is not None
    assert spec.loader is not None
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_simple_workflow_demo_runs_end_to_end(tmp_path) -> None:
    module = _load_example_module(ROOT / "examples" / "simple_workflow_demo.py")

    result, written = module.run_demo(tmp_path / "simple", use_files=True)

    assert result.reference_bundle.species == "rat"
    assert result.reference_bundle.source_metadata.reference == "l6_native"
    assert not result.pred_mat_result.to_frame(copy=False).empty
    assert result.prediction_result.pred_mat_result is result.pred_mat_result
    assert isinstance(result.prediction_result.substrate_list, dict)
    assert not result.scoring_result.profile_scores.empty
    assert result.scoring_result.combined_scores is not None
    assert not result.kinase_activity_result.weighted_activity.empty
    assert not result.kinase_activity_result.ksea_scores.empty
    assert set(written) == {"pred_mat", "weighted_activity", "ksea_scores"}
    assert all(path.exists() for path in written.values())


def test_signalome_workflow_demo_runs_end_to_end(tmp_path) -> None:
    module = _load_example_module(ROOT / "examples" / "signalome_workflow_demo.py")

    signalome_result, map_data, network_data, written = module.run_demo(tmp_path)

    assert sorted(written) == ["map", "network", "signalome"]
    assert signalome_result.modules.to_frame().shape[0] >= 1
    assert map_data.modules().shape[0] >= 1
    assert network_data.nodes().shape[0] >= 1

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
