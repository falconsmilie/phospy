from __future__ import annotations

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd

from phospy import PredMatWorkflow, SignalomeMapData
from phospy.api import PredictionRunConfig, SignalomeWorkflow
from phospy.signalomes.maps import build_signalome_map_data


def make_workflow_inputs() -> tuple[
    pd.DataFrame,
    dict[str, list[str]],
    dict[str, str],
    dict[str, list[str]],
]:
    phospho_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 1.1, 0.9, 1.2, 3.0, 2.9, 3.1, 2.8],
            "sample_2": [2.0, 2.1, 1.9, 2.2, 2.0, 2.1, 1.9, 2.2],
            "sample_3": [3.0, 3.1, 2.9, 3.2, 1.0, 1.1, 0.9, 1.2],
        },
        index=[f"PROTEIN_{i};S{i};" for i in range(1, 9)],
    )
    substrate_map = {
        "KINASE_A": [
            "PROTEIN_1;S1;",
            "PROTEIN_2;S2;",
            "PROTEIN_3;S3;",
            "PROTEIN_4;S4;",
        ],
        "KINASE_B": [
            "PROTEIN_5;S5;",
            "PROTEIN_6;S6;",
            "PROTEIN_7;S7;",
            "PROTEIN_8;S8;",
        ],
    }
    site_sequences = {
        "PROTEIN_1;S1;": "QQAAAAAYY",
        "PROTEIN_2;S2;": "QQAAAAAYY",
        "PROTEIN_3;S3;": "QQAAAAAYY",
        "PROTEIN_4;S4;": "QQAAAAAYY",
        "PROTEIN_5;S5;": "QQTTTTTYY",
        "PROTEIN_6;S6;": "QQTTTTTYY",
        "PROTEIN_7;S7;": "QQTTTTTYY",
        "PROTEIN_8;S8;": "QQTTTTTYY",
    }
    motif_sequences = {
        "KINASE_A": ["QQAAAAAYY", "QQAAAAAYY", "QQAAAAAYY"],
        "KINASE_B": ["QQTTTTTYY", "QQTTTTTYY", "QQTTTTTYY"],
    }
    return phospho_matrix, substrate_map, site_sequences, motif_sequences


def _build_signalome_result():
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    pred_mat_result = PredMatWorkflow(flank_size=2).run(
        phospho_matrix=phospho_matrix,
        substrate_map=substrate_map,
        site_sequences=site_sequences,
        motif_sequences=motif_sequences,
        prediction_config=PredictionRunConfig(
            min_substrates=2,
            min_motif_size=2,
            ensemble_size=3,
            top=4,
            score_threshold=0.75,
            inclusion=3,
            n_iterations=2,
            random_state=17,
        ),
    )
    return SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A", "KINASE_B"],
    )


def _load_example_module(path: Path):
    spec = spec_from_file_location(path.stem, path)
    assert spec is not None
    assert spec.loader is not None
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_signalome_result_to_map_data_returns_canonical_plotting_tables() -> None:
    result = _build_signalome_result()

    map_data = result.to_map_data()

    assert isinstance(map_data, SignalomeMapData)
    assert list(map_data.modules().columns) == [
        "x",
        "y",
        "n_sites",
        "n_proteins",
        "dominant_kinase",
        "dominant_share_percent",
    ]
    assert map_data.modules().loc[1, "x"] == 0.0
    assert map_data.modules().loc[2, "x"] == 1.0
    assert map_data.modules().loc[1, "n_sites"] == 4
    assert map_data.modules().loc[2, "dominant_kinase"] == "KINASE_B"

    assert list(map_data.sites().columns) == [
        "protein_id",
        "module_id",
        "top_kinase",
        "top_kinase_candidates",
        "top_kinase_tie_count",
        "top_kinase_is_ambiguous",
        "top_score",
        "x",
        "y",
        "module_x",
        "module_y",
        "position_in_module",
        "expression_mean",
        "expression_std",
    ]
    assert map_data.sites().loc["PROTEIN_1;S1;", "module_id"] == 1
    assert map_data.sites().loc["PROTEIN_1;S1;", "module_x"] == 0.0
    assert map_data.sites().loc["PROTEIN_5;S5;", "module_x"] == 1.0

    assert list(map_data.kinases().columns) == [
        "x",
        "y",
        "base_x",
        "module_count",
        "total_share_percent",
        "degree",
        "n_substrates",
        "is_kinase_of_interest",
    ]
    assert map_data.kinases().loc["KINASE_A", "base_x"] == 0.0
    assert map_data.kinases().loc["KINASE_B", "base_x"] == 1.0
    assert bool(map_data.kinases().loc["KINASE_A", "is_kinase_of_interest"])

    assert list(map_data.links().columns) == [
        "kinase",
        "module_id",
        "share_percent",
        "kinase_x",
        "kinase_y",
        "module_x",
        "module_y",
        "is_kinase_of_interest",
    ]
    assert map_data.links().to_dict("records") == [
        {
            "kinase": "KINASE_A",
            "module_id": 1,
            "share_percent": 100.0,
            "kinase_x": 0.0,
            "kinase_y": 1.0,
            "module_x": 0.0,
            "module_y": 0.0,
            "is_kinase_of_interest": True,
        },
        {
            "kinase": "KINASE_B",
            "module_id": 2,
            "share_percent": 100.0,
            "kinase_x": 1.0,
            "kinase_y": 1.0,
            "module_x": 1.0,
            "module_y": 0.0,
            "is_kinase_of_interest": True,
        },
    ]


def test_build_signalome_map_data_matches_result_method() -> None:
    result = _build_signalome_result()

    via_method = result.to_map_data()
    via_function = build_signalome_map_data(result)

    pd.testing.assert_frame_equal(via_method.modules(), via_function.modules())
    pd.testing.assert_frame_equal(via_method.sites(), via_function.sites())
    pd.testing.assert_frame_equal(via_method.kinases(), via_function.kinases())
    pd.testing.assert_frame_equal(via_method.links(), via_function.links())


def test_signalome_map_data_to_frames_and_csv_exports_stable_tables(
    tmp_path: Path,
) -> None:
    result = _build_signalome_result()

    map_data = result.to_map_data()
    frames = map_data.to_frames()

    assert list(frames) == [
        "signalome_map_modules",
        "signalome_map_sites",
        "signalome_map_kinases",
        "signalome_map_links",
    ]

    written = map_data.to_csv(tmp_path)
    assert sorted(written) == [
        "signalome_map_kinases",
        "signalome_map_links",
        "signalome_map_modules",
        "signalome_map_sites",
    ]

    reloaded_modules = pd.read_csv(
        written["signalome_map_modules"],
        index_col=0,
    )
    reloaded_sites = pd.read_csv(
        written["signalome_map_sites"],
        index_col=0,
    )
    reloaded_kinases = pd.read_csv(
        written["signalome_map_kinases"],
        index_col=0,
    )
    reloaded_links = pd.read_csv(written["signalome_map_links"])

    reloaded_modules.index.name = map_data.modules().index.name
    reloaded_sites.index.name = map_data.sites().index.name
    reloaded_kinases.index.name = map_data.kinases().index.name

    reloaded_modules = reloaded_modules.astype(map_data.modules().dtypes.to_dict())
    reloaded_sites = reloaded_sites.astype(map_data.sites().dtypes.to_dict())
    reloaded_kinases = reloaded_kinases.astype(map_data.kinases().dtypes.to_dict())
    reloaded_links = reloaded_links.astype(map_data.links().dtypes.to_dict())

    pd.testing.assert_frame_equal(reloaded_modules, map_data.modules())
    pd.testing.assert_frame_equal(reloaded_sites, map_data.sites())
    pd.testing.assert_frame_equal(reloaded_kinases, map_data.kinases())
    pd.testing.assert_frame_equal(reloaded_links, map_data.links())


def test_signalome_map_demo_runs_end_to_end(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    module = _load_example_module(repo_root / "examples" / "signalome_map_demo.py")

    signalome_result, map_data, written = module.run_demo(tmp_path)

    assert set(written) == {
        "signalome_map_modules",
        "signalome_map_sites",
        "signalome_map_kinases",
        "signalome_map_links",
    }
    assert signalome_result.signalome_modules.shape == (2, 2)
    assert map_data.modules().shape[0] == 2
    assert written["signalome_map_links"].exists()

    reloaded_links = pd.read_csv(written["signalome_map_links"])
    reloaded_links = reloaded_links.astype(map_data.links().dtypes.to_dict())
    pd.testing.assert_frame_equal(reloaded_links, map_data.links())
