from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import phospy.signalomes.clustering as signalome_clustering
from phospy import PredMatResult, SignalomeResult
from phospy.api import (
    PredictionRunConfig,
    PredMatWorkflow,
    SignalomeRunConfig,
    SignalomeWorkflow,
)
from phospy.errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    RequestValidationError,
    TableSchemaError,
)
from phospy.signalomes import (
    SignalomeModuleSelectionPolicy,
    build_signalome_result,
)
from phospy.signalomes.analysis import build_kinase_network_view
from phospy.signalomes.assignments import build_expanded_signalomes


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


def _build_pred_mat_workflow_result():
    phospho_matrix, substrate_map, site_sequences, motif_sequences = (
        make_workflow_inputs()
    )
    result = PredMatWorkflow(flank_size=2).run(
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
    return phospho_matrix, result


def test_signalome_workflow_constructs_signalomes_from_scoring_and_prediction_results() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    assert isinstance(result, SignalomeResult)
    assert list(result.pred_mat.columns) == ["KINASE_A", "KINASE_B"]
    assert list(result.signalome_modules.columns) == ["KINASE_A", "KINASE_B"]
    assert set(result.kinase_substrates) == {"KINASE_A", "KINASE_B"}
    assert set(result.site_assignments.columns) == {
        "protein_id",
        "module_id",
        "top_kinase",
        "top_kinase_candidates",
        "top_kinase_tie_count",
        "top_kinase_is_ambiguous",
        "top_score",
    }
    assert set(result.site_assignments.loc[:, "module_id"]) == {1, 2}
    assert result.signalome_modules.loc[1, "KINASE_A"] == 100.0
    assert result.signalome_modules.loc[2, "KINASE_B"] == 100.0

    expanded = result.expanded_signalomes["KINASE_A"]
    assert expanded.kinase == "KINASE_A"
    assert expanded.linked_kinases[0] == "KINASE_A"
    assert expanded.regulated_module_ids == (1,)
    assert set(expanded.expression_matrix.index) == {
        "PROTEIN_1;S1;",
        "PROTEIN_2;S2;",
        "PROTEIN_3;S3;",
        "PROTEIN_4;S4;",
    }


def test_build_kinase_network_view_builds_expected_neighbor_map_and_edges() -> None:
    kinase_network = pd.DataFrame(
        {
            "KINASE_A": [0.0, 0.95, 0.0],
            "KINASE_B": [0.95, 0.0, 0.91],
            "KINASE_C": [0.0, 0.91, 0.0],
        },
        index=["KINASE_A", "KINASE_B", "KINASE_C"],
    )
    kinase_correlation_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 0.95, 0.45],
            "KINASE_B": [0.95, 1.0, 0.91],
            "KINASE_C": [0.45, 0.91, 1.0],
        },
        index=kinase_network.index.copy(),
    )

    network = build_kinase_network_view(
        kinase_network=kinase_network,
        kinase_correlation_matrix=kinase_correlation_matrix,
        kinase_substrates={
            "KINASE_A": ("SITE_1", "SITE_2"),
            "KINASE_B": ("SITE_3",),
            "KINASE_C": (),
        },
    )

    assert network.neighbor_map == {
        "KINASE_A": ("KINASE_B",),
        "KINASE_B": ("KINASE_A", "KINASE_C"),
        "KINASE_C": ("KINASE_B",),
    }
    pd.testing.assert_frame_equal(
        network.node_table,
        pd.DataFrame(
            {
                "degree": [1, 2, 1],
                "n_substrates": [2, 1, 0],
            },
            index=pd.Index(["KINASE_A", "KINASE_B", "KINASE_C"], name="kinase"),
        ).astype({"degree": int, "n_substrates": int}),
    )
    pd.testing.assert_frame_equal(
        network.edge_table,
        pd.DataFrame(
            {
                "source_kinase": ["KINASE_A", "KINASE_B"],
                "target_kinase": ["KINASE_B", "KINASE_C"],
                "correlation": [0.95, 0.91],
            }
        ).astype(
            {
                "source_kinase": str,
                "target_kinase": str,
                "correlation": float,
            }
        ),
    )


def test_build_expanded_signalomes_uses_neighbor_map_and_preserves_site_order() -> None:
    site_assignments = pd.DataFrame(
        {
            "protein_id": ["PROTEIN_1", "PROTEIN_2", "PROTEIN_3", "PROTEIN_4"],
            "module_id": [1, 2, 1, 3],
            "top_kinase": ["KINASE_A", "KINASE_B", "KINASE_B", "KINASE_C"],
            "top_kinase_candidates": [
                '["KINASE_A"]',
                '["KINASE_B"]',
                '["KINASE_B"]',
                '["KINASE_C"]',
            ],
            "top_kinase_tie_count": [1, 1, 1, 1],
            "top_kinase_is_ambiguous": [False, False, False, False],
            "top_score": [0.95, 0.92, 0.91, 0.90],
        },
        index=pd.Index(["SITE_3", "SITE_1", "SITE_4", "SITE_2"], name="site_id"),
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [3.0, 1.0, 4.0, 2.0],
            "sample_2": [3.1, 1.1, 4.1, 2.1],
        },
        index=site_assignments.index.copy(),
    )
    signalome_modules = pd.DataFrame(
        {
            "KINASE_A": [10.0, 80.0, 0.0],
            "KINASE_B": [90.0, 5.0, 0.0],
            "KINASE_C": [0.0, 0.0, 95.0],
        },
        index=pd.Index([1, 2, 3], name="module_id"),
    )

    expanded = build_expanded_signalomes(
        kinases_of_interest=["KINASE_A"],
        kinase_network={"KINASE_A": ("KINASE_B",)},
        kinase_substrates={
            "KINASE_A": ("SITE_2", "SITE_1"),
            "KINASE_B": ("SITE_3", "SITE_4"),
        },
        signalome_modules=signalome_modules,
        site_assignments=site_assignments,
        expression_matrix=expression_matrix,
        min_kinase_module_share_percent=50.0,
    )

    kinase_a = expanded["KINASE_A"]
    assert kinase_a.linked_kinases == ("KINASE_A", "KINASE_B")
    assert kinase_a.regulated_module_ids == (2,)
    assert kinase_a.site_assignments.index.tolist() == ["SITE_1"]
    pd.testing.assert_frame_equal(
        kinase_a.expression_matrix,
        expression_matrix.loc[["SITE_1"]],
    )


def test_signalome_workflow_accepts_canonical_pred_mat_result_input() -> None:
    phospho_matrix, pred_mat_workflow_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_workflow_result.scoring_result,
        prediction_result=pred_mat_workflow_result.pred_mat_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_B"],
    )

    assert result.kinases_of_interest == ("KINASE_B",)
    assert result.expanded_signalomes["KINASE_B"].regulated_module_ids == (2,)


def test_signalome_workflow_rejects_pred_mat_without_candidate_kinases() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    empty_pred_mat = PredMatResult(
        pred_mat_result.pred_mat_result.to_frame(copy=True).iloc[:, 0:0]
    )

    with pytest.raises(
        NoCandidateKinasesError,
        match=(
            "prediction_result does not contain any kinase columns because no "
            "candidate kinases qualified for prediction"
        ),
    ):
        SignalomeWorkflow().run(
            scoring_result=pred_mat_result.scoring_result,
            prediction_result=empty_pred_mat,
            expression_matrix=phospho_matrix,
            kinases_of_interest=["KINASE_A"],
        )


def test_signalome_workflow_accepts_explicit_site_to_protein_mapping() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    renamed_index = [f"SITE_{i}" for i in range(1, phospho_matrix.shape[0] + 1)]
    renamed_expression_matrix = phospho_matrix.copy()
    renamed_expression_matrix.index = renamed_index

    scoring_result = pred_mat_result.scoring_result
    scoring_result.combined_scores.index = renamed_index
    scoring_result.profile_scores.index = renamed_index

    pred_mat = pred_mat_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = renamed_index
    renamed_pred_mat_result = PredMatResult(pred_mat)

    site_to_protein = {
        "SITE_1": "PROTEIN_A",
        "SITE_2": "PROTEIN_A",
        "SITE_3": "PROTEIN_B",
        "SITE_4": "PROTEIN_B",
        "SITE_5": "PROTEIN_C",
        "SITE_6": "PROTEIN_C",
        "SITE_7": "PROTEIN_D",
        "SITE_8": "PROTEIN_D",
    }

    result = SignalomeWorkflow().run(
        scoring_result=scoring_result,
        prediction_result=renamed_pred_mat_result,
        expression_matrix=renamed_expression_matrix,
        kinases_of_interest=["KINASE_A"],
        site_to_protein=site_to_protein,
    )

    assert sorted(result.protein_assignments.index.tolist()) == [
        "PROTEIN_A",
        "PROTEIN_B",
        "PROTEIN_C",
        "PROTEIN_D",
    ]
    assert result.site_assignments.loc["SITE_1", "protein_id"] == "PROTEIN_A"
    assert result.site_assignments.loc["SITE_2", "protein_id"] == "PROTEIN_A"


def test_signalome_workflow_rejects_unsupported_site_identifier_format_without_mapping() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    renamed_index = [f"SITE_{i}" for i in range(1, phospho_matrix.shape[0] + 1)]
    renamed_expression_matrix = phospho_matrix.copy()
    renamed_expression_matrix.index = renamed_index

    scoring_result = pred_mat_result.scoring_result
    scoring_result.combined_scores.index = renamed_index
    scoring_result.profile_scores.index = renamed_index

    pred_mat = pred_mat_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = renamed_index
    renamed_pred_mat_result = PredMatResult(pred_mat)

    with pytest.raises(
        InputCompatibilityError,
        match=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in the supported 'PROTEIN;SITE;...' "
            "format"
        ),
    ):
        SignalomeWorkflow().run(
            scoring_result=scoring_result,
            prediction_result=renamed_pred_mat_result,
            expression_matrix=renamed_expression_matrix,
            kinases_of_interest=["KINASE_A"],
        )


def test_signalome_workflow_rejects_incomplete_site_to_protein_mapping() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    with pytest.raises(
        InputCompatibilityError,
        match="site_to_protein must define a protein ID for every aligned phosphosite row",
    ):
        SignalomeWorkflow().run(
            scoring_result=pred_mat_result.scoring_result,
            prediction_result=pred_mat_result.prediction_result,
            expression_matrix=phospho_matrix,
            kinases_of_interest=["KINASE_A"],
            site_to_protein={"PROTEIN_1;S1;": "PROTEIN_1"},
        )


def test_signalome_workflow_rejects_malformed_supported_site_identifier_without_mapping() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    malformed_index = [
        "PROTEIN_1;banana;",
        *[f"PROTEIN_{i};S{i};" for i in range(2, phospho_matrix.shape[0] + 1)],
    ]
    malformed_expression_matrix = phospho_matrix.copy()
    malformed_expression_matrix.index = malformed_index

    scoring_result = pred_mat_result.scoring_result
    scoring_result.combined_scores.index = malformed_index
    scoring_result.profile_scores.index = malformed_index

    pred_mat = pred_mat_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = malformed_index
    malformed_pred_mat_result = PredMatResult(pred_mat)

    with pytest.raises(
        InputCompatibilityError,
        match=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in the supported 'PROTEIN;SITE;...' "
            "format"
        ),
    ):
        SignalomeWorkflow().run(
            scoring_result=scoring_result,
            prediction_result=malformed_pred_mat_result,
            expression_matrix=malformed_expression_matrix,
            kinases_of_interest=["KINASE_A"],
        )


def test_build_signalome_result_rejects_malformed_supported_site_identifier_without_mapping() -> (
    None
):
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0],
            "KINASE_B": [1.0, 1.0],
        },
        index=["PROTEIN_1;banana;", "PROTEIN_2;S2;"],
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.9, 0.2],
            "KINASE_B": [0.1, 0.8],
        },
        index=scoring_matrix.index.copy(),
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [1.1, 2.1],
        },
        index=scoring_matrix.index.copy(),
    )

    with pytest.raises(
        InputCompatibilityError,
        match=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in the supported 'PROTEIN;SITE;...' "
            "format"
        ),
    ):
        build_signalome_result(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            kinases_of_interest=["KINASE_A"],
            signalome_cutoff=0.5,
            module_count=1,
        )


def test_build_signalome_result_uses_explicit_site_to_protein_mapping_for_grouping() -> (
    None
):
    site_ids = ["SITE_1", "SITE_2", "SITE_3", "SITE_4"]
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0, 4.0, 4.0],
            "KINASE_B": [1.1, 1.1, 4.1, 4.1],
        },
        index=site_ids,
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.93, 0.20, 0.25],
            "KINASE_B": [0.10, 0.12, 0.91, 0.90],
        },
        index=site_ids,
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 1.1, 3.0, 3.1],
            "sample_2": [1.2, 1.0, 2.9, 3.0],
        },
        index=site_ids,
    )

    result = build_signalome_result(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A"],
        site_to_protein={
            "SITE_1": "PROTEIN_1",
            "SITE_2": "PROTEIN_1",
            "SITE_3": "PROTEIN_2",
            "SITE_4": "PROTEIN_2",
        },
        signalome_cutoff=0.5,
        module_count=2,
    )

    assert sorted(result.protein_assignments.index.tolist()) == [
        "PROTEIN_1",
        "PROTEIN_2",
    ]
    assert result.protein_assignments.loc["PROTEIN_1", "site_count"] == 2
    assert result.site_assignments.loc["SITE_1", "protein_id"] == "PROTEIN_1"


def test_build_site_assignments_tracks_tied_top_kinases_deterministically() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0],
            "KINASE_B": [1.0, 1.0],
        },
        index=["PROTEIN_1;S1;", "PROTEIN_2;S2;"],
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_B": [0.8, 0.2],
            "KINASE_A": [0.8, 0.9],
        },
        index=scoring_matrix.index.copy(),
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [1.1, 2.1],
        },
        index=scoring_matrix.index.copy(),
    )

    result = build_signalome_result(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A"],
        signalome_cutoff=0.5,
        module_count=1,
    )

    tied_row = result.site_assignments.loc["PROTEIN_1;S1;"]
    clear_row = result.site_assignments.loc["PROTEIN_2;S2;"]

    assert tied_row["top_kinase"] == "KINASE_A"
    assert tied_row["top_kinase_candidates"] == '["KINASE_A", "KINASE_B"]'
    assert tied_row["top_kinase_tie_count"] == 2
    assert bool(tied_row["top_kinase_is_ambiguous"])

    assert clear_row["top_kinase"] == "KINASE_A"
    assert clear_row["top_kinase_candidates"] == '["KINASE_A"]'
    assert clear_row["top_kinase_tie_count"] == 1
    assert not bool(clear_row["top_kinase_is_ambiguous"])


def test_select_module_count_builds_one_cluster_tree_for_candidate_scoring(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.1, 2.0, 2.1, 3.0, 3.1],
            "KINASE_B": [1.2, 1.0, 2.2, 2.0, 3.2, 3.0],
            "KINASE_C": [0.9, 1.2, 1.9, 2.2, 2.9, 3.2],
        },
        index=[f"PROTEIN_{i};S{i};" for i in range(1, 7)],
    )

    observed_tree_builds: list[int] = []
    original_build_cluster_tree = signalome_clustering.build_cluster_tree

    def counting_build_cluster_tree(scoring_values: object) -> object:
        observed_tree_builds.append(1)
        return original_build_cluster_tree(scoring_values)

    monkeypatch.setattr(
        signalome_clustering,
        "build_cluster_tree",
        counting_build_cluster_tree,
    )

    selected = signalome_clustering.select_module_count(
        scoring_matrix.to_numpy(dtype=float)
    )

    assert selected >= 1
    assert observed_tree_builds == [1]


def test_build_site_assignments_is_stable_when_pred_mat_columns_are_reordered() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0],
            "KINASE_B": [1.0, 1.0],
        },
        index=["PROTEIN_1;S1;", "PROTEIN_2;S2;"],
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [1.1, 2.1],
        },
        index=scoring_matrix.index.copy(),
    )
    pred_mat_left = pd.DataFrame(
        {
            "KINASE_A": [0.8, 0.9],
            "KINASE_B": [0.8, 0.2],
        },
        index=scoring_matrix.index.copy(),
    )
    pred_mat_right = pred_mat_left.loc[:, ["KINASE_B", "KINASE_A"]]

    left = build_signalome_result(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat_left,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A"],
        signalome_cutoff=0.5,
        module_count=1,
    )
    right = build_signalome_result(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat_right,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A"],
        signalome_cutoff=0.5,
        module_count=1,
    )

    pd.testing.assert_frame_equal(left.site_assignments, right.site_assignments)


def test_signalome_result_exposes_canonical_module_assignment_and_network_views() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    assert result.modules.to_frame(copy=False) is result.signalome_modules
    assert list(result.modules.to_relationship_table().columns) == [
        "module_id",
        "kinase",
        "share_percent",
    ]
    assert result.modules.to_relationship_table().to_dict("records") == [
        {"module_id": 1, "kinase": "KINASE_A", "share_percent": 100.0},
        {"module_id": 2, "kinase": "KINASE_B", "share_percent": 100.0},
    ]

    assert result.assignments.sites(copy=False) is result.site_assignments
    assert list(result.assignments.proteins().columns) == ["module_id", "site_count"]
    assert result.assignments.proteins().loc["PROTEIN_1", "module_id"] == 1
    assert result.assignments.proteins().loc["PROTEIN_1", "site_count"] == 1

    assert list(result.network.nodes().columns) == ["degree", "n_substrates"]
    assert result.network.nodes().loc["KINASE_A", "n_substrates"] == 4
    assert list(result.network.edges().columns) == [
        "source_kinase",
        "target_kinase",
        "correlation",
    ]


def test_signalome_result_to_frames_returns_stable_named_outputs() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    frames = result.to_frames()

    assert list(frames) == [
        "signalome_modules",
        "kinase_module_relationships",
        "site_assignments",
        "protein_assignments",
        "kinase_network_nodes",
        "kinase_network_edges",
        "kinase_correlation_matrix",
    ]
    assert "scoring_matrix" not in frames
    assert frames["signalome_modules"].equals(result.signalome_modules)
    assert frames["protein_assignments"].equals(result.protein_assignments)

    frames_with_inputs = result.to_frames(include_inputs=True)
    assert list(frames_with_inputs)[-3:] == [
        "scoring_matrix",
        "pred_mat",
        "expression_matrix",
    ]


def test_signalome_result_to_csv_exports_canonical_tables(tmp_path: Path) -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    written = result.to_csv(tmp_path)

    assert sorted(written) == [
        "kinase_correlation_matrix",
        "kinase_module_relationships",
        "kinase_network_edges",
        "kinase_network_nodes",
        "protein_assignments",
        "signalome_modules",
        "site_assignments",
    ]

    reloaded_signalome_modules = pd.read_csv(
        written["signalome_modules"],
        index_col=0,
    ).astype(float)
    reloaded_site_assignments = pd.read_csv(
        written["site_assignments"],
        index_col=0,
    )
    reloaded_protein_assignments = pd.read_csv(
        written["protein_assignments"],
        index_col=0,
    )

    reloaded_signalome_modules.index.name = result.signalome_modules.index.name
    reloaded_signalome_modules.columns.name = result.signalome_modules.columns.name
    pd.testing.assert_frame_equal(reloaded_signalome_modules, result.signalome_modules)
    pd.testing.assert_frame_equal(reloaded_site_assignments, result.site_assignments)
    pd.testing.assert_frame_equal(
        reloaded_protein_assignments, result.protein_assignments
    )


def test_signalome_workflow_rejects_empty_kinases_of_interest() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    with pytest.raises(RequestValidationError, match="kinases_of_interest"):
        SignalomeWorkflow().run(
            scoring_result=pred_mat_result.scoring_result,
            prediction_result=pred_mat_result.prediction_result,
            expression_matrix=phospho_matrix,
            kinases_of_interest=[],
        )


def test_signalome_workflow_rejects_unknown_kinases_of_interest() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    with pytest.raises(InputCompatibilityError, match="kinases_of_interest"):
        SignalomeWorkflow().run(
            scoring_result=pred_mat_result.scoring_result,
            prediction_result=pred_mat_result.prediction_result,
            expression_matrix=phospho_matrix,
            kinases_of_interest=["KINASE_X"],
        )


def test_build_signalome_result_rejects_non_finite_expression_values() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0],
            "KINASE_B": [1.0, 1.0],
        },
        index=["PROTEIN_1;S1;", "PROTEIN_2;S2;"],
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.9, 0.2],
            "KINASE_B": [0.1, 0.8],
        },
        index=scoring_matrix.index.copy(),
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, float("nan")],
            "sample_2": [1.1, 2.1],
        },
        index=scoring_matrix.index.copy(),
    )

    with pytest.raises(
        TableSchemaError, match="expression_matrix contains non-finite values"
    ):
        build_signalome_result(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            kinases_of_interest=["KINASE_A"],
            signalome_cutoff=0.5,
            module_count=1,
        )


def test_build_signalome_result_reports_domain_error_for_unaligned_inputs() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0],
            "KINASE_B": [1.0, 1.0],
        },
        index=["PROTEIN_1;S1;", "PROTEIN_2;S2;"],
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.9, 0.2],
            "KINASE_B": [0.1, 0.8],
        },
        index=["OTHER_1;S1;", "OTHER_2;S2;"],
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [1.1, 2.1],
        },
        index=scoring_matrix.index.copy(),
    )

    with pytest.raises(
        InputCompatibilityError,
        match=(
            "scoring_matrix, pred_mat, and expression_matrix must share at least "
            "one phosphosite row"
        ),
    ):
        build_signalome_result(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            kinases_of_interest=["KINASE_A"],
            signalome_cutoff=0.5,
            module_count=1,
        )


def test_build_signalome_result_reports_domain_error_for_unaligned_kinase_columns() -> (
    None
):
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0],
            "KINASE_B": [1.0, 1.0],
        },
        index=["PROTEIN_1;S1;", "PROTEIN_2;S2;"],
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_X": [0.9, 0.2],
            "KINASE_Y": [0.1, 0.8],
        },
        index=scoring_matrix.index.copy(),
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [1.1, 2.1],
        },
        index=scoring_matrix.index.copy(),
    )

    with pytest.raises(
        InputCompatibilityError,
        match="scoring_matrix and pred_mat must share at least one kinase column",
    ):
        build_signalome_result(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            kinases_of_interest=["KINASE_A"],
            signalome_cutoff=0.5,
            module_count=1,
        )


def test_build_signalome_result_rejects_non_finite_pred_mat_values() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0],
            "KINASE_B": [1.0, 1.0],
        },
        index=["PROTEIN_1;S1;", "PROTEIN_2;S2;"],
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.9, float("nan")],
            "KINASE_B": [0.1, 0.8],
        },
        index=scoring_matrix.index.copy(),
    )
    expression_matrix = pd.DataFrame(
        {
            "sample_1": [1.0, 2.0],
            "sample_2": [1.1, 2.1],
        },
        index=scoring_matrix.index.copy(),
    )

    with pytest.raises(
        TableSchemaError,
        match="pred_mat contains non-finite values in numeric columns",
    ):
        build_signalome_result(
            scoring_matrix=scoring_matrix,
            pred_mat=pred_mat,
            expression_matrix=expression_matrix,
            kinases_of_interest=["KINASE_A"],
            signalome_cutoff=0.5,
            module_count=1,
        )


def test_build_signalome_result_exposes_module_selection_diagnostics() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    assert result.module_selection_diagnostics.used_automatic_selection
    assert result.module_selection_diagnostics.selected_module_count >= 1
    assert result.module_selection_diagnostics.strategy == "correlation_thresholds"
    assert result.module_selection_diagnostics.reason


def test_signalome_workflow_accepts_explicit_module_selection_policy() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
        config=SignalomeRunConfig(
            module_selection_policy=SignalomeModuleSelectionPolicy(
                strategy="single_module"
            )
        ),
    )

    assert result.module_selection_diagnostics.strategy == "single_module"
    assert result.module_selection_diagnostics.selected_module_count == 1
    assert set(result.site_assignments.loc[:, "module_id"]) == {1}
