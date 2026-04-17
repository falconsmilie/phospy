from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import phospy.signalomes.analysis as signalome_analysis
import phospy.signalomes.assignments as signalome_assignments
import phospy.signalomes.clustering as signalome_clustering
import pytest
from phospy.internal.kinase_workflows import KinaseWorkflow
from phospy.signalomes.analysis import (
    build_kinase_network,
    build_kinase_network_view,
    build_signalome_support_matrix,
)
from phospy.signalomes.assignments import (
    build_expanded_signalomes,
    build_site_assignments,
)
from phospy.signalomes.serialization import serialize_site_assignments_for_export

from phospy.api import (
    PredictionRunConfig,
    SignalomeRunConfig,
    SignalomeWorkflow,
)
from phospy.datasets import (
    AnalysisReadyPhosphoDataset,
    AnalysisReadyPreprocessingProvenance,
    AnalysisReadyRowCounts,
    AnalysisReadySiteMatrixStats,
    DatasetSchema,
)
from phospy.errors import (
    InputCompatibilityError,
    NoCandidateKinasesError,
    RequestValidationError,
    TableSchemaError,
)
from phospy.prediction import PredMatResult
from phospy.signalomes import (
    SignalomeModuleSelectionPolicy,
    SignalomeResult,
    build_signalome_result,
)


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
    result = KinaseWorkflow(flank_size=2).run(
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


def _make_analysis_ready_dataset(
    *,
    phospho_matrix: pd.DataFrame,
    protein_ids: dict[str, str],
    metadata_column: str = "protein_id",
) -> AnalysisReadyPhosphoDataset:
    site_index = pd.Index(
        [str(site_id) for site_id in phospho_matrix.index],
        name="site_id",
    )
    aligned_matrix = phospho_matrix.copy(deep=True)
    aligned_matrix.index = site_index
    site_metadata = pd.DataFrame(
        {
            metadata_column: [protein_ids[str(site_id)] for site_id in site_index],
        },
        index=site_index,
    )
    site_sequences = pd.Series(
        ["SEQUENCE"] * len(site_index),
        index=site_index,
        dtype="string",
        name="centralized_sequence",
    )
    provenance = AnalysisReadyPreprocessingProvenance(
        source="signalome workflow test",
        schema=DatasetSchema(),
        comparisons=None,
        row_counts=AnalysisReadyRowCounts(
            total_unique=0,
            total_filtered=0,
            phospho_filtered=0,
            phospho_corrected=0,
            phospho_matrix_sites=len(site_index),
        ),
        site_matrix_stats=AnalysisReadySiteMatrixStats(
            input_rows=len(site_index),
            dropped_missing_sequence=0,
            dropped_incomplete_values=0,
            deduplicated_site_rows=0,
            retained_rows=len(site_index),
        ),
    )
    return AnalysisReadyPhosphoDataset.from_external(
        phospho_matrix=aligned_matrix,
        site_metadata=site_metadata,
        site_sequences=site_sequences,
        phospho_corrected=pd.DataFrame(index=site_index),
        provenance=provenance,
    )


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
        "top_kinase_candidates",
        "top_kinase_weights",
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
    expected_adjacency = kinase_network.copy(deep=True)
    expected_adjacency.index.name = "kinase"
    expected_adjacency.columns.name = "kinase"
    pd.testing.assert_frame_equal(network.adjacency(), expected_adjacency)
    expected_correlation = kinase_correlation_matrix.copy(deep=True)
    expected_correlation.index.name = "kinase"
    expected_correlation.columns.name = "kinase"
    pd.testing.assert_frame_equal(
        network.correlation_matrix,
        expected_correlation,
    )


def test_build_kinase_network_policies_apply_expected_thresholding() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 2.0, 3.0, 4.0],
            "KINASE_B": [4.0, 3.0, 2.0, 1.0],
        }
    )

    positive_only, _ = build_kinase_network(
        scoring_matrix=scoring_matrix,
        threshold=0.9,
        policy="positive_only",
    )
    absolute_threshold, _ = build_kinase_network(
        scoring_matrix=scoring_matrix,
        threshold=0.9,
        policy="absolute_threshold",
    )
    signed, _ = build_kinase_network(
        scoring_matrix=scoring_matrix,
        threshold=0.9,
        policy="signed",
    )

    assert positive_only.loc["KINASE_A", "KINASE_B"] == 0.0
    assert absolute_threshold.loc["KINASE_A", "KINASE_B"] == 1.0
    assert signed.loc["KINASE_A", "KINASE_B"] == -1.0


def test_build_kinase_network_view_includes_negative_edges_for_signed_policy() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 2.0, 3.0, 4.0],
            "KINASE_B": [4.0, 3.0, 2.0, 1.0],
        }
    )
    signed_network, correlation_matrix = build_kinase_network(
        scoring_matrix=scoring_matrix,
        threshold=0.9,
        policy="signed",
    )

    network = build_kinase_network_view(
        kinase_network=signed_network,
        kinase_correlation_matrix=correlation_matrix,
        kinase_substrates={
            "KINASE_A": ("SITE_1",),
            "KINASE_B": ("SITE_2",),
        },
    )

    assert network.neighbor_map == {
        "KINASE_A": ("KINASE_B",),
        "KINASE_B": ("KINASE_A",),
    }
    assert network.edge_table.to_dict("records") == [
        {
            "source_kinase": "KINASE_A",
            "target_kinase": "KINASE_B",
            "correlation": -1.0,
        }
    ]


def test_build_expanded_signalomes_uses_neighbor_map_and_preserves_site_order() -> None:
    site_assignments = pd.DataFrame(
        {
            "protein_id": ["PROTEIN_1", "PROTEIN_2", "PROTEIN_3", "PROTEIN_4"],
            "module_id": [1, 2, 1, 3],
            "top_kinase_candidates": [
                '["KINASE_A"]',
                '["KINASE_B"]',
                '["KINASE_B"]',
                '["KINASE_C"]',
            ],
            "top_kinase_weights": [
                '{"KINASE_A": 1.0}',
                '{"KINASE_B": 1.0}',
                '{"KINASE_B": 1.0}',
                '{"KINASE_C": 1.0}',
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


def test_build_expanded_signalomes_materializes_views_lazily() -> None:
    site_assignments = pd.DataFrame(
        {
            "protein_id": ["PROTEIN_1", "PROTEIN_2", "PROTEIN_3", "PROTEIN_4"],
            "module_id": [1, 2, 1, 3],
            "top_kinase_candidates": [
                '["KINASE_A"]',
                '["KINASE_B"]',
                '["KINASE_B"]',
                '["KINASE_C"]',
            ],
            "top_kinase_weights": [
                '{"KINASE_A": 1.0}',
                '{"KINASE_B": 1.0}',
                '{"KINASE_B": 1.0}',
                '{"KINASE_C": 1.0}',
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
    assert kinase_a._expression_matrix_cache is None
    assert kinase_a._site_assignments_cache is None

    expression_matrix.loc["SITE_1", "sample_1"] = 123.0
    site_assignments.loc["SITE_1", "top_score"] = 0.123
    assert kinase_a.expression_matrix.loc["SITE_1", "sample_1"] == 123.0
    assert kinase_a.site_assignments.loc["SITE_1", "top_score"] == pytest.approx(0.123)


def test_build_expanded_signalomes_reuses_shared_source_tables() -> None:
    site_assignments = pd.DataFrame(
        {
            "protein_id": ["PROTEIN_1", "PROTEIN_2", "PROTEIN_3", "PROTEIN_4"],
            "module_id": [1, 2, 1, 3],
            "top_kinase_candidates": [
                '["KINASE_A"]',
                '["KINASE_B"]',
                '["KINASE_B"]',
                '["KINASE_C"]',
            ],
            "top_kinase_weights": [
                '{"KINASE_A": 1.0}',
                '{"KINASE_B": 1.0}',
                '{"KINASE_B": 1.0}',
                '{"KINASE_C": 1.0}',
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
        kinases_of_interest=["KINASE_A", "KINASE_B"],
        kinase_network={"KINASE_A": ("KINASE_B",), "KINASE_B": ("KINASE_A",)},
        kinase_substrates={
            "KINASE_A": ("SITE_2", "SITE_1"),
            "KINASE_B": ("SITE_3", "SITE_4"),
        },
        signalome_modules=signalome_modules,
        site_assignments=site_assignments,
        expression_matrix=expression_matrix,
        min_kinase_module_share_percent=5.0,
    )

    kinase_a = expanded["KINASE_A"]
    kinase_b = expanded["KINASE_B"]
    assert kinase_a._expression_matrix_source is expression_matrix
    assert kinase_b._expression_matrix_source is expression_matrix
    assert kinase_a._site_assignments_source is site_assignments
    assert kinase_b._site_assignments_source is site_assignments


def test_signalome_workflow_accepts_canonical_pred_mat_result_input() -> None:
    phospho_matrix, pred_mat_workflow_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_workflow_result.scoring_result,
        prediction_result=pred_mat_workflow_result.prediction_result.pred_mat_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_B"],
    )

    assert result.kinases_of_interest == ("KINASE_B",)
    assert result.expanded_signalomes["KINASE_B"].regulated_module_ids == (2,)


def test_signalome_workflow_rejects_pred_mat_without_candidate_kinases() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    empty_pred_mat = PredMatResult(
        pred_mat_result.prediction_result.pred_mat_result.to_frame(copy=True).iloc[
            :, 0:0
        ]
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

    pred_mat = pred_mat_result.prediction_result.pred_mat_result.to_frame(copy=True)
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
    assert result.site_to_protein_resolution_diagnostics is None


def test_signalome_workflow_run_from_analysis_ready_uses_site_metadata_mapping() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    renamed_index = [f"SITE_{i}" for i in range(1, phospho_matrix.shape[0] + 1)]
    renamed_expression_matrix = phospho_matrix.copy()
    renamed_expression_matrix.index = renamed_index

    scoring_result = pred_mat_result.scoring_result
    scoring_result.combined_scores.index = renamed_index
    scoring_result.profile_scores.index = renamed_index

    pred_mat = pred_mat_result.prediction_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = renamed_index
    renamed_pred_mat_result = PredMatResult(pred_mat)
    protein_ids = {
        "SITE_1": "PROTEIN_A",
        "SITE_2": "PROTEIN_A",
        "SITE_3": "PROTEIN_B",
        "SITE_4": "PROTEIN_B",
        "SITE_5": "PROTEIN_C",
        "SITE_6": "PROTEIN_C",
        "SITE_7": "PROTEIN_D",
        "SITE_8": "PROTEIN_D",
    }
    analysis_ready = _make_analysis_ready_dataset(
        phospho_matrix=renamed_expression_matrix,
        protein_ids=protein_ids,
    )

    result = SignalomeWorkflow().run_from_analysis_ready(
        dataset=analysis_ready,
        scoring_result=scoring_result,
        prediction_result=renamed_pred_mat_result,
        kinases_of_interest=["KINASE_A"],
    )

    assert sorted(result.protein_assignments.index.tolist()) == [
        "PROTEIN_A",
        "PROTEIN_B",
        "PROTEIN_C",
        "PROTEIN_D",
    ]
    assert result.site_assignments.loc["SITE_1", "protein_id"] == "PROTEIN_A"
    assert result.site_assignments.loc["SITE_2", "protein_id"] == "PROTEIN_A"
    diagnostics = result.site_to_protein_resolution_diagnostics
    assert diagnostics is not None
    assert diagnostics.fallback_policy == "strict"
    assert diagnostics.chosen_identifier_column == "protein_id"
    assert diagnostics.fallback_mode == "strict_protein_id"
    assert diagnostics.ambiguous_identifier_count == 0


def test_signalome_workflow_run_from_analysis_ready_rejects_non_dataset_input() -> None:
    _, pred_mat_result = _build_pred_mat_workflow_result()

    with pytest.raises(
        TypeError,
        match="dataset must be an AnalysisReadyPhosphoDataset",
    ) as exc_info:
        SignalomeWorkflow().run_from_analysis_ready(
            dataset=object(),  # type: ignore[arg-type]
            scoring_result=pred_mat_result.scoring_result,
            prediction_result=pred_mat_result.prediction_result,
            kinases_of_interest=["KINASE_A"],
        )

    msg = str(exc_info.value)
    assert "builtins.object" in msg
    assert "SimpleKinaseWorkflow.run(...).analysis_ready_dataset" in msg


def test_signalome_workflow_run_from_analysis_ready_strict_mode_rejects_gene_metadata_fallback() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    renamed_index = [f"SITE_{i}" for i in range(1, phospho_matrix.shape[0] + 1)]
    renamed_expression_matrix = phospho_matrix.copy()
    renamed_expression_matrix.index = renamed_index

    scoring_result = pred_mat_result.scoring_result
    scoring_result.combined_scores.index = renamed_index
    scoring_result.profile_scores.index = renamed_index

    pred_mat = pred_mat_result.prediction_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = renamed_index
    renamed_pred_mat_result = PredMatResult(pred_mat)
    protein_ids = {
        f"SITE_{index}": f"PROTEIN_{index}"
        for index in range(1, phospho_matrix.shape[0] + 1)
    }
    analysis_ready = _make_analysis_ready_dataset(
        phospho_matrix=renamed_expression_matrix,
        protein_ids=protein_ids,
        metadata_column="gene",
    )

    with pytest.raises(
        InputCompatibilityError,
        match="required strict site-to-protein column 'protein_id'",
    ):
        SignalomeWorkflow().run_from_analysis_ready(
            dataset=analysis_ready,
            scoring_result=scoring_result,
            prediction_result=renamed_pred_mat_result,
            kinases_of_interest=["KINASE_A"],
        )


def test_signalome_workflow_run_from_analysis_ready_supports_explicit_gene_fallback_mode() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    renamed_index = [f"SITE_{i}" for i in range(1, phospho_matrix.shape[0] + 1)]
    renamed_expression_matrix = phospho_matrix.copy()
    renamed_expression_matrix.index = renamed_index

    scoring_result = pred_mat_result.scoring_result
    scoring_result.combined_scores.index = renamed_index
    scoring_result.profile_scores.index = renamed_index

    pred_mat = pred_mat_result.prediction_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = renamed_index
    renamed_pred_mat_result = PredMatResult(pred_mat)
    gene_mapping = {
        "SITE_1": "PROTEIN_A",
        "SITE_2": "PROTEIN_A",
        "SITE_3": "PROTEIN_B",
        "SITE_4": "PROTEIN_B",
        "SITE_5": "PROTEIN_C",
        "SITE_6": "PROTEIN_C",
        "SITE_7": "PROTEIN_D",
        "SITE_8": "PROTEIN_D",
    }
    analysis_ready = _make_analysis_ready_dataset(
        phospho_matrix=renamed_expression_matrix,
        protein_ids=gene_mapping,
        metadata_column="gene",
    )

    with pytest.warns(
        UserWarning,
        match="Gene-symbol site-to-protein fallback is enabled",
    ):
        result = SignalomeWorkflow().run_from_analysis_ready(
            dataset=analysis_ready,
            scoring_result=scoring_result,
            prediction_result=renamed_pred_mat_result,
            kinases_of_interest=["KINASE_A"],
            metadata_fallback_policy="metadata",
            metadata_protein_columns=["gene"],
            allow_gene_symbol_fallback=True,
        )

    assert sorted(result.protein_assignments.index.tolist()) == [
        "PROTEIN_A",
        "PROTEIN_B",
        "PROTEIN_C",
        "PROTEIN_D",
    ]
    diagnostics = result.site_to_protein_resolution_diagnostics
    assert diagnostics is not None
    assert diagnostics.fallback_policy == "metadata"
    assert diagnostics.chosen_identifier_column == "gene"
    assert diagnostics.fallback_mode == "metadata_gene_symbol"
    assert diagnostics.gene_symbol_fallback_used
    assert diagnostics.ambiguous_identifier_count == 0


def test_signalome_workflow_run_from_analysis_ready_rejects_ambiguous_metadata_fallback() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    ambiguous_gene_mapping = {
        "PROTEIN_1;S1;": "SHARED_GENE",
        "PROTEIN_2;S2;": "SHARED_GENE",
        "PROTEIN_3;S3;": "PROTEIN_3",
        "PROTEIN_4;S4;": "PROTEIN_4",
        "PROTEIN_5;S5;": "PROTEIN_5",
        "PROTEIN_6;S6;": "PROTEIN_6",
        "PROTEIN_7;S7;": "PROTEIN_7",
        "PROTEIN_8;S8;": "PROTEIN_8",
    }
    analysis_ready = _make_analysis_ready_dataset(
        phospho_matrix=phospho_matrix,
        protein_ids=ambiguous_gene_mapping,
        metadata_column="gene",
    )

    with pytest.raises(
        InputCompatibilityError,
        match="Ambiguous site-to-protein metadata mapping detected",
    ):
        SignalomeWorkflow().run_from_analysis_ready(
            dataset=analysis_ready,
            scoring_result=pred_mat_result.scoring_result,
            prediction_result=pred_mat_result.prediction_result,
            kinases_of_interest=["KINASE_A"],
            metadata_fallback_policy="metadata",
            metadata_protein_columns=["gene"],
            allow_gene_symbol_fallback=True,
        )


def test_signalome_workflow_run_from_analysis_ready_reports_ambiguous_metadata_diagnostics_when_allowed() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    ambiguous_gene_mapping = {
        "PROTEIN_1;S1;": "SHARED_GENE",
        "PROTEIN_2;S2;": "SHARED_GENE",
        "PROTEIN_3;S3;": "PROTEIN_3",
        "PROTEIN_4;S4;": "PROTEIN_4",
        "PROTEIN_5;S5;": "PROTEIN_5",
        "PROTEIN_6;S6;": "PROTEIN_6",
        "PROTEIN_7;S7;": "PROTEIN_7",
        "PROTEIN_8;S8;": "PROTEIN_8",
    }
    analysis_ready = _make_analysis_ready_dataset(
        phospho_matrix=phospho_matrix,
        protein_ids=ambiguous_gene_mapping,
        metadata_column="gene",
    )

    with pytest.warns(UserWarning):
        result = SignalomeWorkflow().run_from_analysis_ready(
            dataset=analysis_ready,
            scoring_result=pred_mat_result.scoring_result,
            prediction_result=pred_mat_result.prediction_result,
            kinases_of_interest=["KINASE_A"],
            metadata_fallback_policy="metadata",
            metadata_protein_columns=["gene"],
            allow_gene_symbol_fallback=True,
            allow_ambiguous_metadata_mapping=True,
        )

    diagnostics = result.site_to_protein_resolution_diagnostics
    assert diagnostics is not None
    assert diagnostics.fallback_policy == "metadata"
    assert diagnostics.fallback_mode == "metadata_gene_symbol"
    assert diagnostics.ambiguous_identifier_count == 1
    assert diagnostics.ambiguous_identifiers == ("SHARED_GENE",)
    assert diagnostics.ambiguous_fallback_allowed


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

    pred_mat = pred_mat_result.prediction_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = renamed_index
    renamed_pred_mat_result = PredMatResult(pred_mat)

    with pytest.raises(
        InputCompatibilityError,
        match=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in canonical 'ENTITY;SITE;' "
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


def test_signalome_workflow_rejects_site_to_protein_mapping_with_null_values() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    site_to_protein = {
        f"PROTEIN_{idx};S{idx};": f"PROTEIN_{idx}"
        for idx in range(1, phospho_matrix.shape[0] + 1)
    }
    site_to_protein["PROTEIN_1;S1;"] = None  # type: ignore[assignment]

    with pytest.raises(
        RequestValidationError,
        match="site_to_protein values must be non-null protein IDs",
    ):
        SignalomeWorkflow().run(
            scoring_result=pred_mat_result.scoring_result,
            prediction_result=pred_mat_result.prediction_result,
            expression_matrix=phospho_matrix,
            kinases_of_interest=["KINASE_A"],
            site_to_protein=site_to_protein,
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

    pred_mat = pred_mat_result.prediction_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = malformed_index
    malformed_pred_mat_result = PredMatResult(pred_mat)

    with pytest.raises(
        InputCompatibilityError,
        match=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in canonical 'ENTITY;SITE;' "
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
            "mapping or phosphosite identifiers in canonical 'ENTITY;SITE;' "
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


def test_signalome_workflow_rejects_non_canonical_multi_token_site_ids_without_mapping() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()
    malformed_index = [
        "P0C1X8;AAK1;S677;AEASLSKSKSATTTPSGSPRTSQQNVSNASE",
        *[f"PROTEIN_{i};S{i};" for i in range(2, phospho_matrix.shape[0] + 1)],
    ]
    malformed_expression_matrix = phospho_matrix.copy()
    malformed_expression_matrix.index = malformed_index

    scoring_result = pred_mat_result.scoring_result
    scoring_result.combined_scores.index = malformed_index
    scoring_result.profile_scores.index = malformed_index

    pred_mat = pred_mat_result.prediction_result.pred_mat_result.to_frame(copy=True)
    pred_mat.index = malformed_index
    malformed_pred_mat_result = PredMatResult(pred_mat)

    with pytest.raises(
        InputCompatibilityError,
        match=(
            "Signalome construction requires either an explicit site_to_protein "
            "mapping or phosphosite identifiers in canonical 'ENTITY;SITE;' "
            "format"
        ),
    ):
        SignalomeWorkflow().run(
            scoring_result=scoring_result,
            prediction_result=malformed_pred_mat_result,
            expression_matrix=malformed_expression_matrix,
            kinases_of_interest=["KINASE_A"],
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


def test_build_site_assignments_tracks_tied_top_kinases_as_weighted_multi_assignments() -> (
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

    assert tied_row["top_kinase_candidates"] == ("KINASE_A", "KINASE_B")
    assert tied_row["top_kinase_weights"] == (("KINASE_A", 0.5), ("KINASE_B", 0.5))
    assert tied_row["top_kinase_tie_count"] == 2
    assert bool(tied_row["top_kinase_is_ambiguous"])

    assert clear_row["top_kinase_candidates"] == ("KINASE_A",)
    assert clear_row["top_kinase_weights"] == (("KINASE_A", 1.0),)
    assert clear_row["top_kinase_tie_count"] == 1
    assert not bool(clear_row["top_kinase_is_ambiguous"])


def test_weighted_top_assignment_policy_propagates_fractional_module_shares() -> None:
    site_ids = ["PROTEIN_1;S1;", "PROTEIN_1;S2;", "PROTEIN_2;S3;", "PROTEIN_2;S4;"]
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0, 4.0, 4.0],
            "KINASE_B": [1.1, 1.1, 4.1, 4.1],
        },
        index=site_ids,
    )
    pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.95, 0.10, 0.10],
            "KINASE_B": [0.95, 0.95, 0.91, 0.91],
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

    weighted = build_signalome_result(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A"],
        assignment_policy="weighted_top",
        signalome_cutoff=0.99,
        module_count=2,
    )
    binary = build_signalome_result(
        scoring_matrix=scoring_matrix,
        pred_mat=pred_mat,
        expression_matrix=expression_matrix,
        kinases_of_interest=["KINASE_A"],
        assignment_policy="cutoff_binary",
        signalome_cutoff=0.99,
        module_count=2,
    )

    assert weighted.signalome_modules.loc[1, "KINASE_A"] == 50.0
    assert weighted.signalome_modules.loc[1, "KINASE_B"] == 50.0
    assert weighted.signalome_modules.loc[2, "KINASE_B"] == 100.0
    assert binary.signalome_modules.loc[1, "KINASE_A"] == 0.0
    assert binary.signalome_modules.loc[1, "KINASE_B"] == 0.0


def test_build_signalome_support_matrix_supports_weighted_top_policy() -> None:
    site_assignments = pd.DataFrame(
        {
            "protein_id": ["PROTEIN_1", "PROTEIN_1", "PROTEIN_2"],
            "module_id": [1, 1, 2],
            "top_kinase_candidates": [
                ("KINASE_A", "KINASE_B"),
                ("KINASE_A",),
                ("KINASE_B",),
            ],
            "top_kinase_weights": [
                (("KINASE_A", 0.5), ("KINASE_B", 0.5)),
                (("KINASE_A", 1.0),),
                (("KINASE_B", 1.0),),
            ],
            "top_kinase_tie_count": [2, 1, 1],
            "top_kinase_is_ambiguous": [True, False, False],
            "top_score": [0.95, 0.96, 0.92],
        },
        index=pd.Index(["SITE_1", "SITE_2", "SITE_3"], name="site_id"),
    )

    support = build_signalome_support_matrix(
        site_assignments=site_assignments,
        kinase_substrates={
            "KINASE_A": ("SITE_2",),
            "KINASE_B": ("SITE_1", "SITE_3"),
        },
        kinases_of_interest=["KINASE_A", "KINASE_B"],
        assignment_policy="weighted_top",
    )

    assert support.loc["KINASE_A"].tolist() == [0.5, 1.0, 0.0]
    assert support.loc["KINASE_B"].tolist() == [0.5, 0.0, 1.0]


def test_build_site_assignments_rejects_missing_site_to_protein_mapping() -> None:
    pred_mat = pd.DataFrame({"KINASE_A": [0.9]}, index=["SITE_1"])
    protein_modules = pd.Series({"PROTEIN_1": 1}, name="module_id")
    site_to_protein = pd.Series({"SITE_2": "PROTEIN_1"})

    with pytest.raises(
        InputCompatibilityError,
        match="site_to_protein must define a protein ID for every pred_mat site",
    ):
        build_site_assignments(
            pred_mat=pred_mat,
            protein_modules=protein_modules,
            site_to_protein=site_to_protein,
        )


def test_build_site_assignments_rejects_rows_without_top_kinase_assignment() -> None:
    pred_mat = pd.DataFrame(
        {"KINASE_A": [float("nan")], "KINASE_B": [float("nan")]},
        index=["SITE_1"],
    )
    protein_modules = pd.Series({"PROTEIN_1": 1}, name="module_id")
    site_to_protein = pd.Series({"SITE_1": "PROTEIN_1"})

    with pytest.raises(
        InputCompatibilityError,
        match=(
            "pred_mat contains phosphosite rows with no top-kinase assignment "
            r"\(all scores missing\). Offending site IDs: SITE_1"
        ),
    ):
        build_site_assignments(
            pred_mat=pred_mat,
            protein_modules=protein_modules,
            site_to_protein=site_to_protein,
        )


def test_extract_top_kinase_assignments_tracks_ties_directly() -> None:
    sorted_pred_mat = pd.DataFrame(
        {
            "KINASE_A": [0.8, 0.9],
            "KINASE_B": [0.8, 0.2],
        },
        index=["SITE_1", "SITE_2"],
    )

    assignments = signalome_assignments._extract_top_kinase_assignments(
        sorted_pred_mat=sorted_pred_mat,
        site_index=pd.Index(sorted_pred_mat.index, name="site_id"),
    )

    assert assignments.top_kinase_candidates == [
        ("KINASE_A", "KINASE_B"),
        ("KINASE_A",),
    ]
    assert assignments.top_kinase_weights == [
        (("KINASE_A", 0.5), ("KINASE_B", 0.5)),
        (("KINASE_A", 1.0),),
    ]
    assert assignments.top_kinase_tie_count.tolist() == [2, 1]


def test_resolve_pre_scoring_module_selection_clamps_explicit_module_count() -> None:
    scoring_values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [1.1, 2.1, 3.1],
            [0.9, 1.9, 2.9],
        ],
        dtype=float,
    )
    profile_degeneracy = signalome_clustering.summarize_profile_degeneracy(
        scoring_values
    )
    policy = signalome_clustering.SignalomeModuleSelectionPolicy()

    selection, _ = signalome_clustering._resolve_pre_scoring_module_selection(
        resolved_policy=policy,
        requested_module_count=5,
        n_sites=scoring_values.shape[0],
        profile_degeneracy=profile_degeneracy,
        correlation_exclusion_note=signalome_clustering.build_correlation_exclusion_note(
            profile_degeneracy
        ),
    )

    assert selection is not None
    assert selection.diagnostics.selected_module_count == 3
    assert selection.diagnostics.requested_module_count == 5
    assert (
        selection.diagnostics.reason
        == "module_count was provided explicitly by the caller"
    )


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


def test_select_module_count_avoids_full_correlation_matrix_for_large_inputs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(signalome_clustering, "MAX_FULL_CORRELATION_SITE_COUNT", 20)
    monkeypatch.setattr(
        signalome_clustering,
        "MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER",
        5,
    )
    scoring_values = pd.DataFrame(
        {
            "KINASE_A": [float(index) / 10.0 for index in range(40)],
            "KINASE_B": [float((index * 7 + 3) % 97) / 11.0 for index in range(40)],
            "KINASE_C": [
                float((index * index + 5) % 101) / 13.0 for index in range(40)
            ],
        }
    ).to_numpy(dtype=float)
    original_corrcoef = signalome_clustering.np.corrcoef
    observed_rows: list[int] = []

    def guarded_corrcoef(values: object) -> object:
        shape = getattr(values, "shape", None)
        if shape is not None and len(shape) >= 1:
            observed_rows.append(int(shape[0]))
            if int(shape[0]) == scoring_values.shape[0]:
                msg = "full correlation matrix materialization should be skipped"
                raise AssertionError(msg)
        return original_corrcoef(values)

    monkeypatch.setattr(signalome_clustering.np, "corrcoef", guarded_corrcoef)

    diagnostics = signalome_clustering.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
    )

    assert diagnostics.selected_module_count >= 1
    assert "sampled within-cluster correlation estimates" in diagnostics.reason
    assert observed_rows
    assert max(observed_rows) < scoring_values.shape[0]


def test_cluster_median_correlation_approximate_is_reproducible_and_row_order_invariant_for_large_clusters() -> (
    None
):
    random_generator = np.random.default_rng(2026)
    cluster_a = random_generator.normal(loc=0.0, scale=0.25, size=(180, 6))
    cluster_b = random_generator.normal(loc=3.5, scale=0.25, size=(180, 6))
    scoring_values = np.vstack([cluster_a, cluster_b]).astype(float, copy=False)
    labels = np.concatenate(
        [
            np.zeros(180, dtype=int),
            np.ones(180, dtype=int),
        ]
    )

    baseline = signalome_clustering.cluster_median_correlation_approximate(
        scoring_values=scoring_values,
        labels=labels,
        label=0,
        max_sites_per_cluster=64,
    )
    repeat = signalome_clustering.cluster_median_correlation_approximate(
        scoring_values=scoring_values,
        labels=labels,
        label=0,
        max_sites_per_cluster=64,
    )

    permutation = np.random.default_rng(7).permutation(scoring_values.shape[0])
    reordered = signalome_clustering.cluster_median_correlation_approximate(
        scoring_values=scoring_values[permutation],
        labels=labels[permutation],
        label=0,
        max_sites_per_cluster=64,
    )

    assert repeat == pytest.approx(baseline, abs=0.0, rel=0.0)
    assert reordered == pytest.approx(baseline, abs=0.0, rel=0.0)


def test_select_module_count_approximate_path_is_row_order_invariant(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(signalome_clustering, "MAX_FULL_CORRELATION_SITE_COUNT", 20)
    monkeypatch.setattr(
        signalome_clustering,
        "MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER",
        32,
    )
    random_generator = np.random.default_rng(42)
    cluster_a = random_generator.normal(loc=0.0, scale=0.3, size=(120, 5))
    cluster_b = random_generator.normal(loc=4.0, scale=0.3, size=(120, 5))
    scoring_values = np.vstack([cluster_a, cluster_b]).astype(float, copy=False)
    permutation = np.random.default_rng(99).permutation(scoring_values.shape[0])

    original = signalome_clustering.select_module_count_with_diagnostics(
        scoring_values=scoring_values
    )
    reordered = signalome_clustering.select_module_count_with_diagnostics(
        scoring_values=scoring_values[permutation]
    )

    assert original.selected_module_count == reordered.selected_module_count
    assert set(original.candidate_scores) == set(reordered.candidate_scores)
    for cluster_count in original.candidate_scores:
        assert original.candidate_scores[cluster_count].min_median_correlation == (
            pytest.approx(
                reordered.candidate_scores[cluster_count].min_median_correlation,
                abs=1e-12,
                rel=0.0,
            )
        )
        assert original.candidate_scores[cluster_count].mean_median_correlation == (
            pytest.approx(
                reordered.candidate_scores[cluster_count].mean_median_correlation,
                abs=1e-12,
                rel=0.0,
            )
        )


def test_select_module_count_excludes_constant_profiles_without_runtime_warnings() -> (
    None
):
    scoring_values = pd.DataFrame(
        {
            "KINASE_A": [1.0, 1.0, 3.0, 4.0],
            "KINASE_B": [1.0, 2.0, 2.0, 5.0],
            "KINASE_C": [1.0, 3.0, 1.0, 6.0],
        }
    ).to_numpy(dtype=float)

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        diagnostics = signalome_clustering.select_module_count_with_diagnostics(
            scoring_values=scoring_values,
        )

    assert diagnostics.selected_module_count >= 1
    assert diagnostics.zero_variance_profile_count == 1
    assert diagnostics.near_constant_profile_count == 0
    assert diagnostics.excluded_from_correlation_count == 1
    assert (
        "Excluded 1 degenerate profile from correlation scoring" in diagnostics.reason
    )


def test_select_module_count_reports_near_constant_profiles_in_diagnostics() -> None:
    scoring_values = pd.DataFrame(
        {
            "KINASE_A": [1.0, 0.0, 2.0, 3.0],
            "KINASE_B": [1.0 + 1e-8, 1.0, 2.5, 2.0],
            "KINASE_C": [1.0 - 1e-8, 2.0, 3.0, 1.0],
        }
    ).to_numpy(dtype=float)

    diagnostics = signalome_clustering.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
    )

    assert diagnostics.selected_module_count >= 1
    assert diagnostics.zero_variance_profile_count == 0
    assert diagnostics.near_constant_profile_count == 1
    assert diagnostics.excluded_from_correlation_count == 1
    assert "near-constant" in diagnostics.reason


def test_select_module_count_falls_back_when_too_few_non_degenerate_profiles_remain() -> (
    None
):
    scoring_values = pd.DataFrame(
        {
            "KINASE_A": [1.0, 2.0, 3.0],
            "KINASE_B": [1.0, 2.0, 3.0 + 1e-8],
            "KINASE_C": [1.0, 2.0, 3.0 - 1e-8],
        }
    ).to_numpy(dtype=float)

    diagnostics = signalome_clustering.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
    )

    assert diagnostics.selected_module_count == 1
    assert diagnostics.zero_variance_profile_count == 2
    assert diagnostics.near_constant_profile_count == 1
    assert diagnostics.excluded_from_correlation_count == 3
    assert diagnostics.candidate_scores == {}
    assert "fewer than two non-degenerate phosphosite profiles remained" in (
        diagnostics.reason
    )


def test_cluster_sites_reuses_cached_candidate_labels_for_automatic_selection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 2.0, 3.0, -1.0, -2.0, -3.0],
            "KINASE_B": [2.0, 4.0, 6.0, -2.0, -4.0, -6.0],
            "KINASE_C": [3.0, 6.0, 9.0, -3.0, -6.0, -9.0],
        },
        index=[f"PROTEIN_{i};S{i};" for i in range(1, 7)],
    )

    observed_fit_calls: list[int] = []
    original_fit_cluster_labels = signalome_clustering.fit_cluster_labels

    def counting_fit_cluster_labels(
        scoring_values: object,
        cluster_count: int,
    ) -> object:
        observed_fit_calls.append(int(cluster_count))
        return original_fit_cluster_labels(scoring_values, cluster_count)

    monkeypatch.setattr(
        signalome_clustering,
        "fit_cluster_labels",
        counting_fit_cluster_labels,
    )

    result = signalome_clustering.cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=None,
    )

    assert result.module_selection_diagnostics.selected_module_count == 2
    assert observed_fit_calls == []


def test_cluster_sites_matches_legacy_two_pass_partition_assignments() -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 2.0, 3.0, -1.0, -2.0, -3.0],
            "KINASE_B": [2.0, 4.0, 6.0, -2.0, -4.0, -6.0],
            "KINASE_C": [3.0, 6.0, 9.0, -3.0, -6.0, -9.0],
        },
        index=[f"PROTEIN_{i};S{i};" for i in range(1, 7)],
    )
    scoring_values = scoring_matrix.to_numpy(dtype=float)
    legacy_diagnostics = signalome_clustering.select_module_count_with_diagnostics(
        scoring_values=scoring_values,
        requested_module_count=None,
    )
    legacy_module_count = max(
        1,
        min(legacy_diagnostics.selected_module_count, scoring_values.shape[0]),
    )
    if legacy_module_count == 1:
        legacy_labels = np.ones(scoring_values.shape[0], dtype=int)
    else:
        legacy_labels = (
            signalome_clustering.fit_cluster_labels(
                scoring_values,
                legacy_module_count,
            )
            + 1
        )

    result = signalome_clustering.cluster_sites_with_diagnostics(
        scoring_matrix=scoring_matrix,
        requested_module_count=None,
    )
    result_labels = result.site_clusters.to_numpy(dtype=int, copy=False)
    legacy_partition = legacy_labels[:, None] == legacy_labels[None, :]
    result_partition = result_labels[:, None] == result_labels[None, :]

    assert result.module_selection_diagnostics == legacy_diagnostics
    assert np.array_equal(result_partition, legacy_partition)


def test_build_kinase_network_zeros_diagonal_with_numpy_fill_diagonal(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    scoring_matrix = pd.DataFrame(
        {
            "KINASE_A": [1.0, 2.0, 3.0, 4.0],
            "KINASE_B": [1.1, 2.1, 3.1, 4.1],
        }
    )
    original_fill_diagonal = signalome_analysis.np.fill_diagonal
    observed_calls: list[float] = []

    def counting_fill_diagonal(array: object, value: float) -> None:
        observed_calls.append(float(value))
        original_fill_diagonal(array, value)

    monkeypatch.setattr(
        signalome_analysis.np,
        "fill_diagonal",
        counting_fill_diagonal,
    )

    network, _ = build_kinase_network(
        scoring_matrix=scoring_matrix,
        threshold=0.0,
        policy="signed",
    )

    assert observed_calls == [0.0]
    assert float(network.loc["KINASE_A", "KINASE_A"]) == 0.0
    assert float(network.loc["KINASE_B", "KINASE_B"]) == 0.0


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

    assert result.modules.to_frame() is not result.signalome_modules
    pd.testing.assert_frame_equal(
        result.modules.to_frame(),
        result.signalome_modules,
    )
    assert list(result.modules.to_relationship_table().columns) == [
        "module_id",
        "kinase",
        "share_percent",
    ]
    assert result.modules.to_relationship_table().to_dict("records") == [
        {"module_id": 1, "kinase": "KINASE_A", "share_percent": 100.0},
        {"module_id": 2, "kinase": "KINASE_B", "share_percent": 100.0},
    ]

    assert result.assignments.sites() is not result.site_assignments
    pd.testing.assert_frame_equal(
        result.assignments.sites(),
        result.site_assignments,
    )
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
        "kinase_adjacency_matrix",
        "kinase_correlation_matrix",
    ]
    assert "scoring_matrix" not in frames
    assert frames["signalome_modules"].equals(result.signalome_modules)
    assert frames["protein_assignments"].equals(result.protein_assignments)
    assert frames["kinase_adjacency_matrix"].equals(result.kinase_adjacency_matrix)
    assert frames["kinase_correlation_matrix"].equals(result.kinase_correlation_matrix)

    frames_with_inputs = result.to_frames(include_inputs=True)
    assert list(frames_with_inputs)[-3:] == [
        "scoring_matrix",
        "pred_mat",
        "expression_matrix",
    ]


def test_signalome_result_expanded_signalomes_materialize_with_parity() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    mutable = result.expanded_signalomes_mutable_unsafe()["KINASE_A"]
    assert mutable._expression_matrix_cache is None
    assert mutable._site_assignments_cache is None

    detached = result.expanded_signalomes["KINASE_A"]
    pd.testing.assert_frame_equal(detached.expression_matrix, mutable.expression_matrix)
    pd.testing.assert_frame_equal(detached.site_assignments, mutable.site_assignments)


def test_signalome_result_expanded_signalomes_are_detached_from_mutable_state() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    mutable = result.expanded_signalomes_mutable_unsafe()["KINASE_A"]
    detached = result.expanded_signalomes["KINASE_A"]

    mutable_expression = mutable.expression_matrix
    mutable_assignments = mutable.site_assignments
    expression_original = float(mutable_expression.iloc[0, 0])
    top_score_col = mutable_assignments.columns.get_loc("top_score")
    assignment_original = float(mutable_assignments.iloc[0, top_score_col])

    detached.expression_matrix.iloc[0, 0] = expression_original + 100.0
    detached.site_assignments.iloc[0, top_score_col] = assignment_original + 0.1

    assert float(mutable.expression_matrix.iloc[0, 0]) == expression_original
    assert float(mutable.site_assignments.iloc[0, top_score_col]) == assignment_original


def test_signalome_result_wrappers_with_pandas_state_are_not_frozen_dataclasses() -> (
    None
):
    from phospy.signalomes import (
        ExpandedSignalome,
        SignalomeAssignments,
        SignalomeKinaseNetwork,
        SignalomeModules,
    )

    assert ExpandedSignalome.__dataclass_params__.frozen is False
    assert SignalomeModules.__dataclass_params__.frozen is False
    assert SignalomeAssignments.__dataclass_params__.frozen is False
    assert SignalomeKinaseNetwork.__dataclass_params__.frozen is False


def test_signalome_result_to_frames_supports_safe_copy_default_and_explicit_mutable_unsafe() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    detached_frames = result.to_frames()
    mutable_frames = result.to_mutable_frames_unsafe()

    assert (
        detached_frames["signalome_modules"] is not mutable_frames["signalome_modules"]
    )
    assert detached_frames["site_assignments"] is not mutable_frames["site_assignments"]


def test_signalome_result_to_owned_frames_returns_shared_state() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    owned_frames = result.to_owned_frames()
    mutable_frames = result.to_mutable_frames_unsafe()

    assert owned_frames["signalome_modules"] is mutable_frames["signalome_modules"]
    assert owned_frames["site_assignments"] is mutable_frames["site_assignments"]


def test_signalome_nested_result_ownership_accessors_are_explicit() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    modules = result.modules
    detached_modules = modules.to_frame()
    owned_modules = modules.to_owned_frame()
    mutable_modules = modules.to_mutable_frame_unsafe()
    assert detached_modules is not owned_modules
    assert owned_modules is mutable_modules

    assignments = result.assignments
    detached_sites = assignments.to_site_assignments()
    owned_sites = assignments.to_owned_site_assignments()
    mutable_sites = assignments.to_mutable_site_assignments_unsafe()
    assert detached_sites is not owned_sites
    assert owned_sites is mutable_sites


def test_signalome_result_default_accessors_are_detached_from_owned_state() -> None:
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    mutable_frames = result.to_mutable_frames_unsafe()

    module_copy = result.signalome_modules
    module_original = float(mutable_frames["signalome_modules"].iloc[0, 0])
    module_copy.iloc[0, 0] = module_original + 99.0
    assert float(mutable_frames["signalome_modules"].iloc[0, 0]) == module_original

    site_copy = result.site_assignments
    module_col = site_copy.columns.get_loc("module_id")
    original_site_module = int(mutable_frames["site_assignments"].iloc[0, module_col])
    site_copy.iloc[0, module_col] = original_site_module + 99
    assert (
        int(mutable_frames["site_assignments"].iloc[0, module_col])
        == original_site_module
    )


def test_signalome_result_to_mutable_frames_unsafe_exposes_explicit_mutable_state() -> (
    None
):
    phospho_matrix, pred_mat_result = _build_pred_mat_workflow_result()

    result = SignalomeWorkflow().run(
        scoring_result=pred_mat_result.scoring_result,
        prediction_result=pred_mat_result.prediction_result,
        expression_matrix=phospho_matrix,
        kinases_of_interest=["KINASE_A"],
    )

    mutable_frames = result.to_mutable_frames_unsafe()
    mutable_frames["signalome_modules"].iloc[0, 0] = 321.0
    assert (
        float(result.to_mutable_frames_unsafe()["signalome_modules"].iloc[0, 0])
        == 321.0
    )


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
        "kinase_adjacency_matrix",
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
    pd.testing.assert_frame_equal(
        reloaded_site_assignments,
        serialize_site_assignments_for_export(result.site_assignments),
    )
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
