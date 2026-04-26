from __future__ import annotations

import json
from pathlib import Path

import pytest

from phospy import (
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflowRequest,
)
from phospy.api.results import SignalomeWorkflowResult
from phospy.io.publishers.workflows import publish_signalome_workflow
from tests.support.rewrite_fixture_data import build_rat_l6_dataset

pytestmark = pytest.mark.integration


def test_publish_signalome_workflow_writes_supported_lane_output_layout(
    tmp_path: Path,
) -> None:
    result = _build_signalome_result()
    output_root = tmp_path / "out"

    written = publish_signalome_workflow(result, output_root, output_format="csv")

    assert written["dataset.phospho"] == output_root / "dataset" / "phospho.csv"
    assert (
        written["kinase.scoring.profile_scores"]
        == output_root / "kinase" / "scoring" / "profile_scores.csv"
    )
    assert (
        written["kinase.scoring.rank_weighted_fusion_scores"]
        == output_root / "kinase" / "scoring" / "rank_weighted_fusion_scores.csv"
    )
    assert "kinase.scoring.motif_scores" not in written
    assert "kinase.scoring.score_fusion_weights" not in written
    assert (
        written["signalome.module_assignments"]
        == output_root / "signalome" / "module_assignments.csv"
    )
    assert (
        written["signalome.signalome_modules"]
        == output_root / "signalome" / "signalome_modules.csv"
    )
    assert (
        written["signalome.kinase_network.nodes"]
        == output_root / "signalome" / "kinase_network_nodes.csv"
    )
    assert (
        written["signalome.kinase_network.edges"]
        == output_root / "signalome" / "kinase_network_edges.csv"
    )
    assert (
        written["signalome.kinase_network.candidate_correlations"]
        == output_root / "signalome" / "kinase_network_candidate_correlations.csv"
    )
    assert written["signalome.manifest"] == output_root / "signalome" / "manifest.json"
    assert (
        written["signalome.expanded_signalome"]
        == output_root / "signalome" / "expanded_signalome.csv"
    )

    assert (output_root / "dataset" / "manifest.json").exists()
    assert (output_root / "kinase" / "manifest.json").exists()
    assert (output_root / "signalome" / "module_assignments.csv").exists()
    assert (output_root / "signalome" / "signalome_modules.csv").exists()
    assert (output_root / "signalome" / "kinase_network_nodes.csv").exists()
    assert (output_root / "signalome" / "kinase_network_edges.csv").exists()
    assert (
        output_root / "signalome" / "kinase_network_candidate_correlations.csv"
    ).exists()
    assert (output_root / "signalome" / "expanded_signalome.csv").exists()

    manifest = json.loads(
        (output_root / "signalome" / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["expanded_signalome_present"] is True
    assert manifest["kinase_network_nodes_present"] is True
    assert manifest["output_format"] == "csv"
    assert manifest["reference_organism"] == "rat"
    assert manifest["module_selection_strategy"] in {
        "correlation_thresholds",
        "explicit_module_count",
    }
    assert int(manifest["selected_module_count"]) >= 1
    assert manifest["used_automatic_module_selection"] is True
    preconditioning = manifest["score_preconditioning_diagnostics"]
    assert preconditioning["policy"] == "allow_and_report"
    assert int(preconditioning["input_row_count"]) >= 0
    assert int(preconditioning["dropped_all_missing_row_count"]) >= 0
    assert int(preconditioning["retained_row_count"]) >= 0
    network_correlation_diagnostics = manifest["network_correlation_diagnostics"]
    assert int(network_correlation_diagnostics["total_candidate_correlations"]) >= 0
    assert int(network_correlation_diagnostics["finite_correlations"]) >= 0
    assert int(network_correlation_diagnostics["undefined_correlations"]) >= 0
    assert int(network_correlation_diagnostics["edges_created"]) >= 0
    assert manifest["provenance"]["workflow_name"] == "signalome_workflow"

    dataset_manifest = json.loads(
        (output_root / "dataset" / "manifest.json").read_text(encoding="utf-8")
    )
    assert dataset_manifest["provenance"]["workflow_name"] == "dataset_builder"

    kinase_manifest = json.loads(
        (output_root / "kinase" / "manifest.json").read_text(encoding="utf-8")
    )
    assert kinase_manifest["provenance"]["workflow_name"] == "kinase_workflow"


def _build_signalome_result():
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=6,
                deterministic_max_selected_kinases=12,
                adaptive_ensemble_runs=12,
            ),
            activity_config=None,
        )
    )
    return SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(substrate_support_cutoff=0.5),
        )
    )


def test_publish_signalome_workflow_omits_expanded_when_result_field_is_none(
    tmp_path: Path,
) -> None:
    result = _build_signalome_result()
    without_expanded = SignalomeWorkflowResult._from_owned(
        dataset=result.dataset,
        kinase_result=result.kinase_result,
        module_assignments=result.module_assignments,
        signalome_modules=result.signalome_modules,
        kinase_network=result.kinase_network,
        module_selection_diagnostics=result.module_selection_diagnostics,
        expanded_signalome=None,
    )

    written = publish_signalome_workflow(
        without_expanded,
        tmp_path / "out",
        output_format="csv",
    )

    assert "signalome.expanded_signalome" not in written
