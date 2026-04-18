from __future__ import annotations

import json
from pathlib import Path

import pytest

from phospy import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    ReferencePreset,
    SignalomeConfig,
    SignalomeWorkflow,
    SignalomeWorkflowRequest,
)
from phospy.io.publishing import publish_signalome_workflow
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
    assert written["signalome.manifest"] == output_root / "signalome" / "manifest.json"
    assert "signalome.expanded_signalome" not in written

    assert (output_root / "dataset" / "manifest.json").exists()
    assert (output_root / "kinase" / "manifest.json").exists()
    assert (output_root / "signalome" / "module_assignments.csv").exists()
    assert (output_root / "signalome" / "signalome_modules.csv").exists()
    assert (output_root / "signalome" / "kinase_network_nodes.csv").exists()
    assert (output_root / "signalome" / "kinase_network_edges.csv").exists()
    assert not (output_root / "signalome" / "expanded_signalome.csv").exists()

    manifest = json.loads(
        (output_root / "signalome" / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest == {
        "expanded_signalome_present": False,
        "kinase_network_nodes_present": True,
        "output_format": "csv",
        "reference_organism": "rat",
    }


def _build_signalome_result():
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=1),
            prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
            activity_config=None,
        )
    )
    return SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(signalome_cutoff=0.5),
        )
    )
