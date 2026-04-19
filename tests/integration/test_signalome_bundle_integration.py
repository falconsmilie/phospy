from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
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
from phospy.io.bundles.signalome import (
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
    SignalomeWorkflowConfigSnapshot,
    load_signalome_workflow_bundle,
    save_signalome_workflow_bundle,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset

pytestmark = pytest.mark.integration


def test_signalome_bundle_round_trip_preserves_outputs_and_config(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    config_snapshot = SignalomeWorkflowConfigSnapshot.from_request(request)
    bundle_root = tmp_path / "signalome_bundle"

    written = save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=config_snapshot,
        output_format="csv",
    )

    assert written["manifest"] == bundle_root / "manifest.json"
    loaded = load_signalome_workflow_bundle(bundle_root)

    assert loaded.manifest_version == SIGNALOME_BUNDLE_MANIFEST_VERSION
    assert loaded.config_snapshot == config_snapshot
    _assert_signalome_result_equal(loaded.result, result)


def test_signalome_bundle_manifest_v1_is_explicit_and_handles_optional_outputs(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / "signalome_bundle"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))

    assert manifest["manifest_version"] == SIGNALOME_BUNDLE_MANIFEST_VERSION
    assert manifest["bundle_type"] == "signalome_workflow_result"
    assert manifest["table_format"] == "csv"
    assert manifest["config_snapshot"] == "config/snapshot.json"
    assert manifest["upstream_kinase_outputs"]["activity"]["enabled"] is False
    assert manifest["upstream_kinase_outputs"]["scoring"]["tables"] == {
        "combined_scores": "scoring/combined_scores.csv",
        "motif_scores": None,
        "profile_scores": "scoring/profile_scores.csv",
        "weights": None,
    }
    assert manifest["signalome_outputs"]["tables"] == {
        "expanded_signalome": None,
        "kinase_network_edges": "signalome/kinase_network_edges.csv",
        "kinase_network_nodes": "signalome/kinase_network_nodes.csv",
        "module_assignments": "signalome/module_assignments.csv",
        "signalome_modules": "signalome/signalome_modules.csv",
    }
    assert manifest["signalome_outputs"]["metadata"] == {
        "expanded_signalome_present": False,
        "kinase_network_nodes_present": True,
    }

    loaded = load_signalome_workflow_bundle(bundle_root)
    assert loaded.result.expanded_signalome is None


def test_signalome_config_snapshot_accepts_legacy_cutoff_payload() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        {"signalome_config": {"signalome_cutoff": 0.6}}
    )
    assert snapshot.signalome_config.substrate_support_cutoff == pytest.approx(0.6)
    assert snapshot.signalome_config.network_correlation_threshold == pytest.approx(0.6)


def _build_signalome_request_and_result():
    dataset = build_rat_l6_dataset(n_sites=260)
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=6, ensemble_size=12),
            activity_config=None,
        )
    )
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )
    result = SignalomeWorkflow().run(request)
    return request, result


def _assert_signalome_result_equal(left, right) -> None:
    assert left.dataset.organism == right.dataset.organism
    assert left.dataset.transformation_state == right.dataset.transformation_state
    assert (
        left.kinase_result.references.organism
        == right.kinase_result.references.organism
    )

    pd.testing.assert_frame_equal(
        left.dataset.phospho,
        right.dataset.phospho,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.dataset.site_metadata,
        right.dataset.site_metadata,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.kinase_result.references.kinase_substrate_map,
        right.kinase_result.references.kinase_substrate_map,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.kinase_result.references.site_sequences,
        right.kinase_result.references.site_sequences,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.kinase_result.scoring_result.profile_scores,
        right.kinase_result.scoring_result.profile_scores,
        check_dtype=False,
        check_names=False,
    )
    _assert_optional_frame_equal(
        left.kinase_result.scoring_result.motif_scores,
        right.kinase_result.scoring_result.motif_scores,
    )
    _assert_optional_frame_equal(
        left.kinase_result.scoring_result.combined_scores,
        right.kinase_result.scoring_result.combined_scores,
    )
    _assert_optional_frame_equal(
        left.kinase_result.scoring_result.weights,
        right.kinase_result.scoring_result.weights,
    )
    pd.testing.assert_frame_equal(
        left.kinase_result.prediction_result.pred_mat,
        right.kinase_result.prediction_result.pred_mat,
        check_dtype=False,
        check_names=False,
    )
    assert left.kinase_result.activity_result is None
    assert right.kinase_result.activity_result is None

    pd.testing.assert_frame_equal(
        left.module_assignments.table,
        right.module_assignments.table,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.signalome_modules.table,
        right.signalome_modules.table,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.kinase_network.edges,
        right.kinase_network.edges,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.kinase_network.nodes,
        right.kinase_network.nodes,
        check_dtype=False,
        check_names=False,
    )
    assert left.expanded_signalome is None
    assert right.expanded_signalome is None


def _assert_optional_frame_equal(left, right) -> None:
    if left is None or right is None:
        assert left is right
        return
    pd.testing.assert_frame_equal(
        left,
        right,
        check_dtype=False,
        check_names=False,
    )
