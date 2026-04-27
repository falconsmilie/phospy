from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
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
from phospy.api.configs import SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
from phospy.api.results import SignalomeWorkflowResult
from phospy.io.bundles.signalome import (
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
    SignalomeWorkflowConfigSnapshot,
    load_signalome_workflow_bundle,
    save_signalome_workflow_bundle,
)
from phospy.provenance.serialization import to_payload as provenance_to_payload
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
    assert loaded.result.provenance is not None
    assert result.provenance is not None
    assert provenance_to_payload(loaded.result.provenance) == provenance_to_payload(
        result.provenance
    )
    assert loaded.result.kinase_result.provenance is not None
    assert result.kinase_result.provenance is not None
    assert provenance_to_payload(
        loaded.result.kinase_result.provenance
    ) == provenance_to_payload(result.kinase_result.provenance)
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
        "rank_weighted_fusion_scores": "scoring/rank_weighted_fusion_scores.csv",
        "motif_scores": None,
        "profile_scores": "scoring/profile_scores.csv",
        "score_fusion_weights": None,
    }
    assert manifest["signalome_outputs"]["tables"] == {
        "expanded_signalome": "signalome/expanded_signalome.csv",
        "kinase_network_candidate_correlations": "signalome/kinase_network_candidate_correlations.csv",
        "kinase_network_edges": "signalome/kinase_network_edges.csv",
        "kinase_network_nodes": "signalome/kinase_network_nodes.csv",
        "module_assignments": "signalome/module_assignments.csv",
        "signalome_modules": "signalome/signalome_modules.csv",
    }
    signalome_metadata = manifest["signalome_outputs"]["metadata"]
    assert signalome_metadata["expanded_signalome_present"] is True
    assert signalome_metadata["kinase_network_nodes_present"] is True
    diagnostics_payload = signalome_metadata["module_selection_diagnostics"]
    assert diagnostics_payload["strategy"] in {
        "correlation_thresholds",
        "explicit_module_count",
    }
    assert int(diagnostics_payload["selected_module_count"]) >= 1
    assert isinstance(diagnostics_payload["candidate_scores"], dict)
    preconditioning_payload = signalome_metadata["score_preconditioning_diagnostics"]
    assert preconditioning_payload["policy"] == "allow_and_report"
    assert int(preconditioning_payload["input_row_count"]) >= 0
    assert int(preconditioning_payload["dropped_all_missing_row_count"]) >= 0
    assert int(preconditioning_payload["retained_row_count"]) >= 0
    network_correlation_payload = signalome_metadata["network_correlation_diagnostics"]
    assert int(network_correlation_payload["total_candidate_correlations"]) >= 0
    assert int(network_correlation_payload["finite_correlations"]) >= 0
    assert int(network_correlation_payload["undefined_correlations"]) >= 0
    assert int(network_correlation_payload["edges_created"]) >= 0
    assert int(network_correlation_payload["edges_skipped_non_finite_correlation"]) >= 0
    assert "provenance" in manifest
    provenance = manifest["provenance"]
    assert provenance["workflow_name"] == "signalome_workflow"
    assert provenance["environment"]["package_name"] == "phospy"
    assert "signalome_config" in provenance["workflow_parameters"]
    signalome_config = provenance["workflow_parameters"]["signalome_config"]
    assert signalome_config["cluster_tree_backend"] == "exact"
    assert signalome_config["candidate_scoring_backend"] == "full"
    assert signalome_config["max_exact_cluster_tree_sites"] == 2000
    assert signalome_config["max_full_correlation_sites"] == 2000
    assert "scale_guard" in provenance["workflow_parameters"]
    scale_guard = provenance["workflow_parameters"]["scale_guard"]
    assert scale_guard["site_count"] >= 1
    assert scale_guard["cluster_tree_backend"] == "exact"
    assert scale_guard["candidate_scoring_backend"] == "full"
    assert scale_guard["max_exact_cluster_tree_sites"] == 2000
    assert scale_guard["max_full_correlation_sites"] == 2000
    assert scale_guard["scale_guard_passed"] is True
    assert scale_guard["exact_cluster_tree_built"] is True
    assert scale_guard["candidate_scoring_mode"] == "full"
    assert scale_guard["candidate_scoring_evaluated"] is True
    assert scale_guard["candidate_scoring_skip_reason"] is None
    assert scale_guard["candidate_scoring_sampling"] is None
    assert (
        scale_guard["candidate_scoring_applies_to"]
        == "candidate_module_count_evaluation_only"
    )
    assert scale_guard["final_module_assignment_backend"] == "exact_cluster_tree_cut"
    assert scale_guard["final_module_assignment_uses_candidate_scoring"] is False
    assert "module_selection_diagnostics" in provenance["workflow_parameters"]
    assert provenance["workflow_parameters"]["upstream_kinase_provenance"] is not None
    output_names = {entry["name"] for entry in provenance["output_tables"]}
    assert "outputs.signalome.module_assignments" in output_names
    assert "outputs.signalome.signalome_modules" in output_names
    assert "outputs.signalome.kinase_network.edges" in output_names
    assert "outputs.signalome.expanded_signalome" in output_names

    loaded = load_signalome_workflow_bundle(bundle_root)
    assert loaded.result.expanded_signalome is not None


def test_signalome_config_snapshot_accepts_compatibility_cutoff_payload() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        {"signalome_config": {"signalome_cutoff": 0.6}}
    )
    assert snapshot.signalome_config.substrate_support_cutoff == pytest.approx(0.6)
    assert snapshot.signalome_config.network_correlation_threshold == pytest.approx(0.6)
    assert snapshot.signalome_config.network_policy == "signed"
    assert (
        snapshot.signalome_config.assignment_policy
        == SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    )
    assert snapshot.signalome_config.score_preconditioning_policy == "allow_and_report"
    assert snapshot.signalome_config.cluster_tree_backend == "exact"
    assert snapshot.signalome_config.candidate_scoring_backend == "full"
    assert snapshot.signalome_config.max_exact_cluster_tree_sites == 2000
    assert snapshot.signalome_config.max_full_correlation_sites == 2000


def test_signalome_bundle_manifest_tracks_absent_expanded_output_when_none(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    without_expanded = SignalomeWorkflowResult._from_owned(
        dataset=result.dataset,
        kinase_result=result.kinase_result,
        module_assignments=result.module_assignments,
        signalome_modules=result.signalome_modules,
        kinase_network=result.kinase_network,
        module_selection_diagnostics=result.module_selection_diagnostics,
        expanded_signalome=None,
    )
    bundle_root = tmp_path / "signalome_bundle_no_expanded"

    save_signalome_workflow_bundle(
        without_expanded,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))

    assert manifest["signalome_outputs"]["tables"]["expanded_signalome"] is None
    assert manifest["signalome_outputs"]["metadata"]["expanded_signalome_present"] is (
        False
    )
    loaded = load_signalome_workflow_bundle(bundle_root)
    assert loaded.result.expanded_signalome is None


def test_signalome_bundle_loads_legacy_manifest_without_provenance(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / "signalome_bundle_legacy"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.pop("provenance", None)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    loaded = load_signalome_workflow_bundle(bundle_root)
    assert loaded.result.provenance is None


def _build_signalome_request_and_result():
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
    request = SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )
    result = SignalomeWorkflow().run(request)
    return request, result


def _assert_signalome_result_equal(left, right) -> None:
    assert left.dataset.organism == right.dataset.organism
    assert left.dataset.intensity_scale_state == right.dataset.intensity_scale_state
    assert left.dataset.processing_state == right.dataset.processing_state
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
        left.kinase_result.scoring_result.rank_weighted_fusion_scores,
        right.kinase_result.scoring_result.rank_weighted_fusion_scores,
    )
    _assert_optional_frame_equal(
        left.kinase_result.scoring_result.score_fusion_weights,
        right.kinase_result.scoring_result.score_fusion_weights,
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
    _assert_optional_frame_equal(
        left.kinase_network.candidate_correlations,
        right.kinase_network.candidate_correlations,
    )
    assert (
        left.kinase_network.correlation_diagnostics
        == right.kinase_network.correlation_diagnostics
    )
    assert left.module_selection_diagnostics == right.module_selection_diagnostics
    assert (
        left.score_preconditioning_diagnostics
        == right.score_preconditioning_diagnostics
    )
    _assert_optional_frame_equal(
        left.expanded_signalome,
        right.expanded_signalome,
    )


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
