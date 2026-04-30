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
    SignalomeWorkflowRequest,
)
from phospy.api.results import SignalomeWorkflowResult
from phospy.errors import (
    PhosPyInputError,
    PhosPyValidationError,
    WorkflowValidationError,
)
from phospy.io.bundles.signalome import (
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
    SignalomeWorkflowConfigSnapshot,
    load_signalome_workflow_bundle,
    save_signalome_workflow_bundle,
)
from phospy.provenance.serialization import to_payload as provenance_to_payload
from phospy.signalomes.clustering import (
    SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE,
)
from tests.support.rewrite_fixture_data import build_rat_l6_dataset
from tests.support.signalome_config import build_signalome_config

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
        "protein_site_context": "signalome/protein_site_context.csv",
        "site_membership": "signalome/site_membership.csv",
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
    alignment_payload = signalome_metadata["alignment_diagnostics"]
    assert int(alignment_payload["dataset_sites"]["provided_count"]) >= 0
    assert int(alignment_payload["dataset_sites"]["retained_count"]) >= 0
    assert int(alignment_payload["dataset_sites"]["dropped_count"]) >= 0
    assert isinstance(alignment_payload["dataset_sites"]["dropped_reasons"], dict)
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
    assert signalome_config["clustering"]["tree_engine"] == "exact"
    assert signalome_config["clustering"]["candidate_scoring_policy"] == "full"
    assert signalome_config["performance"]["max_exact_tree_sites"] == 2000
    assert signalome_config["performance"]["max_full_candidate_scoring_sites"] == 2000
    assert "scale_guard" in provenance["workflow_parameters"]
    scale_guard = provenance["workflow_parameters"]["scale_guard"]
    assert scale_guard["site_count"] >= 1
    assert scale_guard["input_protein_count"] >= 1
    assert scale_guard["input_kinase_count"] >= 1
    assert scale_guard["tree_engine"] == "exact"
    assert scale_guard["tree_generation_backend"] in {
        "scipy_hierarchical_tree",
        "exact_python_tree",
    }
    assert scale_guard["tree_generation_mode"] == "full_exact_tree_construction"
    assert scale_guard["tree_generation_is_approximate"] is False
    assert (
        scale_guard["tree_generation_scope"]
        == "module_count_selection_and_final_assignment"
    )
    assert scale_guard["tree_generation_guard_triggered"] is False
    assert scale_guard["candidate_scoring_policy"] == "full"
    assert scale_guard["candidate_scoring_requested_policy"] == "full"
    assert scale_guard["candidate_scoring_strategy"] in {
        "full",
        "sampled",
        "not_evaluated",
    }
    assert scale_guard["candidate_scoring_is_approximate"] is False
    assert scale_guard["candidate_scoring_guard_triggered"] is False
    assert scale_guard["candidate_scoring_sampled_site_total"] is None
    assert scale_guard["candidate_scoring_sampled_pair_count"] is None
    assert scale_guard["max_exact_tree_sites"] == 2000
    assert scale_guard["max_full_candidate_scoring_sites"] == 2000
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
    assert (
        scale_guard["final_module_assignment_backend"]
        == SIGNALOME_FINAL_MODULE_ASSIGNMENT_BACKEND_EXACT_CLUSTER_TREE
    )
    assert scale_guard["final_module_assignment_uses_candidate_scoring"] is False
    score_semantics = provenance["workflow_parameters"]["signalome_score_semantics"]
    assert score_semantics["downstream_score_source"] in {
        "rank_weighted_fusion_scores",
        "profile_scores",
    }
    assert score_semantics["candidate_scoring_mode"] == "full"
    assert score_semantics["candidate_scoring_is_approximate"] is False
    assert score_semantics["candidate_scoring_sampled_site_total"] is None
    assert score_semantics["candidate_scoring_sampled_pair_count"] is None
    assert (
        score_semantics["candidate_scoring_scope"]
        == "candidate_module_count_evaluation_only"
    )
    assert score_semantics["tree_generation_mode"] == "full_exact_tree_construction"
    assert score_semantics["tree_generation_is_approximate"] is False
    assert (
        score_semantics["tree_generation_scope"]
        == "module_count_selection_and_final_assignment"
    )
    assert score_semantics["tree_generation_backend"] in {
        "scipy_hierarchical_tree",
        "exact_python_tree",
    }
    assert score_semantics["input_sizes"]["site_count"] >= 1
    assert score_semantics["input_sizes"]["protein_count"] >= 1
    assert score_semantics["input_sizes"]["kinase_count"] >= 1
    assert score_semantics["scale_guard_status"] == {
        "exact_tree_guard_triggered": False,
        "candidate_scoring_guard_triggered": False,
        "passed": True,
    }
    assert (
        score_semantics["network_policy"]
        == signalome_config["output"]["network_policy"]
    )
    assert score_semantics["clustering_engine"] == scale_guard["clustering_engine"]
    assert "probabilities" in score_semantics["scientific_interpretation_limits"]
    assert "causal" in score_semantics["scientific_interpretation_limits"]
    thresholds = score_semantics["thresholds_and_limits"]
    assert thresholds["network_correlation_threshold"] == pytest.approx(
        signalome_config["output"]["network_correlation_threshold"]
    )
    assert thresholds["max_exact_tree_sites"] == 2000
    assert thresholds["max_full_candidate_scoring_sites"] == 2000
    assert "module_selection_diagnostics" in provenance["workflow_parameters"]
    assert provenance["workflow_parameters"]["upstream_kinase_provenance"] is not None
    output_names = {entry["name"] for entry in provenance["output_tables"]}
    assert "outputs.signalome.module_assignments" in output_names
    assert "outputs.signalome.signalome_modules" in output_names
    assert "outputs.signalome.kinase_network.edges" in output_names
    assert "outputs.signalome.expanded_signalome" in output_names
    assert "outputs.signalome.site_membership" in output_names
    assert "outputs.signalome.protein_site_context" in output_names

    loaded = load_signalome_workflow_bundle(bundle_root)
    assert loaded.result.expanded_signalome is not None
    assert loaded.result.site_membership is not None
    assert loaded.result.protein_site_context is not None


def test_signalome_config_snapshot_rejects_removed_signalome_cutoff_alias() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="contains unsupported field\\(s\\): signalome_cutoff",
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            {"signalome_config": {"signalome_cutoff": 0.6}}
        )


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
        provenance=result.provenance,
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


def test_signalome_bundle_rejects_manifest_without_provenance(
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

    with pytest.raises(
        PhosPyInputError,
        match=(
            "bundle manifest is missing required field\\(s\\): provenance.*"
            "Regenerate this bundle with the current PhosPy version"
        ),
    ):
        load_signalome_workflow_bundle(bundle_root)


def test_signalome_bundle_rejects_manifest_with_null_provenance(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / "signalome_bundle_null_provenance"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["provenance"] = None
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "bundle manifest.provenance is required.*"
            "Regenerate this bundle with the current PhosPy version"
        ),
    ):
        load_signalome_workflow_bundle(bundle_root)


def test_signalome_bundle_rejects_manifest_without_candidate_correlations_marker(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / "signalome_bundle_missing_candidate_correlations"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["signalome_outputs"]["tables"].pop(
        "kinase_network_candidate_correlations", None
    )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "bundle manifest.signalome_outputs.tables is missing required field\\(s\\): "
            "kinase_network_candidate_correlations.*"
            "Regenerate this bundle with the current PhosPy version"
        ),
    ):
        load_signalome_workflow_bundle(bundle_root)


@pytest.mark.parametrize(
    ("site_membership_present", "protein_site_context_present"),
    (
        (True, False),
        (False, True),
        (False, False),
    ),
)
def test_signalome_bundle_round_trip_handles_optional_sidecar_presence(
    tmp_path: Path,
    site_membership_present: bool,
    protein_site_context_present: bool,
) -> None:
    request, result = _build_signalome_request_and_result()
    assert result.site_membership is not None
    assert result.protein_site_context is not None
    variant = _with_sidecars(
        result,
        site_membership=(
            result.site_membership.copy(deep=True) if site_membership_present else None
        ),
        protein_site_context=(
            result.protein_site_context.copy(deep=True)
            if protein_site_context_present
            else None
        ),
    )
    bundle_root = tmp_path / (
        "signalome_bundle_sidecars_"
        f"site_{int(site_membership_present)}_protein_{int(protein_site_context_present)}"
    )

    save_signalome_workflow_bundle(
        variant,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    loaded = load_signalome_workflow_bundle(bundle_root)

    _assert_optional_frame_equal(
        loaded.result.site_membership,
        variant.site_membership,
    )
    _assert_optional_frame_equal(
        loaded.result.protein_site_context,
        variant.protein_site_context,
    )


def test_signalome_bundle_rejects_manifest_without_site_membership_marker(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / "signalome_bundle_missing_site_membership"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["signalome_outputs"]["tables"].pop("site_membership", None)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "bundle manifest.signalome_outputs.tables is missing required field\\(s\\): "
            "site_membership.*Regenerate this bundle with the current PhosPy version"
        ),
    ):
        load_signalome_workflow_bundle(bundle_root)


def test_signalome_bundle_rejects_malformed_site_membership_manifest_entry(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / "signalome_bundle_malformed_site_membership_entry"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["signalome_outputs"]["tables"]["site_membership"] = {
        "path": "signalome/site_membership.csv"
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    with pytest.raises(
        PhosPyInputError,
        match="bundle manifest.signalome_outputs.tables.site_membership must be a string",
    ):
        load_signalome_workflow_bundle(bundle_root)


def test_signalome_bundle_rejects_missing_declared_site_membership_payload(
    tmp_path: Path,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / "signalome_bundle_missing_site_membership_payload"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    site_membership_path = Path(
        manifest["signalome_outputs"]["tables"]["site_membership"]
    )
    (bundle_root / site_membership_path).unlink()

    with pytest.raises(
        PhosPyInputError,
        match="input file does not exist: .*site_membership\\.csv",
    ):
        load_signalome_workflow_bundle(bundle_root)


@pytest.mark.parametrize(
    ("table_key", "invalid_table", "pattern"),
    (
        (
            "site_membership",
            pd.DataFrame({"invalid": [1]}),
            "signalome_result.site_membership.*missing required columns",
        ),
        (
            "protein_site_context",
            pd.DataFrame({"invalid": [1]}),
            "signalome_result.protein_site_context.*missing required columns",
        ),
    ),
)
def test_signalome_bundle_rejects_invalid_sidecar_table_schema(
    tmp_path: Path,
    table_key: str,
    invalid_table: pd.DataFrame,
    pattern: str,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / f"signalome_bundle_invalid_{table_key}_schema"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    table_path = Path(manifest["signalome_outputs"]["tables"][table_key])
    invalid_table.to_csv(bundle_root / table_path)

    with pytest.raises(WorkflowValidationError, match=pattern):
        load_signalome_workflow_bundle(bundle_root)


@pytest.mark.parametrize(
    ("table_key", "invalid_table", "pattern"),
    (
        (
            "module_assignments",
            pd.DataFrame({"invalid": [1]}),
            "signalome_result.module_assignments.table.*missing required columns",
        ),
        (
            "signalome_modules",
            pd.DataFrame(
                {"K1": [60.0], "K2": [20.0]}, index=pd.Index([1], name="module_id")
            ),
            "signalome_result.signalome_modules.table row totals must be approximately 0.0 or 100.0",
        ),
        (
            "kinase_network_edges",
            pd.DataFrame(
                {"source": ["K1"], "target_kinase": ["K2"], "correlation": [0.8]}
            ),
            "signalome_result.kinase_network.edges.*missing required columns",
        ),
        (
            "kinase_network_candidate_correlations",
            pd.DataFrame(
                {
                    "source_kinase": ["K1"],
                    "target_kinase": ["K2"],
                    "correlation": [0.8],
                    "correlation_status": ["invalid_status"],
                    "valid_observations": [4],
                    "correlation_reason": [None],
                }
            ),
            "signalome_result.kinase_network.candidate_correlations.correlation_status contains unsupported values",
        ),
    ),
)
def test_signalome_bundle_rejects_invalid_public_output_table_schema(
    tmp_path: Path,
    table_key: str,
    invalid_table: pd.DataFrame,
    pattern: str,
) -> None:
    request, result = _build_signalome_request_and_result()
    bundle_root = tmp_path / f"signalome_bundle_invalid_{table_key}_schema"

    save_signalome_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=SignalomeWorkflowConfigSnapshot.from_request(request),
    )
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    table_path = Path(manifest["signalome_outputs"]["tables"][table_key])
    invalid_table.to_csv(bundle_root / table_path)

    with pytest.raises(PhosPyValidationError, match=pattern):
        load_signalome_workflow_bundle(bundle_root)


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
        config=build_signalome_config(substrate_support_cutoff=0.5),
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
    assert left.alignment_diagnostics == right.alignment_diagnostics
    _assert_optional_frame_equal(
        left.expanded_signalome,
        right.expanded_signalome,
    )
    _assert_optional_frame_equal(
        left.site_membership,
        right.site_membership,
    )
    _assert_optional_frame_equal(
        left.protein_site_context,
        right.protein_site_context,
    )


def _with_sidecars(
    result: SignalomeWorkflowResult,
    *,
    site_membership: pd.DataFrame | None,
    protein_site_context: pd.DataFrame | None,
) -> SignalomeWorkflowResult:
    return SignalomeWorkflowResult._from_owned(
        dataset=result.dataset,
        kinase_result=result.kinase_result,
        module_assignments=result.module_assignments,
        signalome_modules=result.signalome_modules,
        kinase_network=result.kinase_network,
        module_selection_diagnostics=result.module_selection_diagnostics,
        score_preconditioning_diagnostics=result.score_preconditioning_diagnostics,
        alignment_diagnostics=result.alignment_diagnostics,
        expanded_signalome=result.expanded_signalome,
        site_membership=site_membership,
        protein_site_context=protein_site_context,
        provenance=result.provenance,
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
