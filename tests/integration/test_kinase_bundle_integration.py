from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    KinaseWorkflow,
)
from phospy.api import (
    DatasetBuildRequest,
    DatasetIntensityTransformConfig,
    DatasetPreprocessingConfig,
    DatasetTotalProteinCorrectionConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.io.bundles.kinase import (
    KINASE_BUNDLE_MANIFEST_VERSION,
    KinaseWorkflowConfigSnapshot,
    load_kinase_workflow_bundle,
    save_kinase_workflow_bundle,
)

pytestmark = pytest.mark.integration


def test_kinase_bundle_round_trip_preserves_outputs_and_config(
    tmp_path: Path,
) -> None:
    request = _build_request(activity=True)
    result = KinaseWorkflow().run(request)
    config_snapshot = KinaseWorkflowConfigSnapshot.from_request(request)
    bundle_root = tmp_path / "kinase_bundle"

    written = save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=config_snapshot,
        output_format="csv",
    )

    assert (bundle_root / "manifest.json") == written["manifest"]
    loaded = load_kinase_workflow_bundle(bundle_root)

    assert loaded.manifest_version == KINASE_BUNDLE_MANIFEST_VERSION
    assert loaded.config_snapshot == config_snapshot
    assert loaded.result.provenance == result.provenance
    _assert_kinase_result_equal(loaded.result, result)


def test_kinase_bundle_round_trip_preserves_total_protein_correction_state(
    tmp_path: Path,
) -> None:
    request = _build_request_with_subtract_log_total(activity=False)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle_total_correction"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    loaded = load_kinase_workflow_bundle(bundle_root)

    correction = loaded.result.dataset.processing_state.total_protein_correction
    assert loaded.result.dataset.intensity_scale_state.label == "log2"
    assert (
        loaded.result.dataset.intensity_scale_state.quantity.value
        == "phospho_total_log_ratio"
    )
    assert correction.policy == "subtract_log_total"
    assert correction.applied is True
    assert correction.formula == "log2_phospho - log2_total"
    assert correction.requires_log_scale is True
    assert correction.input_scale == "log2"
    assert correction.output_scale == "log2_ratio"
    assert correction.quantitative_meaning == "phospho_total_log_ratio"
    assert correction.diagnostics is not None
    assert correction.diagnostics.get("matched_rows") == len(
        result.dataset.phospho.index
    )
    assert isinstance(correction.diagnostics.get("input_phospho_hash"), str)
    assert isinstance(correction.diagnostics.get("output_phospho_hash"), str)
    pd.testing.assert_frame_equal(
        loaded.result.dataset.phospho,
        result.dataset.phospho,
        check_dtype=False,
        check_names=False,
    )


def test_kinase_bundle_manifest_v1_is_explicit(tmp_path: Path) -> None:
    request = _build_request(activity=True)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )

    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["manifest_version"] == KINASE_BUNDLE_MANIFEST_VERSION
    assert manifest["bundle_type"] == "kinase_workflow_result"
    assert manifest["table_format"] == "csv"
    assert manifest["config_snapshot"] == "config/snapshot.json"
    assert manifest["dataset"]["tables"] == {
        "phospho": "dataset/phospho.csv",
        "sample_metadata": "dataset/sample_metadata.csv",
        "site_metadata": "dataset/site_metadata.csv",
        "total": "dataset/total.csv",
    }
    assert manifest["resolved_references"]["tables"] == {
        "kinase_substrate_map": "references/kinase_substrate_map.csv",
        "site_sequences": "references/site_sequences.csv",
    }
    assert manifest["outputs"]["scoring"]["tables"] == {
        "rank_weighted_fusion_scores": "scoring/rank_weighted_fusion_scores.csv",
        "motif_scores": None,
        "profile_scores": "scoring/profile_scores.csv",
        "score_fusion_weights": None,
    }
    assert manifest["outputs"]["prediction"]["tables"] == {
        "pred_mat": "prediction/pred_mat.csv",
        "substrate_list": "prediction/substrate_list.csv",
    }
    assert manifest["outputs"]["activity"] == {
        "enabled": True,
        "tables": {
            "weighted_activity": "activity/weighted_activity.csv",
            "thresholded_substrate_mean_activity": "activity/thresholded_substrate_mean_activity.csv",
            "thresholded_substrate_counts": "activity/thresholded_substrate_counts.csv",
            "target_counts": "activity/target_counts.csv",
            "target_table": "activity/target_table.csv",
        },
    }
    assert "provenance" in manifest
    provenance = manifest["provenance"]
    assert provenance["environment"]["package_name"] == "phospy"
    assert provenance["environment"]["python_version"]
    dependency_versions = provenance["environment"]["dependency_versions"]
    assert {"numpy", "pandas", "scikit-learn"}.issubset(set(dependency_versions.keys()))
    input_names = {entry["name"] for entry in provenance["input_tables"]}
    output_names = {entry["name"] for entry in provenance["output_tables"]}
    assert "dataset.phospho" in input_names
    assert "references.kinase_substrate_map" in input_names
    assert "outputs.scoring.profile_scores" in output_names
    assert "outputs.prediction.pred_mat" in output_names
    assert provenance["workflow_name"] == "kinase_workflow"
    assert "prediction_config" in provenance["workflow_parameters"]


def test_kinase_bundle_round_trip_supports_disabled_activity(
    tmp_path: Path,
) -> None:
    request = _build_request(activity=False)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    loaded = load_kinase_workflow_bundle(bundle_root)

    assert result.activity_result is None
    assert loaded.result.activity_result is None
    assert loaded.config_snapshot.activity_config is None
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["outputs"]["activity"] == {
        "enabled": False,
        "tables": {
            "weighted_activity": None,
            "thresholded_substrate_mean_activity": None,
            "thresholded_substrate_counts": None,
            "target_counts": None,
            "target_table": None,
        },
    }


def test_kinase_bundle_loads_legacy_manifest_without_provenance(
    tmp_path: Path,
) -> None:
    request = _build_request(activity=False)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle_legacy"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.pop("provenance", None)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    loaded = load_kinase_workflow_bundle(bundle_root)
    assert loaded.result.provenance is None


def test_kinase_bundle_loads_legacy_minimal_total_correction_state(
    tmp_path: Path,
) -> None:
    request = _build_request_with_subtract_log_total(activity=False)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle_total_correction_legacy"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    correction_payload = manifest["dataset"]["metadata"]["processing_state"][
        "total_protein_correction"
    ]
    manifest["dataset"]["metadata"]["processing_state"]["total_protein_correction"] = {
        "policy": correction_payload["policy"],
        "applied": correction_payload["applied"],
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    loaded = load_kinase_workflow_bundle(bundle_root)
    correction = loaded.result.dataset.processing_state.total_protein_correction
    assert (
        loaded.result.dataset.intensity_scale_state.quantity.value
        == "phospho_total_log_ratio"
    )
    assert correction.policy == "subtract_log_total"
    assert correction.applied is True
    assert correction.formula is None
    assert correction.requires_log_scale is None
    assert correction.input_scale is None
    assert correction.output_scale is None
    assert correction.quantitative_meaning is None
    assert correction.diagnostics is None
    pd.testing.assert_frame_equal(
        loaded.result.dataset.phospho,
        result.dataset.phospho,
        check_dtype=False,
        check_names=False,
    )


def test_kinase_bundle_loads_legacy_manifest_without_quantitative_meaning_field(
    tmp_path: Path,
) -> None:
    request = _build_request_with_subtract_log_total(activity=False)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle_legacy_quantity"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["dataset"]["metadata"]["intensity_scale_state"].pop("quantity", None)
    manifest["dataset"]["metadata"]["processing_state"]["intensity_scale"].pop(
        "quantity", None
    )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    loaded = load_kinase_workflow_bundle(bundle_root)
    assert loaded.result.dataset.intensity_scale_state.label == "log2"
    assert (
        loaded.result.dataset.intensity_scale_state.quantity.value
        == "phospho_total_log_ratio"
    )


def _build_request(*, activity: bool) -> KinaseWorkflowRequest:
    phospho = pd.DataFrame(
        {
            "sample_a": [1.20, 0.85, 1.05],
            "sample_b": [1.05, 0.70, 1.25],
        },
        index=["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
                "MDFGLCKEGIKDGATMKLCKRERANWQPWQ",
            ],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {
            "group": ["treated", "control"],
        },
        index=phospho.columns.copy(),
    )
    total = pd.DataFrame(
        {
            "sample_a": [1.0, 1.1],
            "sample_b": [1.2, 1.0],
        },
        index=["MAPK14", "AKT1"],
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            organism=Organism.RAT,
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "GSK3B;S9;",
                    "AKT1;T308;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": site_metadata.loc[:, "site_sequence"],
            },
            index=pd.Index(site_metadata.index.copy(), name="site_id"),
        ),
    )
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=(
            KinaseActivityConfig(enabled=True, threshold=0.5, min_substrates=2)
            if activity
            else None
        ),
    )


def _build_request_with_subtract_log_total(*, activity: bool) -> KinaseWorkflowRequest:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0, 7.0, 3.0],
            "sample_b": [31.0, 15.0, 7.0],
        },
        index=["MAPK14;Y182;", "AKT1;T308;", "GSK3B;S9;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "MDFGLCKEGIKDGATMKLCKRERANWQPWQ",
                "RARTSSFAEPGGGGGGGGGPGGSASPARPAR",
            ],
        },
        index=phospho.index.copy(),
    )
    sample_metadata = pd.DataFrame(
        {
            "group": ["treated", "control"],
        },
        index=phospho.columns.copy(),
    )
    total = pd.DataFrame(
        {
            "sample_a": [3.0, 1.0, 1.0],
            "sample_b": [7.0, 3.0, 1.0],
        },
        index=["MAPK14", "AKT1", "GSK3B"],
    )
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            organism=Organism.RAT,
            preprocessing_config=DatasetPreprocessingConfig(
                intensity_transform=DatasetIntensityTransformConfig(
                    policy="log2",
                    pseudocount=1.0,
                ),
                total_protein_correction=DatasetTotalProteinCorrectionConfig(
                    policy="subtract_log_total"
                ),
            ),
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "GSK3B;S9;",
                    "AKT1;T308;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": site_metadata.loc[:, "site_sequence"],
            },
            index=pd.Index(site_metadata.index.copy(), name="site_id"),
        ),
    )
    return KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=(
            KinaseActivityConfig(enabled=True, threshold=0.5, min_substrates=2)
            if activity
            else None
        ),
    )


def _assert_kinase_result_equal(left, right) -> None:
    assert left.dataset.organism == right.dataset.organism
    assert left.dataset.intensity_scale_state == right.dataset.intensity_scale_state
    assert left.dataset.processing_state == right.dataset.processing_state
    assert left.references.organism == right.references.organism

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
        left.dataset.sample_metadata,
        right.dataset.sample_metadata,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.dataset.total,
        right.dataset.total,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.references.kinase_substrate_map,
        right.references.kinase_substrate_map,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.references.site_sequences,
        right.references.site_sequences,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.scoring_result.profile_scores,
        right.scoring_result.profile_scores,
        check_dtype=False,
        check_names=False,
    )
    _assert_optional_frame_equal(
        left.scoring_result.motif_scores,
        right.scoring_result.motif_scores,
    )
    _assert_optional_frame_equal(
        left.scoring_result.rank_weighted_fusion_scores,
        right.scoring_result.rank_weighted_fusion_scores,
    )
    _assert_optional_frame_equal(
        left.scoring_result.score_fusion_weights,
        right.scoring_result.score_fusion_weights,
    )
    pd.testing.assert_frame_equal(
        left.prediction_result.pred_mat,
        right.prediction_result.pred_mat,
        check_dtype=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        left.prediction_result.substrate_list,
        right.prediction_result.substrate_list,
        check_dtype=False,
        check_names=False,
    )
    if left.activity_result is None or right.activity_result is None:
        assert left.activity_result is right.activity_result
        return
    pd.testing.assert_frame_equal(
        left.activity_result.weighted_activity,
        right.activity_result.weighted_activity,
        check_dtype=False,
        check_names=False,
        check_index_type=False,
    )
    pd.testing.assert_frame_equal(
        left.activity_result.thresholded_substrate_mean_activity,
        right.activity_result.thresholded_substrate_mean_activity,
        check_dtype=False,
        check_names=False,
        check_index_type=False,
    )
    pd.testing.assert_series_equal(
        left.activity_result.thresholded_substrate_counts,
        right.activity_result.thresholded_substrate_counts,
        check_dtype=False,
        check_names=False,
        check_index_type=False,
    )
    pd.testing.assert_series_equal(
        left.activity_result.target_counts,
        right.activity_result.target_counts,
        check_dtype=False,
        check_names=False,
        check_index_type=False,
    )
    pd.testing.assert_frame_equal(
        left.activity_result.target_table,
        right.activity_result.target_table,
        check_dtype=False,
        check_names=False,
        check_index_type=False,
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
