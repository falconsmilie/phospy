from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
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
    DatasetTotalProteinCorrectionIdentityConfig,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.configs import (
    DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED,
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
)
from phospy.contracts.results import KinaseWorkflowResult
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles.kinase import (
    KINASE_BUNDLE_MANIFEST_VERSION,
    KinaseWorkflowConfigSnapshot,
    load_kinase_workflow_bundle,
    save_kinase_workflow_bundle,
)
from phospy.provenance.models import RunProvenance
from phospy.science.activities.method_contracts import (
    kinase_activity_method_quantitative_input_contract,
)
from phospy.science.activities.methods import (
    KseaZScoreActivityMethod,
    SimplifiedWeightedSubstrateActivityMethod,
    SsgseaSubstrateEnrichmentActivityMethod,
)
from phospy.science.activities.models import (
    KinaseActivityInputs,
    KinaseActivityResult,
    PredMatOverlapSummary,
)
from phospy.science.activities.semantics import (
    ActivityAggregationMetadata,
    ActivityAggregationRecord,
    ActivityInputMatrix,
    ActivityProfileAxis,
    ActivityQuantitativeSemantics,
)
from phospy.science.quantitative_method_contracts import (
    ResolvedMethodQuantitativeInputContract,
)
from tests.support.site_keys import site_key_index_from_display_ids

pytestmark = pytest.mark.integration


def _collect_keys(value: object) -> set[str]:
    if isinstance(value, dict):
        keys = {str(key) for key in value}
        for item in value.values():
            keys.update(_collect_keys(item))
        return keys
    if isinstance(value, list):
        keys: set[str] = set()
        for item in value:
            keys.update(_collect_keys(item))
        return keys
    return set()


def _sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _entry_path(entry: object) -> str | None:
    if entry is None:
        return None
    assert isinstance(entry, dict)
    path = entry["path"]
    assert isinstance(path, str)
    return path


def _table_paths(tables: dict[str, object]) -> dict[str, str | None]:
    return {key: _entry_path(value) for key, value in tables.items()}


def _assert_file_entry(
    entry: object,
    *,
    bundle_root: Path,
    relative_path: str,
    logical_type: str,
    shape: tuple[int, int] | None = None,
) -> None:
    assert isinstance(entry, dict)
    path = bundle_root / relative_path
    assert entry["path"] == relative_path
    assert entry["logical_type"] == logical_type
    assert entry["byte_size"] == path.stat().st_size
    assert entry["sha256"] == _sha256_path(path)
    if shape is not None:
        assert entry["shape"] == {"rows": shape[0], "columns": shape[1]}


def _iter_file_entries(value: object):
    if isinstance(value, dict):
        if {"path", "sha256", "byte_size"}.issubset(value):
            yield value
            return
        for item in value.values():
            yield from _iter_file_entries(item)
        return
    if isinstance(value, list):
        for item in value:
            yield from _iter_file_entries(item)


def _assert_manifest_covers_bundle_payload_files(
    manifest: dict[str, object],
    *,
    bundle_root: Path,
) -> None:
    entries = list(_iter_file_entries(manifest))
    declared_paths = {str(entry["path"]) for entry in entries}
    observed_paths = {
        path.relative_to(bundle_root).as_posix()
        for path in bundle_root.rglob("*")
        if path.is_file() and path.name != "manifest.json"
    }
    assert declared_paths == observed_paths
    for entry in entries:
        path = bundle_root / str(entry["path"])
        assert entry["byte_size"] == path.stat().st_size
        assert entry["sha256"] == _sha256_path(path)


def _refresh_manifest_table_entry(
    entry: dict[str, object],
    *,
    bundle_root: Path,
    table: pd.DataFrame,
) -> None:
    path = bundle_root / str(entry["path"])
    entry["byte_size"] = path.stat().st_size
    entry["sha256"] = _sha256_path(path)
    entry["shape"] = {"rows": int(table.shape[0]), "columns": int(table.shape[1])}


def _save_basic_kinase_bundle(
    tmp_path: Path,
    *,
    bundle_name: str,
    activity: bool = False,
) -> tuple[Path, dict[str, object]]:
    request = _build_request(activity=activity)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / bundle_name
    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    return bundle_root, manifest


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
    assert loaded.result.dataset.provenance is not None
    assert (
        loaded.result.dataset.provenance.workflow_name
        == "analysis_ready_dataset_direct_construction"
    )
    construction = loaded.result.dataset.provenance.workflow_parameters.get(
        "construction"
    )
    assert isinstance(construction, dict)
    assert construction["method"] == "AnalysisReadyPhosphoDataset.from_trusted_tables"
    assert loaded.result.dataset.trusted_construction_assertions is not None
    assert loaded.result.dataset.trusted_construction_assertions.all_required_assertions_present
    _assert_kinase_result_equal(loaded.result, result)
    assert loaded.result.activity_result is not None
    with pytest.warns(
        DeprecationWarning,
        match="KinaseActivityResult.activity_scores.*activity_matrix",
    ):
        activity_scores = loaded.result.activity_result.activity_scores
    pd.testing.assert_frame_equal(
        activity_scores,
        loaded.result.activity_result.activity_matrix,
    )
    with pytest.warns(
        DeprecationWarning,
        match="KinaseActivityResult.weighted_activity.*activity_matrix",
    ):
        weighted_activity = loaded.result.activity_result.weighted_activity
    pd.testing.assert_frame_equal(
        weighted_activity,
        loaded.result.activity_result.activity_matrix,
    )


def test_kinase_bundle_round_trip_preserves_adaptive_prediction_seed(
    tmp_path: Path,
) -> None:
    base_request = _build_request(activity=False)
    seeded_request = KinaseWorkflowRequest(
        dataset=base_request.dataset,
        references=base_request.references,
        scoring_config=base_request.scoring_config,
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            mode="adaptive_ensemble",
            n_iterations=2,
            random_state=19,
        ),
        activity_config=base_request.activity_config,
    )
    result = KinaseWorkflow().run(seeded_request)
    bundle_root = tmp_path / "kinase_bundle_adaptive_seed"
    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(seeded_request),
    )

    loaded = load_kinase_workflow_bundle(bundle_root)
    assert loaded.config_snapshot.prediction_config.mode == "adaptive_ensemble"
    assert loaded.config_snapshot.prediction_config.random_state == 19
    assert loaded.result.provenance is not None
    assert (
        loaded.result.provenance.workflow_parameters["prediction_config"][
            "random_state"
        ]
        == 19
    )
    assert loaded.result.provenance.random_state == 19


def test_kinase_bundle_round_trip_preserves_stage_table_fingerprints(
    tmp_path: Path,
) -> None:
    request = _build_request(activity=False)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle_stage_provenance"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    loaded = load_kinase_workflow_bundle(bundle_root)

    assert result.provenance is not None
    assert loaded.result.provenance is not None
    original_stage = next(
        stage
        for stage in result.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )
    restored_stage = next(
        stage
        for stage in loaded.result.provenance.preprocessing_stages
        if stage.stage == "missing_data"
    )
    assert original_stage.schema_version >= 2
    assert original_stage.consumed_input_tables
    assert original_stage.produced_output_tables
    assert restored_stage.schema_version == original_stage.schema_version
    assert restored_stage.consumed_input_tables == original_stage.consumed_input_tables
    assert (
        restored_stage.produced_output_tables == original_stage.produced_output_tables
    )
    assert restored_stage.backend == original_stage.backend
    assert restored_stage.random_seed == original_stage.random_seed
    assert restored_stage.determinism == original_stage.determinism


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
    assert correction.diagnostics.get("diagnostics_schema_version") == 1
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


def test_kinase_bundle_round_trip_preserves_mixed_total_protein_quantitative_meaning(
    tmp_path: Path,
) -> None:
    request = _build_request_with_subtract_log_total_and_uncorrected_rows(
        activity=False
    )
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_bundle_total_correction_mixed"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    loaded = load_kinase_workflow_bundle(bundle_root)

    mixed_meaning = "mixed_phospho_total_log_ratio_and_phosphosite_log_abundance"
    assert loaded.result.dataset.intensity_scale_state.quantity is not None
    assert loaded.result.dataset.intensity_scale_state.quantity.value == mixed_meaning
    correction = loaded.result.dataset.processing_state.total_protein_correction
    assert correction.quantitative_meaning == mixed_meaning
    assert correction.diagnostics is not None
    assert correction.diagnostics.get("uncorrected_row_count") == 1
    assert correction.diagnostics.get("unmatched_policy") == "allow_uncorrected"

    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    correction_payload = manifest["dataset"]["metadata"]["processing_state"][
        "total_protein_correction"
    ]
    assert correction_payload["quantitative_meaning"] == mixed_meaning
    assert correction_payload["diagnostics"]["quantitative_meaning"] == mixed_meaning


def test_kinase_bundle_manifest_v3_is_explicit_and_content_addressed(
    tmp_path: Path,
) -> None:
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
    _assert_file_entry(
        manifest["config_snapshot"],
        bundle_root=bundle_root,
        relative_path="config/snapshot.json",
        logical_type="config_snapshot",
    )
    assert _table_paths(manifest["dataset"]["tables"]) == {
        "phospho": "dataset/phospho.csv",
        "sample_metadata": "dataset/sample_metadata.csv",
        "site_metadata": "dataset/site_metadata.csv",
        "total": "dataset/total.csv",
    }
    _assert_file_entry(
        manifest["dataset"]["tables"]["phospho"],
        bundle_root=bundle_root,
        relative_path="dataset/phospho.csv",
        logical_type="dataset.phospho",
        shape=result.dataset.phospho.shape,
    )
    _assert_manifest_covers_bundle_payload_files(manifest, bundle_root=bundle_root)
    assert _table_paths(manifest["resolved_references"]["tables"]) == {
        "kinase_substrate_map": "references/kinase_substrate_map.csv",
        "site_sequences": "references/site_sequences.csv",
    }
    assert _table_paths(manifest["outputs"]["scoring"]["tables"]) == {
        "rank_weighted_fusion_scores": "scoring/rank_weighted_fusion_scores.csv",
        "motif_scores": None,
        "profile_scores": "scoring/profile_scores.csv",
        "score_fusion_weights": None,
    }
    assert _table_paths(manifest["outputs"]["prediction"]["tables"]) == {
        "pred_mat": "prediction/pred_mat.csv",
        "substrate_list": "prediction/substrate_list.csv",
    }
    assert manifest["outputs"]["activity"]["enabled"] is True
    assert manifest["outputs"]["activity"]["method"] == {
        "activity_method_id": "simplified_weighted_substrate_activity_v1",
        "activity_method_family": "heuristic_weighted_substrate_score",
        "activity_method_label": "simplified weighted substrate activity-like score",
        "is_ksea": False,
        "is_phosr_kinase_activity_equivalent": False,
    }
    assert manifest["outputs"]["activity"]["summary"] == (
        None
        if result.activity_result is None
        or result.activity_result.method_summary is None
        else result.activity_result.method_summary.to_payload()
    )
    assert result.activity_result is not None
    assert manifest["outputs"]["activity"]["input_semantics"] == (
        result.activity_result.input_semantics.to_payload()
    )
    assert manifest["outputs"]["activity"]["profile_metadata"] == (
        result.activity_result.profile_metadata.to_payload()
    )
    assert _table_paths(manifest["outputs"]["activity"]["tables"]) == {
        "weighted_activity": "activity/weighted_activity.csv",
        "thresholded_substrate_mean_activity": "activity/thresholded_substrate_mean_activity.csv",
        "thresholded_substrate_counts": "activity/thresholded_substrate_counts.csv",
        "activity_substrate_counts": None,
        "target_counts": "activity/target_counts.csv",
        "target_table": "activity/target_table.csv",
        "statistics_table": None,
    }
    correction_diagnostics = manifest["dataset"]["metadata"]["processing_state"][
        "total_protein_correction"
    ]["diagnostics"]
    assert correction_diagnostics["diagnostics_schema_version"] == 1
    assert correction_diagnostics["quantitative_meaning"] == "phosphosite_abundance"
    assert "provenance" in manifest
    provenance = manifest["provenance"]
    provenance_keys = _collect_keys(provenance)
    assert "hash_algorithm" not in provenance_keys
    assert "hash_value" not in provenance_keys
    assert "is_deterministic" not in provenance_keys
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
        "method": None,
        "summary": None,
        "input_semantics": None,
        "profile_metadata": None,
        "tables": {
            "weighted_activity": None,
            "thresholded_substrate_mean_activity": None,
            "thresholded_substrate_counts": None,
            "activity_substrate_counts": None,
            "target_counts": None,
            "target_table": None,
            "statistics_table": None,
        },
    }


@pytest.mark.parametrize(
    "case",
    [
        pytest.param(
            {
                "name": "weighted_sample_abundance",
                "activity_result": "weighted_sample",
                "method": KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
                "scale": "linear",
                "meaning": "phosphosite_abundance",
            },
            id="weighted-sample-abundance",
        ),
        pytest.param(
            {
                "name": "weighted_condition_summary_abundance",
                "activity_result": "weighted_condition_summary",
                "method": KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
                "scale": "linear",
                "meaning": "phosphosite_abundance",
            },
            id="weighted-condition-summary-abundance",
        ),
        pytest.param(
            {
                "name": "ksea_sample_log_abundance",
                "activity_result": "ksea_sample_log_abundance",
                "method": KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
                "scale": "log2",
                "meaning": "phosphosite_log_abundance",
            },
            id="ksea-sample-log-abundance",
        ),
        pytest.param(
            {
                "name": "ksea_contrast_log_fold_change",
                "activity_result": "ksea_contrast_log_fold_change",
                "method": KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
                "scale": "log2",
                "meaning": "contrast_log2_fold_change",
            },
            id="ksea-contrast-log-fold-change",
        ),
        pytest.param(
            {
                "name": "ksea_standardised_effect",
                "activity_result": "ksea_standardised_effect",
                "method": KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
                "scale": "log2",
                "meaning": "differential_effect_size",
            },
            id="ksea-standardised-effect",
        ),
        pytest.param(
            {
                "name": "ssgsea_contrast_log_fold_change",
                "activity_result": "ssgsea_contrast_log_fold_change",
                "method": KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
                "scale": "log2",
                "meaning": "contrast_log2_fold_change",
            },
            id="ssgsea-contrast-log-fold-change",
        ),
        pytest.param(
            {
                "name": "ssgsea_standardised_effect",
                "activity_result": "ssgsea_standardised_effect",
                "method": KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT,
                "scale": "log2",
                "meaning": "differential_effect_size",
            },
            id="ssgsea-standardised-effect",
        ),
    ],
)
def test_kinase_bundle_round_trip_preserves_exact_activity_semantics(
    tmp_path: Path,
    case: dict[str, object],
) -> None:
    base_request = _build_request(activity=False)
    base_result = KinaseWorkflow().run(base_request)
    activity_result_factory = _activity_result_factory(str(case["activity_result"]))
    activity_result = activity_result_factory()
    assert isinstance(activity_result, KinaseActivityResult)
    method = str(case["method"])
    result = _replace_activity_result_with_semantic_provenance(
        base_result,
        activity_result=activity_result,
        method=method,
        resolved_scale=str(case["scale"]),
        resolved_meaning=str(case["meaning"]),
    )
    request = replace(
        base_request,
        activity_config=_activity_config_for_method(method),
    )
    bundle_root = tmp_path / f"kinase_bundle_{case['name']}"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    loaded = load_kinase_workflow_bundle(bundle_root)

    assert loaded.result.activity_result is not None
    _assert_activity_result_semantics_round_tripped(
        loaded.result.activity_result,
        activity_result,
    )
    assert loaded.result.provenance is not None
    _assert_activity_provenance_matches_result(loaded.result)


@pytest.mark.parametrize(
    ("scenario", "activity_factory_name", "mutation", "pattern"),
    [
        pytest.param(
            "missing_input_semantics",
            "weighted_sample",
            "missing_input_semantics",
            "bundle manifest.outputs.activity is missing required field\\(s\\): "
            "input_semantics",
            id="missing-input-semantics",
        ),
        pytest.param(
            "missing_profile_metadata",
            "weighted_sample",
            "missing_profile_metadata",
            "bundle manifest.outputs.activity is missing required field\\(s\\): "
            "profile_metadata",
            id="missing-profile-metadata",
        ),
        pytest.param(
            "axis_quantitative_semantics_mismatch",
            "weighted_sample",
            "axis_quantitative_semantics_mismatch",
            "bundle manifest.outputs.activity.input_semantics is invalid: "
            "activity input semantics are inconsistent",
            id="axis-quantitative-mismatch",
        ),
        pytest.param(
            "profile_ids_reordered",
            "weighted_sample",
            "profile_ids_reordered",
            "bundle manifest.outputs.activity.profile_metadata.profile_ids "
            "must exactly match activity/weighted_activity table columns in order",
            id="profile-ids-reordered",
        ),
        pytest.param(
            "unknown_profile_ids",
            "weighted_sample",
            "unknown_profile_ids",
            "bundle manifest.outputs.activity.profile_metadata.profile_ids "
            "must exactly match activity/weighted_activity table columns in order",
            id="unknown-profile-ids",
        ),
        pytest.param(
            "contradictory_condition_ids",
            "weighted_sample",
            "contradictory_condition_ids",
            "bundle manifest.outputs.activity.profile_metadata.condition_ids "
            "must be empty when profile_metadata.axis is 'sample'",
            id="contradictory-condition-ids",
        ),
        pytest.param(
            "missing_condition_summary_aggregation",
            "weighted_condition_summary",
            "missing_condition_summary_aggregation",
            "bundle manifest.outputs.activity.profile_metadata is invalid: "
            "condition-summary activity input requires explicit "
            "ActivityAggregationMetadata",
            id="missing-condition-summary-aggregation",
        ),
        pytest.param(
            "duplicate_aggregation_records",
            "weighted_condition_summary",
            "duplicate_aggregation_records",
            "bundle manifest.outputs.activity.profile_metadata is invalid: "
            "activity_aggregation_metadata.records profile_id values must be unique",
            id="duplicate-aggregation-records",
        ),
        pytest.param(
            "provenance_manifest_axis_conflict",
            "weighted_sample",
            "provenance_manifest_axis_conflict",
            "resolved_activity_profile_axis must agree with "
            "bundle manifest.outputs.activity.input_semantics.profile_axis",
            id="provenance-axis-conflict",
        ),
        pytest.param(
            "provenance_manifest_quantitative_conflict",
            "weighted_sample",
            "provenance_manifest_quantitative_conflict",
            "resolved_activity_quantitative_semantics must agree with "
            "bundle manifest.outputs.activity.input_semantics.quantitative_semantics",
            id="provenance-quantitative-conflict",
        ),
        pytest.param(
            "semantic_metadata_supplied_while_disabled",
            None,
            "semantic_metadata_supplied_while_disabled",
            "bundle manifest.outputs.activity.input_semantics must be null when "
            "activity is disabled",
            id="disabled-semantics",
        ),
        pytest.param(
            "table_tamper_with_digest_refresh",
            "weighted_sample",
            "table_tamper_with_digest_refresh",
            "bundle manifest.outputs.activity.profile_metadata.profile_ids "
            "must exactly match activity/weighted_activity table columns in order",
            id="table-tamper-with-digest-refresh",
        ),
    ],
)
def test_kinase_bundle_loader_rejects_contradictory_activity_semantics(
    tmp_path: Path,
    scenario: str,
    activity_factory_name: str | None,
    mutation: str,
    pattern: str,
) -> None:
    if activity_factory_name is None:
        bundle_root, manifest = _save_basic_kinase_bundle(
            tmp_path,
            bundle_name=f"kinase_bundle_semantic_{scenario}",
            activity=False,
        )
    else:
        activity_factory = _activity_result_factory(activity_factory_name)
        bundle_root, manifest = _save_semantic_activity_bundle(
            tmp_path,
            bundle_name=f"kinase_bundle_semantic_{scenario}",
            activity_result=activity_factory(),
            method=KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
            resolved_scale="linear",
            resolved_meaning="phosphosite_abundance",
        )

    _mutate_activity_manifest(
        manifest,
        bundle_root=bundle_root,
        mutation=mutation,
    )
    (bundle_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    with pytest.raises(PhosPyInputError, match=pattern):
        load_kinase_workflow_bundle(bundle_root)


def test_kinase_bundle_v2_manifest_is_rejected_as_legacy_semantics_schema(
    tmp_path: Path,
) -> None:
    bundle_root, manifest = _save_basic_kinase_bundle(
        tmp_path,
        bundle_name="kinase_bundle_v2_rejected",
        activity=True,
    )
    manifest["manifest_version"] = 2
    activity_payload = manifest["outputs"]["activity"]
    activity_payload.pop("input_semantics", None)
    activity_payload.pop("profile_metadata", None)
    (bundle_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "bundle manifest.manifest_version=2 is a legacy kinase bundle "
            "schema.*activity input semantics.*profile identity.*condition-summary "
            "aggregation metadata.*bundle must be regenerated"
        ),
    ):
        load_kinase_workflow_bundle(bundle_root)


def test_kinase_bundle_loader_rejects_table_tampering(tmp_path: Path) -> None:
    bundle_root, manifest = _save_basic_kinase_bundle(
        tmp_path,
        bundle_name="kinase_bundle_table_tamper",
    )
    phospho_entry = manifest["dataset"]["tables"]["phospho"]
    assert isinstance(phospho_entry, dict)
    phospho_path = bundle_root / str(phospho_entry["path"])
    phospho_bytes = bytearray(phospho_path.read_bytes())
    phospho_bytes[-1] = phospho_bytes[-1] ^ 1
    phospho_path.write_bytes(bytes(phospho_bytes))

    with pytest.raises(
        PhosPyInputError,
        match="digest mismatch: path=dataset/phospho\\.csv",
    ):
        load_kinase_workflow_bundle(bundle_root)


def test_kinase_bundle_loader_rejects_json_tampering(tmp_path: Path) -> None:
    bundle_root, manifest = _save_basic_kinase_bundle(
        tmp_path,
        bundle_name="kinase_bundle_json_tamper",
    )
    config_entry = manifest["config_snapshot"]
    assert isinstance(config_entry, dict)
    config_path = bundle_root / str(config_entry["path"])
    config_bytes = bytearray(config_path.read_bytes())
    config_bytes[-1] = config_bytes[-1] ^ 1
    config_path.write_bytes(bytes(config_bytes))

    with pytest.raises(
        PhosPyInputError,
        match="digest mismatch: path=config/snapshot\\.json",
    ):
        load_kinase_workflow_bundle(bundle_root)


def test_kinase_bundle_loader_rejects_missing_declared_file(tmp_path: Path) -> None:
    bundle_root, manifest = _save_basic_kinase_bundle(
        tmp_path,
        bundle_name="kinase_bundle_missing_file",
    )
    prediction_entry = manifest["outputs"]["prediction"]["tables"]["pred_mat"]
    assert isinstance(prediction_entry, dict)
    (bundle_root / str(prediction_entry["path"])).unlink()

    with pytest.raises(
        PhosPyInputError,
        match="declared file is missing: path=prediction/pred_mat\\.csv",
    ):
        load_kinase_workflow_bundle(bundle_root)


def test_kinase_bundle_loader_rejects_extra_stale_file(tmp_path: Path) -> None:
    bundle_root, _ = _save_basic_kinase_bundle(
        tmp_path,
        bundle_name="kinase_bundle_extra_file",
    )
    (bundle_root / "stale.csv").write_text("stale\n", encoding="utf-8")

    with pytest.raises(
        PhosPyInputError,
        match="undeclared file\\(s\\): stale\\.csv",
    ):
        load_kinase_workflow_bundle(bundle_root)


def test_kinase_bundle_interrupted_write_does_not_publish_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.io.bundles._kinase.writer as kinase_writer

    request = _build_request(activity=False)
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / "kinase_interrupted"
    original_write_json = kinase_writer.write_json

    def fail_manifest_write(path: Path, payload: object, *, label: str) -> None:
        if label == "bundle manifest":
            raise PhosPyInputError("simulated interrupted manifest write")
        original_write_json(path, payload, label=label)

    monkeypatch.setattr(kinase_writer, "write_json", fail_manifest_write)

    with pytest.raises(PhosPyInputError, match="simulated interrupted manifest write"):
        save_kinase_workflow_bundle(
            result,
            bundle_root,
            config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
        )

    assert not bundle_root.exists()
    assert not list(tmp_path.glob(".kinase_interrupted.tmp-*"))


def test_kinase_bundle_overwrite_policy_is_explicit_and_transactional(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.io.bundles._kinase.writer as kinase_writer

    request = _build_request(activity=False)
    result = KinaseWorkflow().run(request)
    config_snapshot = KinaseWorkflowConfigSnapshot.from_request(request)
    bundle_root = tmp_path / "kinase_overwrite"
    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=config_snapshot,
    )

    with pytest.raises(PhosPyInputError, match="Pass overwrite=True"):
        save_kinase_workflow_bundle(
            result,
            bundle_root,
            config_snapshot=config_snapshot,
        )

    original_write_json = kinase_writer.write_json

    def fail_manifest_write(path: Path, payload: object, *, label: str) -> None:
        if label == "bundle manifest":
            raise PhosPyInputError("simulated overwrite interruption")
        original_write_json(path, payload, label=label)

    monkeypatch.setattr(kinase_writer, "write_json", fail_manifest_write)
    with pytest.raises(PhosPyInputError, match="simulated overwrite interruption"):
        save_kinase_workflow_bundle(
            result,
            bundle_root,
            config_snapshot=config_snapshot,
            overwrite=True,
        )
    loaded_after_failed_overwrite = load_kinase_workflow_bundle(bundle_root)
    _assert_kinase_result_equal(loaded_after_failed_overwrite.result, result)

    monkeypatch.setattr(kinase_writer, "write_json", original_write_json)
    stale_file = bundle_root / "stale.csv"
    stale_file.write_text("stale\n", encoding="utf-8")
    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=config_snapshot,
        overwrite=True,
    )

    assert not stale_file.exists()
    loaded = load_kinase_workflow_bundle(bundle_root)
    _assert_kinase_result_equal(loaded.result, result)


def test_kinase_bundle_v1_manifest_is_rejected_with_migration_message(
    tmp_path: Path,
) -> None:
    bundle_root, manifest = _save_basic_kinase_bundle(
        tmp_path,
        bundle_name="kinase_bundle_v1_rejected",
    )
    manifest["manifest_version"] = 1
    (bundle_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    with pytest.raises(
        PhosPyInputError,
        match=(
            "Legacy kinase bundle schemas are no longer supported.*"
            "Regenerate the bundle with the current PhosPy version.*"
            "unsupported bundle manifest version '1'; expected "
            f"{KINASE_BUNDLE_MANIFEST_VERSION}"
        ),
    ):
        load_kinase_workflow_bundle(bundle_root)


@pytest.mark.parametrize(
    (
        "scenario_name",
        "requires_total_correction_payload",
        "mutation_kind",
        "pattern",
    ),
    [
        pytest.param(
            "missing_provenance",
            False,
            "missing_provenance",
            (
                "Legacy kinase bundle schemas are no longer supported.*"
                "Regenerate the bundle with the current PhosPy version.*"
                "bundle manifest is missing required field\\(s\\): provenance.*"
            ),
            id="missing-provenance",
        ),
        pytest.param(
            "null_provenance",
            False,
            "null_provenance",
            (
                "Legacy kinase bundle schemas are no longer supported.*"
                "Regenerate the bundle with the current PhosPy version.*"
                "bundle manifest.provenance is required.*"
            ),
            id="null-provenance",
        ),
        pytest.param(
            "missing_activity_enabled",
            False,
            "missing_activity_enabled",
            (
                "Legacy kinase bundle schemas are no longer supported.*"
                "Regenerate the bundle with the current PhosPy version.*"
                "bundle manifest.outputs.activity is missing required field\\(s\\): "
                "enabled"
            ),
            id="missing-activity-enabled",
        ),
        pytest.param(
            "legacy_provenance_hash_alias",
            False,
            "legacy_provenance_hash_alias",
            (
                "Legacy kinase bundle schemas are no longer supported.*"
                "Regenerate the bundle with the current PhosPy version.*"
                "Legacy provenance schemas are no longer supported"
            ),
            id="legacy-provenance-hash-alias",
        ),
        pytest.param(
            "legacy_minimal_total_correction",
            True,
            "legacy_minimal_total_correction",
            (
                "Legacy kinase bundle schemas are no longer supported.*"
                "Regenerate the bundle with the current PhosPy version.*"
                "dataset.metadata.processing_state.total_protein_correction."
                "requires_log_scale is required"
            ),
            id="total-correction-legacy-minimal",
        ),
        pytest.param(
            "missing_total_correction_quantitative_meaning",
            True,
            "missing_total_correction_quantitative_meaning",
            (
                "Legacy kinase bundle schemas are no longer supported.*"
                "Regenerate the bundle with the current PhosPy version.*"
                "dataset.metadata.processing_state.total_protein_correction."
                "quantitative_meaning is required"
            ),
            id="total-correction-missing-quantitative-meaning",
        ),
        pytest.param(
            "missing_total_correction_diagnostics_schema_version",
            True,
            "missing_total_correction_diagnostics_schema_version",
            (
                "Legacy kinase bundle schemas are no longer supported.*"
                "Regenerate the bundle with the current PhosPy version.*"
                "dataset.metadata.processing_state.total_protein_correction."
                "diagnostics.diagnostics_schema_version is required"
            ),
            id="total-correction-missing-diagnostics-schema-version",
        ),
    ],
)
def test_kinase_bundle_contract_rejection_matrix(
    tmp_path: Path,
    scenario_name: str,
    requires_total_correction_payload: bool,
    mutation_kind: str,
    pattern: str,
) -> None:
    # Bundle contract matrix keeps missing/legacy payload failures explicit.
    request = (
        _build_request_with_subtract_log_total(activity=False)
        if requires_total_correction_payload
        else _build_request(activity=False)
    )
    result = KinaseWorkflow().run(request)
    bundle_root = tmp_path / f"kinase_bundle_contract_{scenario_name}"

    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    manifest_path = bundle_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    if mutation_kind == "missing_provenance":
        manifest.pop("provenance", None)
    elif mutation_kind == "null_provenance":
        manifest["provenance"] = None
    elif mutation_kind == "missing_activity_enabled":
        manifest["outputs"]["activity"].pop("enabled", None)
    elif mutation_kind == "legacy_provenance_hash_alias":
        table_payload = manifest["provenance"]["input_tables"][0]
        table_payload["hash_algorithm"] = "sha256"
        table_payload["hash_value"] = table_payload["tolerance_hash_value"]
    elif mutation_kind == "legacy_minimal_total_correction":
        correction_payload = manifest["dataset"]["metadata"]["processing_state"][
            "total_protein_correction"
        ]
        manifest["dataset"]["metadata"]["processing_state"][
            "total_protein_correction"
        ] = {
            "policy": correction_payload["policy"],
            "applied": correction_payload["applied"],
        }
    elif mutation_kind == "missing_total_correction_quantitative_meaning":
        manifest["dataset"]["metadata"]["processing_state"][
            "total_protein_correction"
        ].pop("quantitative_meaning", None)
    elif mutation_kind == "missing_total_correction_diagnostics_schema_version":
        correction_payload = manifest["dataset"]["metadata"]["processing_state"][
            "total_protein_correction"
        ]
        diagnostics_payload = dict(correction_payload["diagnostics"])
        diagnostics_payload.pop("diagnostics_schema_version", None)
        correction_payload["diagnostics"] = diagnostics_payload
    else:
        raise AssertionError(f"Unknown mutation scenario: {mutation_kind}")

    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    with pytest.raises(PhosPyInputError, match=pattern):
        load_kinase_workflow_bundle(bundle_root)


def _activity_result_factory(name: str):
    factories = {
        "weighted_sample": _weighted_sample_activity_result,
        "weighted_condition_summary": _weighted_condition_summary_activity_result,
        "ksea_sample_log_abundance": _ksea_sample_log_abundance_activity_result,
        "ksea_contrast_log_fold_change": _ksea_contrast_log_fold_change_activity_result,
        "ksea_standardised_effect": _ksea_standardised_effect_activity_result,
        "ssgsea_contrast_log_fold_change": _ssgsea_contrast_log_fold_change_activity_result,
        "ssgsea_standardised_effect": _ssgsea_standardised_effect_activity_result,
    }
    return factories[name]


def _weighted_sample_activity_result() -> KinaseActivityResult:
    pred_mat = _activity_pred_mat()
    matrix = _sample_abundance_activity_matrix()
    activity_input = ActivityInputMatrix.sample_level_abundance(
        matrix,
        _assume_owned=True,
    )
    return SimplifiedWeightedSubstrateActivityMethod(
        threshold=0.5,
        min_substrates=2,
        top_n_substrates=2,
    ).run(
        _activity_inputs(
            pred_mat=pred_mat,
            matrix=matrix,
            activity_input=activity_input,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=2,
        )
    )


def _weighted_condition_summary_activity_result() -> KinaseActivityResult:
    pred_mat = _activity_pred_mat()
    matrix = pd.DataFrame(
        {
            "treated_mean": [2.0, 4.0, 6.0],
            "control_mean": [1.0, 3.0, 5.0],
        },
        index=pred_mat.index.copy(),
        dtype=float,
    )
    aggregation_metadata = ActivityAggregationMetadata(
        aggregation_method="mean",
        records=(
            ActivityAggregationRecord(
                profile_id="treated_mean",
                source_profile_ids=("treated_rep1", "treated_rep2"),
            ),
            ActivityAggregationRecord(
                profile_id="control_mean",
                source_profile_ids=("control_rep1", "control_rep2"),
            ),
        ),
    )
    activity_input = ActivityInputMatrix.condition_summary_abundance(
        matrix,
        aggregation_metadata=aggregation_metadata,
        _assume_owned=True,
    )
    return SimplifiedWeightedSubstrateActivityMethod(
        threshold=0.5,
        min_substrates=2,
        top_n_substrates=2,
    ).run(
        _activity_inputs(
            pred_mat=pred_mat,
            matrix=matrix,
            activity_input=activity_input,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=2,
        )
    )


def _ksea_sample_log_abundance_activity_result() -> KinaseActivityResult:
    pred_mat = _activity_pred_mat()
    matrix = np.log2(_sample_abundance_activity_matrix())
    activity_input = ActivityInputMatrix.sample_level_abundance(
        matrix,
        _assume_owned=True,
    )
    return _ksea_activity_result(
        pred_mat=pred_mat,
        matrix=matrix,
        activity_input=activity_input,
    )


def _ksea_contrast_log_fold_change_activity_result() -> KinaseActivityResult:
    pred_mat = _activity_pred_mat()
    matrix = _contrast_activity_matrix()
    activity_input = ActivityInputMatrix.contrast_log_fold_change(
        matrix,
        _assume_owned=True,
    )
    return _ksea_activity_result(
        pred_mat=pred_mat,
        matrix=matrix,
        activity_input=activity_input,
    )


def _ksea_standardised_effect_activity_result() -> KinaseActivityResult:
    pred_mat = _activity_pred_mat()
    matrix = _standardised_effect_activity_matrix()
    activity_input = ActivityInputMatrix.standardised_effect(
        matrix,
        _assume_owned=True,
    )
    return _ksea_activity_result(
        pred_mat=pred_mat,
        matrix=matrix,
        activity_input=activity_input,
    )


def _ksea_activity_result(
    *,
    pred_mat: pd.DataFrame,
    matrix: pd.DataFrame,
    activity_input: ActivityInputMatrix,
) -> KinaseActivityResult:
    return KseaZScoreActivityMethod(
        evidence_threshold=0.5,
        min_substrates=2,
        adjust_p_values=True,
    ).run(
        _activity_inputs(
            pred_mat=pred_mat,
            matrix=matrix,
            activity_input=activity_input,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=1,
        )
    )


def _ssgsea_contrast_log_fold_change_activity_result() -> KinaseActivityResult:
    matrix = _contrast_activity_matrix()
    activity_input = ActivityInputMatrix.contrast_log_fold_change(
        matrix,
        _assume_owned=True,
    )
    return _ssgsea_activity_result(activity_input)


def _ssgsea_standardised_effect_activity_result() -> KinaseActivityResult:
    matrix = _standardised_effect_activity_matrix()
    activity_input = ActivityInputMatrix.standardised_effect(
        matrix,
        _assume_owned=True,
    )
    return _ssgsea_activity_result(activity_input)


def _ssgsea_activity_result(
    activity_input: ActivityInputMatrix,
) -> KinaseActivityResult:
    return SsgseaSubstrateEnrichmentActivityMethod(
        min_substrates=2,
        permutation_count=0,
        adjust_p_values=True,
    ).run(
        activity_input=activity_input,
        kinase_substrate_membership=_ssgsea_membership(),
    )


def _activity_inputs(
    *,
    pred_mat: pd.DataFrame,
    matrix: pd.DataFrame,
    activity_input: ActivityInputMatrix,
    threshold: float,
    min_substrates: int,
    top_n_substrates: int,
) -> KinaseActivityInputs:
    return KinaseActivityInputs(
        pred_mat=pred_mat,
        phospho_matrix=matrix,
        threshold=threshold,
        min_substrates=min_substrates,
        top_n_substrates=top_n_substrates,
        overlap_summary=PredMatOverlapSummary(
            overlap_count=int(pred_mat.index.intersection(matrix.index).size),
            pred_mat_rows=int(pred_mat.index.size),
            phospho_rows=int(matrix.index.size),
        ),
        activity_input=activity_input,
    )


def _activity_pred_mat() -> pd.DataFrame:
    site_index = _activity_site_index()
    pred_mat = pd.DataFrame(
        {
            "K1": [0.9, 0.8, 0.1],
            "K2": [0.1, 0.85, 0.9],
        },
        index=site_index,
        dtype=float,
    )
    pred_mat.columns.name = "kinase"
    return pred_mat


def _sample_abundance_activity_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 4.0],
            "sample_b": [2.0, 3.0, 8.0],
        },
        index=_activity_pred_mat().index.copy(),
        dtype=float,
    )


def _contrast_activity_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "treated_vs_control": [1.25, 0.75, -0.5],
            "drug_vs_vehicle": [-0.75, 0.5, 1.5],
        },
        index=_activity_pred_mat().index.copy(),
        dtype=float,
    )


def _standardised_effect_activity_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "effect_a": [1.1, 0.4, -0.8],
            "effect_b": [-0.6, 0.7, 1.3],
        },
        index=_activity_pred_mat().index.copy(),
        dtype=float,
    )


def _ssgsea_membership() -> pd.DataFrame:
    site_ids = tuple(str(value) for value in _activity_site_index())
    return pd.DataFrame(
        {
            "kinase": ["K1", "K1", "K2", "K2"],
            "substrate_site": [site_ids[0], site_ids[1], site_ids[1], site_ids[2]],
        }
    )


def _activity_site_index() -> pd.Index:
    return site_key_index_from_display_ids(
        ("S1;S1;", "S2;S2;", "S3;S3;"),
        protein_namespace="gene_symbol",
    )


def _save_semantic_activity_bundle(
    tmp_path: Path,
    *,
    bundle_name: str,
    activity_result: KinaseActivityResult,
    method: str,
    resolved_scale: str,
    resolved_meaning: str,
) -> tuple[Path, dict[str, object]]:
    base_request = _build_request(activity=False)
    base_result = KinaseWorkflow().run(base_request)
    result = _replace_activity_result_with_semantic_provenance(
        base_result,
        activity_result=activity_result,
        method=method,
        resolved_scale=resolved_scale,
        resolved_meaning=resolved_meaning,
    )
    request = replace(
        base_request,
        activity_config=_activity_config_for_method(method),
    )
    bundle_root = tmp_path / bundle_name
    save_kinase_workflow_bundle(
        result,
        bundle_root,
        config_snapshot=KinaseWorkflowConfigSnapshot.from_request(request),
    )
    manifest = json.loads((bundle_root / "manifest.json").read_text(encoding="utf-8"))
    return bundle_root, manifest


def _replace_activity_result_with_semantic_provenance(
    base_result: KinaseWorkflowResult,
    *,
    activity_result: KinaseActivityResult,
    method: str,
    resolved_scale: str,
    resolved_meaning: str,
) -> KinaseWorkflowResult:
    assert base_result.provenance is not None
    return KinaseWorkflowResult._from_owned(
        dataset=base_result.dataset,
        references=base_result.references,
        scoring_result=base_result.scoring_result,
        prediction_result=base_result.prediction_result,
        eligibility_report=base_result.eligibility_report,
        site_attrition_summary=base_result.site_attrition_summary,
        attrition_provenance=base_result.attrition_provenance,
        activity_result=activity_result,
        provenance=_provenance_with_activity_semantics(
            base_result.provenance,
            activity_result=activity_result,
            method=method,
            resolved_scale=resolved_scale,
            resolved_meaning=resolved_meaning,
        ),
        substrate_contributions=base_result.substrate_contributions,
        caveats=base_result.caveats,
    )


def _provenance_with_activity_semantics(
    provenance: RunProvenance,
    *,
    activity_result: KinaseActivityResult,
    method: str,
    resolved_scale: str,
    resolved_meaning: str,
) -> RunProvenance:
    workflow_parameters = dict(provenance.workflow_parameters)
    activity_config_raw = workflow_parameters.get("activity_config")
    activity_config = (
        dict(activity_config_raw) if isinstance(activity_config_raw, dict) else {}
    )
    activity_config["method"] = method
    activity_config["method_input_contract"] = ResolvedMethodQuantitativeInputContract(
        contract=kinase_activity_method_quantitative_input_contract(method),
        resolved_scale=resolved_scale,
        resolved_meaning=resolved_meaning,
        resolved_activity_profile_axis=(activity_result.input_semantics.profile_axis),
        resolved_activity_quantitative_semantics=(
            activity_result.input_semantics.quantitative_semantics
        ),
        enforcement_context="kinase bundle semantic round-trip test",
    ).to_payload()
    activity_config["activity_method"] = activity_result.activity_method.to_payload()
    activity_config["activity_method_summary"] = (
        None
        if activity_result.method_summary is None
        else activity_result.method_summary.to_payload()
    )
    workflow_parameters["activity_config"] = activity_config
    return replace(provenance, workflow_parameters=workflow_parameters)


def _activity_config_for_method(method: str) -> KinaseActivityConfig:
    if method == KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY:
        return KinaseActivityConfig(
            enabled=True,
            method=method,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=2,
            ksea_min_substrates=2,
            ksea_evidence_threshold=0.5,
            ssgsea_min_substrates=2,
        )
    if method == KINASE_ACTIVITY_METHOD_KSEA_ZSCORE:
        return KinaseActivityConfig(
            enabled=True,
            method=method,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=2,
            ksea_min_substrates=2,
            ksea_evidence_threshold=0.5,
            ssgsea_min_substrates=2,
        )
    if method == KINASE_ACTIVITY_METHOD_SSGSEA_SUBSTRATE_ENRICHMENT:
        return KinaseActivityConfig(
            enabled=True,
            method=method,
            threshold=0.5,
            min_substrates=2,
            top_n_substrates=2,
            ksea_min_substrates=2,
            ksea_evidence_threshold=0.5,
            ssgsea_min_substrates=2,
            ssgsea_permutations=0,
            ssgsea_random_seed=0,
        )
    raise AssertionError(f"Unsupported test activity method: {method}")


def _assert_activity_result_semantics_round_tripped(
    loaded: KinaseActivityResult,
    original: KinaseActivityResult,
) -> None:
    assert loaded.input_semantics == original.input_semantics
    assert loaded.profile_metadata == original.profile_metadata
    assert tuple(str(column) for column in loaded.activity_matrix.columns) == (
        original.profile_metadata.profile_ids
    )
    assert loaded.activity_matrix.columns.name == _expected_activity_axis_name(
        original.input_semantics.profile_axis
    )
    aggregation = loaded.profile_metadata.aggregation_metadata
    expected_aggregation = original.profile_metadata.aggregation_metadata
    if expected_aggregation is None:
        assert aggregation is None
    else:
        assert aggregation is not None
        assert aggregation.aggregation_method == expected_aggregation.aggregation_method
        assert aggregation.records == expected_aggregation.records
        assert tuple(
            record.source_profile_ids for record in aggregation.records
        ) == tuple(record.source_profile_ids for record in expected_aggregation.records)

    pd.testing.assert_frame_equal(
        loaded.activity_matrix,
        original.activity_matrix,
        check_dtype=False,
        check_index_type=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        loaded.thresholded_substrate_mean_activity,
        original.thresholded_substrate_mean_activity,
        check_dtype=False,
        check_index_type=False,
        check_column_type=False,
        check_names=False,
    )
    pd.testing.assert_series_equal(
        loaded.thresholded_substrate_counts,
        original.thresholded_substrate_counts,
        check_dtype=False,
        check_index_type=False,
        check_names=False,
    )
    _assert_optional_frame_equal(
        loaded.activity_substrate_counts,
        original.activity_substrate_counts,
    )
    pd.testing.assert_series_equal(
        loaded.target_counts,
        original.target_counts,
        check_dtype=False,
        check_index_type=False,
        check_names=False,
    )
    pd.testing.assert_frame_equal(
        loaded.target_table,
        original.target_table,
        check_dtype=False,
        check_index_type=False,
        check_names=False,
    )
    _assert_optional_frame_equal(loaded.statistics_table, original.statistics_table)
    assert loaded.activity_method == original.activity_method
    assert loaded.method_summary == original.method_summary


def _expected_activity_axis_name(axis: ActivityProfileAxis | str) -> str | None:
    axis_value = axis.value if isinstance(axis, ActivityProfileAxis) else str(axis)
    if axis_value == ActivityProfileAxis.SAMPLE.value:
        return None
    if axis_value == ActivityProfileAxis.CONDITION_SUMMARY.value:
        return "condition"
    return "profile_id"


def _assert_activity_provenance_matches_result(
    result: KinaseWorkflowResult,
) -> None:
    assert result.provenance is not None
    assert result.activity_result is not None
    activity_config = result.provenance.workflow_parameters["activity_config"]
    assert isinstance(activity_config, dict)
    method_contract = activity_config["method_input_contract"]
    assert isinstance(method_contract, dict)
    profile_axis = result.activity_result.input_semantics.profile_axis
    quantity = result.activity_result.input_semantics.quantitative_semantics
    profile_axis_value = (
        profile_axis.value
        if isinstance(profile_axis, ActivityProfileAxis)
        else str(profile_axis)
    )
    quantity_value = (
        quantity.value
        if isinstance(quantity, ActivityQuantitativeSemantics)
        else str(quantity)
    )
    assert method_contract["resolved_activity_profile_axis"] == (profile_axis_value)
    assert method_contract["resolved_activity_quantitative_semantics"] == (
        quantity_value
    )


def _mutate_activity_manifest(
    manifest: dict[str, object],
    *,
    bundle_root: Path,
    mutation: str,
) -> None:
    activity_payload = manifest["outputs"]["activity"]
    assert isinstance(activity_payload, dict)
    if mutation == "missing_input_semantics":
        activity_payload.pop("input_semantics", None)
        return
    if mutation == "missing_profile_metadata":
        activity_payload.pop("profile_metadata", None)
        return
    if mutation == "semantic_metadata_supplied_while_disabled":
        activity_payload["input_semantics"] = {
            "profile_axis": "sample",
            "quantitative_semantics": "sample_level_abundance",
        }
        activity_payload["profile_metadata"] = {
            "axis": "sample",
            "profile_ids": ["sample_a"],
            "sample_ids": ["sample_a"],
            "condition_ids": [],
            "contrast_ids": [],
            "aggregation_metadata": None,
        }
        return

    input_semantics = activity_payload["input_semantics"]
    profile_metadata = activity_payload["profile_metadata"]
    assert isinstance(input_semantics, dict)
    assert isinstance(profile_metadata, dict)

    if mutation == "axis_quantitative_semantics_mismatch":
        input_semantics["profile_axis"] = "sample"
        input_semantics["quantitative_semantics"] = "contrast_log_fold_change"
    elif mutation == "profile_ids_reordered":
        reordered = list(reversed(profile_metadata["profile_ids"]))
        profile_metadata["profile_ids"] = reordered
        profile_metadata["sample_ids"] = reordered
    elif mutation == "unknown_profile_ids":
        profile_metadata["profile_ids"] = ["unknown_a", "unknown_b"]
        profile_metadata["sample_ids"] = ["unknown_a", "unknown_b"]
    elif mutation == "contradictory_condition_ids":
        profile_metadata["condition_ids"] = ["condition_that_should_not_exist"]
    elif mutation == "missing_condition_summary_aggregation":
        profile_metadata["aggregation_metadata"] = None
    elif mutation == "duplicate_aggregation_records":
        aggregation = profile_metadata["aggregation_metadata"]
        assert isinstance(aggregation, dict)
        records = aggregation["records"]
        assert isinstance(records, list)
        assert isinstance(records[0], dict)
        assert isinstance(records[1], dict)
        records[1]["profile_id"] = records[0]["profile_id"]
    elif mutation == "provenance_manifest_axis_conflict":
        method_contract = _activity_method_contract_payload(manifest)
        method_contract["resolved_activity_profile_axis"] = "contrast"
    elif mutation == "provenance_manifest_quantitative_conflict":
        method_contract = _activity_method_contract_payload(manifest)
        method_contract["resolved_activity_quantitative_semantics"] = (
            "contrast_log_fold_change"
        )
    elif mutation == "table_tamper_with_digest_refresh":
        tables = activity_payload["tables"]
        assert isinstance(tables, dict)
        entry = tables["weighted_activity"]
        assert isinstance(entry, dict)
        table_path = bundle_root / str(entry["path"])
        table = pd.read_csv(table_path, index_col=0)
        table = table.loc[:, list(reversed(table.columns))]
        table.to_csv(table_path)
        _refresh_manifest_table_entry(entry, bundle_root=bundle_root, table=table)
    else:
        raise AssertionError(f"Unknown mutation scenario: {mutation}")


def _activity_method_contract_payload(
    manifest: dict[str, object],
) -> dict[str, object]:
    provenance = manifest["provenance"]
    assert isinstance(provenance, dict)
    workflow_parameters = provenance["workflow_parameters"]
    assert isinstance(workflow_parameters, dict)
    activity_config = workflow_parameters["activity_config"]
    assert isinstance(activity_config, dict)
    method_contract = activity_config["method_input_contract"]
    assert isinstance(method_contract, dict)
    return method_contract


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
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
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
            input_intensity_scale="linear",
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
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
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


def _build_request_with_subtract_log_total_and_uncorrected_rows(
    *,
    activity: bool,
) -> KinaseWorkflowRequest:
    phospho = pd.DataFrame(
        {
            "sample_a": [15.0, 7.0],
            "sample_b": [31.0, 15.0],
        },
        index=["MAPK14;Y182;", "AKT1;T308;"],
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
            "protein_id": ["MAPK14", "AKT1"],
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
            "sample_a": [3.0],
            "sample_b": [7.0],
        },
        index=["MAPK14"],
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
                    policy="subtract_log_total",
                    identity=DatasetTotalProteinCorrectionIdentityConfig(
                        unmatched_policy=DATASET_TOTAL_PROTEIN_CORRECTION_UNMATCHED_POLICY_ALLOW_UNCORRECTED
                    ),
                ),
            ),
        )
    )
    references = ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6"],
                "substrate_site": [
                    "MAPK14;Y182;",
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
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            allow_mixed_total_protein_quantitative_meaning=True,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
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
            "protein_id": ["MAPK14", "AKT1", "GSK3B"],
            "site": ["Y182", "T308", "S9"],
            "site_sequence": [
                "LDFGLARHTDDEMTGYVATRWYRAPEIMLNW",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95] * phospho.shape[0],
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
        scoring_config=KinaseScoringConfig(
            reliability_profile="custom",
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
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
    assert left.attrition_provenance == right.attrition_provenance
    assert left.caveats == right.caveats

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
    assert left.activity_result.input_semantics == right.activity_result.input_semantics
    assert (
        left.activity_result.profile_metadata == right.activity_result.profile_metadata
    )
    pd.testing.assert_frame_equal(
        left.activity_result.activity_matrix,
        right.activity_result.activity_matrix,
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
        check_column_type=False,
    )
    pd.testing.assert_series_equal(
        left.activity_result.thresholded_substrate_counts,
        right.activity_result.thresholded_substrate_counts,
        check_dtype=False,
        check_names=False,
        check_index_type=False,
    )
    _assert_optional_frame_equal(
        left.activity_result.activity_substrate_counts,
        right.activity_result.activity_substrate_counts,
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
    _assert_optional_frame_equal(
        left.activity_result.statistics_table,
        right.activity_result.statistics_table,
    )
    assert left.activity_result.method_summary == right.activity_result.method_summary


def _assert_optional_frame_equal(left, right) -> None:
    if left is None or right is None:
        assert left is right
        return
    pd.testing.assert_frame_equal(
        left,
        right,
        check_dtype=False,
        check_names=False,
        check_column_type=False,
    )
