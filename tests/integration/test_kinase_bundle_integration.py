from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyDatasetBuilder,
    DatasetBuildRequest,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
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
    _assert_kinase_result_equal(loaded.result, result)


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
        "combined_scores": None,
        "motif_scores": None,
        "profile_scores": "scoring/profile_scores.csv",
        "weights": None,
    }
    assert manifest["outputs"]["prediction"]["tables"] == {
        "pred_mat": "prediction/pred_mat.csv",
        "substrate_list": "prediction/substrate_list.csv",
    }
    assert manifest["outputs"]["activity"] == {
        "enabled": True,
        "tables": {"activity_scores": "activity/activity_scores.csv"},
    }


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
        "tables": {"activity_scores": None},
    }


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
        scoring_config=KinaseScoringConfig(min_substrates=1),
        prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=2),
        activity_config=(
            KinaseActivityConfig(enabled=True, threshold=0.5) if activity else None
        ),
    )


def _assert_kinase_result_equal(left, right) -> None:
    assert left.dataset.organism == right.dataset.organism
    assert left.dataset.transformation_state == right.dataset.transformation_state
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
        left.scoring_result.combined_scores,
        right.scoring_result.combined_scores,
    )
    _assert_optional_frame_equal(
        left.scoring_result.weights,
        right.scoring_result.weights,
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
        left.activity_result.activity_scores,
        right.activity_result.activity_scores,
        check_dtype=False,
        check_names=False,
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
