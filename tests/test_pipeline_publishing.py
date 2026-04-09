from __future__ import annotations

import json
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError
from pathlib import Path

import pandas as pd
import pytest

from phospy import KinaseActivityAnalyzer, PhosRPipeline, PredMatResult
from phospy.constants import RUN_MANIFEST_FILENAME
from phospy.core_processing import CorePreprocessingConfig
from phospy.dataset import PhosphoDataset
from phospy.pipeline import _PipelineRequestLoader
from phospy.publishing import OutputPublisher, RunManifestWriter, package_version
from phospy.validation.requests import CorePipelineRequest
from phospy.validation.schemas import PredMatSchema

EXAMPLE_COMPARISONS = (("group1", "group4"), ("group2", "group5"), ("group3", "group6"))


def make_total_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "genes": ["Prkaca", "Prkaca", "Btk", "Lyn"],
            "group1": [1.0, 5.0, 2.0, 3.0],
            "group2": [1.0, 5.0, 2.0, 3.0],
            "group3": [1.0, 5.0, 2.0, 3.0],
            "group4": [1.0, 5.0, 2.0, 3.0],
            "group5": [1.0, 5.0, 2.0, 3.0],
            "group6": [1.0, 5.0, 2.0, 3.0],
        }
    )


def make_phospho_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "uid": ["u1", "u2", "u3", "u4"],
            "gene_names": ["PRKACA", "BTK", "LYN", "PRKACA"],
            "gene_p_site": ["PRKACA_S339", "BTK_Y551", "LYN_Y397", "PRKACA_S339"],
            "localization_prob": [0.95, 0.95, 0.95, 0.95],
            "centralized_sequence": ["AAAAAA", "BBBBBB", "CCCCCC", "DDDDDD"],
            "p_group1": [8.0, 6.0, 7.0, 9.0],
            "p_group2": [7.0, 5.0, 6.0, 8.0],
            "p_group3": [6.0, 4.0, 5.0, 7.0],
            "p_group4": [5.0, 3.0, 4.0, 6.0],
            "p_group5": [4.0, 2.0, 3.0, 5.0],
            "p_group6": [3.0, 1.0, 2.0, 4.0],
        }
    )


def make_pred_mat() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "PRKACA": [0.9, 0.8, 0.7],
            "BTK": [0.2, 0.85, 0.75],
        },
        index=["PRKACA;S339;", "BTK;Y551;", "LYN;Y397;"],
    )


def _build_core_result():
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    return dataset.preprocessing.run()


def test_output_publisher_replaces_existing_directory_via_crash_recoverable_swap(
    tmp_path: Path,
) -> None:
    target_dir = tmp_path / "published"
    staging_dir = tmp_path / "staging"
    target_dir.mkdir()
    staging_dir.mkdir()
    (target_dir / "old.txt").write_text("old", encoding="utf-8")
    (staging_dir / "new.txt").write_text("new", encoding="utf-8")

    OutputPublisher().publish(staging_dir=staging_dir, target_dir=target_dir)

    assert not staging_dir.exists()
    assert target_dir.exists()
    assert not (target_dir / "old.txt").exists()
    assert (target_dir / "new.txt").read_text(encoding="utf-8") == "new"
    assert not list(tmp_path.glob(".published.backup-*"))
    assert not (tmp_path / ".published.publish-state.json").exists()


def test_output_publisher_restores_backup_on_publish_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    target_dir = tmp_path / "published"
    staging_dir = tmp_path / "staging"
    target_dir.mkdir()
    staging_dir.mkdir()
    (target_dir / "old.txt").write_text("old", encoding="utf-8")
    (staging_dir / "new.txt").write_text("new", encoding="utf-8")

    def replace_with_failure(source: Path, target: Path) -> None:
        if source == staging_dir and target == target_dir:
            raise OSError("boom")
        source.replace(target)

    monkeypatch.setattr(
        OutputPublisher,
        "_replace_directory",
        staticmethod(replace_with_failure),
    )

    with pytest.raises(OSError, match="boom"):
        OutputPublisher().publish(staging_dir=staging_dir, target_dir=target_dir)

    assert target_dir.exists()
    assert (target_dir / "old.txt").read_text(encoding="utf-8") == "old"
    assert staging_dir.exists()
    assert (staging_dir / "new.txt").read_text(encoding="utf-8") == "new"
    assert not list(tmp_path.glob(".published.backup-*"))
    assert not (tmp_path / ".published.publish-state.json").exists()


def test_output_publisher_recovers_interrupted_replacement_before_next_publish(
    tmp_path: Path,
) -> None:
    target_dir = tmp_path / "published"
    staging_dir = tmp_path / "staging"
    backup_dir = tmp_path / ".published.backup-stale"
    marker_path = tmp_path / ".published.publish-state.json"

    backup_dir.mkdir()
    staging_dir.mkdir()
    (backup_dir / "old.txt").write_text("old", encoding="utf-8")
    (staging_dir / "new.txt").write_text("new", encoding="utf-8")
    marker_path.write_text(
        json.dumps(
            {
                "target_dir": str(target_dir),
                "backup_dir": str(backup_dir),
                "created_at_utc": "2026-04-02T08:00:00+00:00",
            }
        )
        + "\n",
        encoding="utf-8",
    )

    OutputPublisher().publish(staging_dir=staging_dir, target_dir=target_dir)

    assert target_dir.exists()
    assert (target_dir / "new.txt").read_text(encoding="utf-8") == "new"
    assert not backup_dir.exists()
    assert not marker_path.exists()


def test_output_publisher_cleans_stale_backup_after_completed_replacement(
    tmp_path: Path,
) -> None:
    target_dir = tmp_path / "published"
    staging_dir = tmp_path / "staging"
    backup_dir = tmp_path / ".published.backup-stale"
    marker_path = tmp_path / ".published.publish-state.json"

    target_dir.mkdir()
    backup_dir.mkdir()
    staging_dir.mkdir()
    (target_dir / "current.txt").write_text("current", encoding="utf-8")
    (backup_dir / "old.txt").write_text("old", encoding="utf-8")
    (staging_dir / "new.txt").write_text("new", encoding="utf-8")
    marker_path.write_text(
        json.dumps(
            {
                "target_dir": str(target_dir),
                "backup_dir": str(backup_dir),
                "created_at_utc": "2026-04-02T08:00:00+00:00",
            }
        )
        + "\n",
        encoding="utf-8",
    )

    OutputPublisher().publish(staging_dir=staging_dir, target_dir=target_dir)

    assert target_dir.exists()
    assert (target_dir / "new.txt").read_text(encoding="utf-8") == "new"
    assert not backup_dir.exists()
    assert not marker_path.exists()


def test_output_publisher_quarantines_corrupt_marker_and_continues_publish(
    tmp_path: Path,
) -> None:
    target_dir = tmp_path / "published"
    staging_dir = tmp_path / "staging"
    marker_path = tmp_path / ".published.publish-state.json"

    target_dir.mkdir()
    staging_dir.mkdir()
    (target_dir / "old.txt").write_text("old", encoding="utf-8")
    (staging_dir / "new.txt").write_text("new", encoding="utf-8")
    marker_path.write_text("{not valid json", encoding="utf-8")

    OutputPublisher().publish(staging_dir=staging_dir, target_dir=target_dir)

    assert target_dir.exists()
    assert (target_dir / "new.txt").read_text(encoding="utf-8") == "new"
    assert not (target_dir / "old.txt").exists()
    assert not marker_path.exists()

    quarantined_path = tmp_path / ".published.publish-state.json.corrupt"
    assert quarantined_path.exists()
    assert quarantined_path.read_text(encoding="utf-8") == "{not valid json"


def test_output_publisher_quarantines_malformed_marker_object_and_continues_publish(
    tmp_path: Path,
) -> None:
    target_dir = tmp_path / "published"
    staging_dir = tmp_path / "staging"
    marker_path = tmp_path / ".published.publish-state.json"

    target_dir.mkdir()
    staging_dir.mkdir()
    (target_dir / "old.txt").write_text("old", encoding="utf-8")
    (staging_dir / "new.txt").write_text("new", encoding="utf-8")
    marker_path.write_text(
        json.dumps({"target_dir": str(target_dir)}) + "\n",
        encoding="utf-8",
    )

    OutputPublisher().publish(staging_dir=staging_dir, target_dir=target_dir)

    assert target_dir.exists()
    assert (target_dir / "new.txt").read_text(encoding="utf-8") == "new"
    assert not (target_dir / "old.txt").exists()
    assert not marker_path.exists()

    quarantined_path = tmp_path / ".published.publish-state.json.corrupt"
    assert quarantined_path.exists()
    assert json.loads(quarantined_path.read_text(encoding="utf-8")) == {
        "target_dir": str(target_dir)
    }


def test_pipeline_delegates_manifest_and_publish_to_publishing_layer(
    tmp_path: Path,
) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    pred_path = tmp_path / "predMat.csv"
    outdir = tmp_path / "published"

    make_total_df().to_csv(total_path, sep="	", index=False)
    make_phospho_df().to_csv(phospho_path, sep="	", index=False)
    make_pred_mat().to_csv(pred_path)

    manifest_calls: list[dict[str, object]] = []
    publish_calls: list[dict[str, Path]] = []

    class RecordingManifestWriter(RunManifestWriter):
        def write(self, *, outdir, core, kinase_activity, preprocessing_config) -> None:
            manifest_calls.append(
                {
                    "outdir": outdir,
                    "has_kinase_activity": kinase_activity is not None,
                    "min_observed": preprocessing_config.min_observed,
                    "site_matrix_rows": int(core.site_matrix.matrix.shape[0]),
                }
            )
            super().write(
                outdir=outdir,
                core=core,
                kinase_activity=kinase_activity,
                preprocessing_config=preprocessing_config,
            )

    class RecordingOutputPublisher(OutputPublisher):
        def publish(self, *, staging_dir: Path, target_dir: Path) -> None:
            publish_calls.append(
                {
                    "staging_dir": staging_dir,
                    "target_dir": target_dir,
                }
            )
            super().publish(staging_dir=staging_dir, target_dir=target_dir)

    pipeline = PhosRPipeline.from_files(
        total_path=total_path,
        phospho_path=phospho_path,
        pred_mat_path=pred_path,
    )
    pipeline.manifest_writer = RecordingManifestWriter(
        package_version_resolver=lambda: "1.2.3-test",
        clock=lambda: datetime(2026, 4, 1, 9, 30, tzinfo=timezone.utc),
    )
    pipeline.output_publisher = RecordingOutputPublisher()

    pipeline.run(outdir=outdir)

    assert len(manifest_calls) == 1
    assert manifest_calls[0]["has_kinase_activity"] is True
    assert manifest_calls[0]["min_observed"] == 4
    assert manifest_calls[0]["site_matrix_rows"] == 3

    assert len(publish_calls) == 1
    assert publish_calls[0]["target_dir"] == outdir
    assert publish_calls[0]["staging_dir"].name.startswith(".published.tmp-")

    manifest = json.loads((outdir / RUN_MANIFEST_FILENAME).read_text(encoding="utf-8"))
    assert manifest["package_version"] == "1.2.3-test"
    assert manifest["generated_at_utc"] == "2026-04-01T09:30:00+00:00"


def test_run_manifest_writer_serializes_expected_metadata(tmp_path: Path) -> None:
    core = _build_core_result()
    kinase_activity = KinaseActivityAnalyzer().run(
        pred_mat=make_pred_mat(),
        phospho_matrix=core.site_matrix.matrix,
    )
    outdir = tmp_path / "manifest-output"
    writer = RunManifestWriter(
        package_version_resolver=lambda: "9.9.9",
        clock=lambda: datetime(2026, 3, 31, 12, 0, tzinfo=timezone.utc),
    )

    writer.write(
        outdir=outdir,
        core=core,
        kinase_activity=kinase_activity,
        preprocessing_config=CorePreprocessingConfig(),
    )

    manifest = json.loads((outdir / RUN_MANIFEST_FILENAME).read_text(encoding="utf-8"))
    assert manifest == {
        "status": "success",
        "generated_at_utc": "2026-03-31T12:00:00+00:00",
        "package_version": "9.9.9",
        "has_kinase_activity": True,
        "core_rows": {
            "total_unique": 3,
            "total_filtered": 3,
            "phospho_filtered": 4,
            "phospho_corrected": 4,
            "site_matrix": 3,
        },
        "preprocessing_config": {
            "localization_threshold": 0.75,
            "min_observed": 4,
            "max_unmatched_fraction": 0.0,
            "total_sentinel": 10.0,
            "phospho_sentinel": 12.0,
        },
    }


def test_package_version_returns_unknown_when_distribution_metadata_is_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def missing_distribution(_: str) -> str:
        raise PackageNotFoundError("phospy")

    monkeypatch.setattr("phospy.publishing.version", missing_distribution)

    assert package_version() == "unknown"


def test_package_version_propagates_unexpected_metadata_errors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def blow_up(_: str) -> str:
        raise RuntimeError("metadata backend failed")

    monkeypatch.setattr("phospy.publishing.version", blow_up)

    with pytest.raises(RuntimeError, match="metadata backend failed"):
        package_version()


def test_pipeline_request_loader_builds_dataset_and_config(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    total_path = tmp_path / "total.tsv"
    phospho_path = tmp_path / "phospho.tsv"
    pred_path = tmp_path / "predMat.csv"
    make_total_df().to_csv(total_path, sep="	", index=False)
    make_phospho_df().to_csv(phospho_path, sep="	", index=False)
    make_pred_mat().to_csv(pred_path)

    pred_calls: list[str] = []
    original_pred_validate = PredMatSchema.validate

    def counting_pred_validate(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
        pred_calls.append(context)
        return original_pred_validate(df, context=context)

    monkeypatch.setattr(
        PredMatSchema,
        "validate",
        staticmethod(counting_pred_validate),
    )

    request = CorePipelineRequest.validate_request(
        total_path=total_path,
        phospho_path=phospho_path,
        pred_mat_path=pred_path,
    )

    inputs = _PipelineRequestLoader().load(request)

    assert inputs.pred_mat is not None
    assert inputs.preprocessing_config.min_observed == 4
    assert list(inputs.dataset.total_df_copy.columns) == list(make_total_df().columns)
    assert pred_calls == [f"pred_mat ({pred_path})"]


def test_pipeline_passes_pred_mat_result_to_validation_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    dataset = PhosphoDataset(
        total_df=make_total_df(),
        phospho_df=make_phospho_df(),
        comparisons=EXAMPLE_COMPARISONS,
    )
    pred_mat_result = PredMatResult(make_pred_mat())
    captured: list[object] = []
    original_validate_request = PhosRPipeline.__init__.__globals__[
        "validate_pipeline_construction_request"
    ]

    def capturing_validate_request(*, dataset, pred_mat=None, **kwargs):
        captured.append(pred_mat)
        return original_validate_request(
            dataset=dataset,
            pred_mat=pred_mat,
            **kwargs,
        )

    monkeypatch.setitem(
        PhosRPipeline.__init__.__globals__,
        "validate_pipeline_construction_request",
        capturing_validate_request,
    )

    pipeline = PhosRPipeline(dataset, pred_mat=pred_mat_result)

    assert captured == [pred_mat_result]
    pd.testing.assert_frame_equal(pipeline.pred_mat, pred_mat_result.data_frame)
