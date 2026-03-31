from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pytest

from phospy import KinaseActivityAnalyzer
from phospy.core_processing import CorePreprocessingConfig
from phospy.dataset import PhosphoDataset
from phospy.publishing import OutputPublisher, RunManifestWriter

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
    return dataset.process_core()


def test_output_publisher_replaces_existing_directory_atomically(
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
            raise RuntimeError("boom")
        source.replace(target)

    monkeypatch.setattr(
        OutputPublisher,
        "_replace_directory",
        staticmethod(replace_with_failure),
    )

    with pytest.raises(RuntimeError, match="boom"):
        OutputPublisher().publish(staging_dir=staging_dir, target_dir=target_dir)

    assert target_dir.exists()
    assert (target_dir / "old.txt").read_text(encoding="utf-8") == "old"
    assert staging_dir.exists()
    assert (staging_dir / "new.txt").read_text(encoding="utf-8") == "new"
    assert not list(tmp_path.glob(".published.backup-*"))


def test_run_manifest_writer_serializes_expected_metadata(tmp_path: Path) -> None:
    core = _build_core_result()
    kinase_activity = KinaseActivityAnalyzer().analyze(
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

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
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
