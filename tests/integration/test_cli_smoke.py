from __future__ import annotations

from pathlib import Path

import pytest

from phospy.cli import main as cli_main
from tests.support.rewrite_fixture_data import load_rat_l6_phospho, site_metadata_for

pytestmark = pytest.mark.integration


def test_cli_dataset_build_from_files_writes_expected_outputs(tmp_path: Path) -> None:
    phospho = load_rat_l6_phospho().head(64).copy(deep=True)
    site_metadata = site_metadata_for(phospho)
    phospho_path = tmp_path / "phospho.csv"
    site_metadata_path = tmp_path / "site_metadata.csv"
    outdir = tmp_path / "out"
    phospho.to_csv(phospho_path)
    site_metadata.to_csv(site_metadata_path)

    exit_code = cli_main(
        [
            "dataset-build",
            "--phospho",
            str(phospho_path),
            "--site-metadata",
            str(site_metadata_path),
            "--organism",
            "rat",
            "--outdir",
            str(outdir),
        ]
    )

    assert exit_code == 0
    assert (outdir / "dataset" / "phospho.csv").exists()
    assert (outdir / "dataset" / "site_metadata.csv").exists()
    assert (outdir / "dataset" / "manifest.json").exists()


def test_cli_simple_kinase_from_files_writes_supported_lane_outputs(
    tmp_path: Path,
) -> None:
    phospho = load_rat_l6_phospho().head(260).copy(deep=True)
    site_metadata = site_metadata_for(phospho)
    phospho_path = tmp_path / "phospho.csv"
    site_metadata_path = tmp_path / "site_metadata.csv"
    outdir = tmp_path / "out"
    phospho.to_csv(phospho_path)
    site_metadata.to_csv(site_metadata_path)

    exit_code = cli_main(
        [
            "simple-kinase",
            "--phospho",
            str(phospho_path),
            "--site-metadata",
            str(site_metadata_path),
            "--organism",
            "rat",
            "--reference",
            "auto",
            "--prediction-top-k",
            "6",
            "--prediction-ensemble-size",
            "8",
            "--outdir",
            str(outdir),
        ]
    )

    assert exit_code == 0
    assert (outdir / "dataset" / "phospho.csv").exists()
    assert (outdir / "simple_kinase" / "scoring" / "profile_scores.csv").exists()
    assert (outdir / "simple_kinase" / "prediction" / "pred_mat.csv").exists()
    assert (
        outdir / "simple_kinase" / "references" / "kinase_substrate_map.csv"
    ).exists()
    assert (outdir / "simple_kinase" / "activity" / "activity_scores.csv").exists()
    assert (outdir / "simple_kinase" / "manifest.json").exists()


def test_cli_reports_rewrite_taxonomy_for_input_format_errors(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    phospho = load_rat_l6_phospho().head(8).copy(deep=True)
    site_metadata = site_metadata_for(phospho)
    phospho_path = tmp_path / "phospho.invalid"
    site_metadata_path = tmp_path / "site_metadata.csv"
    phospho.to_csv(phospho_path)
    site_metadata.to_csv(site_metadata_path)

    exit_code = cli_main(
        [
            "dataset-build",
            "--phospho",
            str(phospho_path),
            "--site-metadata",
            str(site_metadata_path),
            "--organism",
            "rat",
            "--outdir",
            str(tmp_path / "out"),
        ]
    )
    captured = capsys.readouterr()

    assert exit_code == 2
    assert "UnsupportedInputFormatError" in captured.err


def test_cli_simple_kinase_rejects_non_rat_bundled_resolution(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    phospho = load_rat_l6_phospho().head(32).copy(deep=True)
    site_metadata = site_metadata_for(phospho)
    phospho_path = tmp_path / "phospho.csv"
    site_metadata_path = tmp_path / "site_metadata.csv"
    outdir = tmp_path / "out"
    phospho.to_csv(phospho_path)
    site_metadata.to_csv(site_metadata_path)

    exit_code = cli_main(
        [
            "simple-kinase",
            "--phospho",
            str(phospho_path),
            "--site-metadata",
            str(site_metadata_path),
            "--organism",
            "human",
            "--reference",
            "auto",
            "--outdir",
            str(outdir),
        ]
    )
    captured = capsys.readouterr()

    assert exit_code == 2
    assert "UnsupportedOrganismError" in captured.err
    assert "supported bundled organisms: rat" in captured.err
