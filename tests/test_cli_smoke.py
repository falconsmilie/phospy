from __future__ import annotations

import sys
from pathlib import Path

import pytest

import phospy.cli as cli
from phospy.api import DatasetLoadOptions, KinaseActivityConfig
from phospy.internal.constants import (
    CORE_OUTPUT_ARTIFACT_BASENAMES,
    KINASE_OUTPUT_FILENAMES,
)
from phospy.preprocessing import CorePreprocessingConfig


def _run_cli(monkeypatch, *args: str) -> None:
    monkeypatch.setattr(sys, "argv", ["phospy", *args])
    cli.main()


def _example_data_dir() -> Path:
    return Path(__file__).resolve().parents[1] / "examples" / "data"


def test_build_cli_configs_maps_scalar_args_to_typed_configs() -> None:
    parser = cli.build_parser()
    args = parser.parse_args(
        [
            "--total",
            "total.tsv",
            "--phospho",
            "phospho.tsv",
            "--outdir",
            "out",
            "--phospho-encoding",
            "utf-16le",
            "--localization-threshold",
            "0.8",
            "--min-observed",
            "5",
            "--total-sentinel",
            "11.0",
            "--phospho-sentinel",
            "13.0",
            "--max-unmatched-fraction",
            "0.25",
            "--kinase-activity-threshold",
            "0.7",
            "--kinase-activity-min-substrates",
            "4",
            "--kinase-activity-top-n-substrates",
            "9",
        ]
    )

    dataset_options, preprocessing_config, activity_config = cli.build_cli_configs(args)

    assert dataset_options == DatasetLoadOptions(phospho_encoding="utf-16le")
    assert preprocessing_config == CorePreprocessingConfig(
        localization_threshold=0.8,
        min_observed=5,
        total_sentinel=11.0,
        phospho_sentinel=13.0,
        max_unmatched_fraction=0.25,
    )
    assert activity_config == KinaseActivityConfig(
        threshold=0.7,
        min_substrates=4,
        top_n_substrates=9,
    )


def test_cli_main_runs_end_to_end(tmp_path, monkeypatch) -> None:
    data_dir = _example_data_dir()
    outdir = tmp_path / "cli-output"
    _run_cli(
        monkeypatch,
        "--total",
        str(data_dir / "total.tsv"),
        "--phospho",
        str(data_dir / "phospho.tsv"),
        "--pred-mat",
        str(data_dir / "predMat.csv"),
        "--outdir",
        str(outdir),
        "--phospho-encoding",
        "utf-16le",
        "--max-unmatched-fraction",
        "0.1",
        "--kinase-activity-min-substrates",
        "1",
        "--kinase-activity-top-n-substrates",
        "1",
    )

    expected = {
        *(f"{basename}.csv" for basename in CORE_OUTPUT_ARTIFACT_BASENAMES),
        *KINASE_OUTPUT_FILENAMES,
    }
    assert expected.issubset({path.name for path in outdir.iterdir()})


def test_cli_main_reports_missing_paths_without_traceback(
    tmp_path,
    monkeypatch,
    capsys,
) -> None:
    data_dir = _example_data_dir()
    outdir = tmp_path / "cli-output"

    with pytest.raises(SystemExit) as exc_info:
        _run_cli(
            monkeypatch,
            "--total",
            str(tmp_path / "missing_total.tsv"),
            "--phospho",
            str(data_dir / "phospho.tsv"),
            "--outdir",
            str(outdir),
        )

    assert exc_info.value.code == cli.CLI_EXIT_USER_ERROR
    stderr = capsys.readouterr().err
    assert "Path does not exist" in stderr
    assert "Traceback" not in stderr


def test_cli_main_reports_malformed_files_without_traceback(
    tmp_path,
    monkeypatch,
    capsys,
) -> None:
    data_dir = _example_data_dir()
    outdir = tmp_path / "cli-output"
    malformed_phospho = tmp_path / "phospho-malformed.tsv"
    malformed_phospho.write_bytes(b"\xff\xfe\x00\x00")

    with pytest.raises(SystemExit) as exc_info:
        _run_cli(
            monkeypatch,
            "--total",
            str(data_dir / "total.tsv"),
            "--phospho",
            str(malformed_phospho),
            "--outdir",
            str(outdir),
        )

    assert exc_info.value.code == cli.CLI_EXIT_USER_ERROR
    stderr = capsys.readouterr().err
    assert "unable to read file" in stderr
    assert "Traceback" not in stderr


def test_cli_main_reports_schema_failures_without_traceback(
    tmp_path,
    monkeypatch,
    capsys,
) -> None:
    data_dir = _example_data_dir()
    outdir = tmp_path / "cli-output"
    invalid_total = tmp_path / "total-invalid.tsv"
    invalid_total.write_text(
        (
            "protein\tgroup1\tgroup2\tgroup3\tgroup4\tgroup5\tgroup6\n"
            "MAPK1\t1\t1\t1\t1\t1\t1\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(SystemExit) as exc_info:
        _run_cli(
            monkeypatch,
            "--total",
            str(invalid_total),
            "--phospho",
            str(data_dir / "phospho.tsv"),
            "--outdir",
            str(outdir),
        )

    assert exc_info.value.code == cli.CLI_EXIT_USER_ERROR
    stderr = capsys.readouterr().err
    assert "missing required columns: genes" in stderr
    assert "Traceback" not in stderr


def test_cli_main_reports_runtime_input_incompatibility_without_traceback(
    tmp_path,
    monkeypatch,
    capsys,
) -> None:
    data_dir = _example_data_dir()
    outdir = tmp_path / "cli-output"
    incompatible_pred_mat = tmp_path / "predMat-incompatible.csv"
    incompatible_pred_mat.write_text(
        ",KINASE_A\nNO_SHARED_SITE,0.9\n",
        encoding="utf-8",
    )

    with pytest.raises(SystemExit) as exc_info:
        _run_cli(
            monkeypatch,
            "--total",
            str(data_dir / "total.tsv"),
            "--phospho",
            str(data_dir / "phospho.tsv"),
            "--phospho-encoding",
            "utf-16le",
            "--pred-mat",
            str(incompatible_pred_mat),
            "--outdir",
            str(outdir),
        )

    assert exc_info.value.code == cli.CLI_EXIT_USER_ERROR
    stderr = capsys.readouterr().err
    assert "Pipeline runtime compatibility failed after preprocessing" in stderr
    assert "no overlapping phosphosite IDs" in stderr
    assert "Traceback" not in stderr
