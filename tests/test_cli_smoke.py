from __future__ import annotations

import sys
from pathlib import Path

import phospy.cli as cli
from phospy.api import DatasetLoadOptions, KinaseActivityConfig
from phospy.internal.constants import (
    CORE_OUTPUT_ARTIFACT_BASENAMES,
    KINASE_OUTPUT_FILENAMES,
)
from phospy.preprocessing import CorePreprocessingConfig


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
    repo_root = Path(__file__).resolve().parents[1]
    outdir = tmp_path / "cli-output"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "phospy",
            "--total",
            str(repo_root / "examples" / "data" / "total.tsv"),
            "--phospho",
            str(repo_root / "examples" / "data" / "phospho.tsv"),
            "--pred-mat",
            str(repo_root / "examples" / "data" / "predMat.csv"),
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
        ],
    )

    cli.main()

    expected = {
        *(f"{basename}.csv" for basename in CORE_OUTPUT_ARTIFACT_BASENAMES),
        *KINASE_OUTPUT_FILENAMES,
    }
    assert expected.issubset({path.name for path in outdir.iterdir()})
