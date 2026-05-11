"""Typed command objects for CLI execution orchestration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.api.requests import DatasetBuildRequest
from phospy.references.models import ReferencePreset


@dataclass(frozen=True, slots=True)
class OutputTarget:
    """Filesystem publishing target for CLI commands."""

    outdir: Path
    output_format: str


@dataclass(frozen=True, slots=True)
class DatasetBuildCommand:
    """Command payload for dataset-build execution."""

    request: DatasetBuildRequest
    output: OutputTarget


@dataclass(frozen=True, slots=True)
class KinaseCommand:
    """Command payload for kinase workflow execution."""

    dataset_request: DatasetBuildRequest
    references: ReferencePreset
    scoring_config: KinaseScoringConfig
    prediction_config: KinasePredictionConfig
    activity_config: KinaseActivityConfig | None
    output: OutputTarget


@dataclass(frozen=True, slots=True)
class SignalomeCommand:
    """Command payload for signalome workflow execution."""

    kinase: KinaseCommand
    signalome_config: SignalomeConfig
    output: OutputTarget


CliCommand = DatasetBuildCommand | KinaseCommand | SignalomeCommand
