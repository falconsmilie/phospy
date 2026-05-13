"""CLI-oriented output publishing adapters."""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.io.cli_commands import OutputTarget
from phospy.io.publishers.workflows import (
    publish_dataset,
    publish_kinase_workflow,
    publish_signalome_workflow,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset


class CliOutputPublisher(Protocol):
    """Publishing contract for command runner collaboration."""

    def publish_dataset(
        self,
        dataset: AnalysisReadyPhosphoDataset,
        output: OutputTarget,
    ) -> dict[str, Path]: ...

    def publish_kinase(
        self,
        result: KinaseWorkflowResult,
        output: OutputTarget,
    ) -> dict[str, Path]: ...

    def publish_signalome(
        self,
        result: SignalomeWorkflowResult,
        output: OutputTarget,
    ) -> dict[str, Path]: ...


class WorkflowOutputPublisher:
    """Default publisher for CLI command execution outputs."""

    def publish_dataset(
        self,
        dataset: AnalysisReadyPhosphoDataset,
        output: OutputTarget,
    ) -> dict[str, Path]:
        return publish_dataset(
            dataset,
            output.outdir,
            output_format=output.output_format,
        )

    def publish_kinase(
        self,
        result: KinaseWorkflowResult,
        output: OutputTarget,
    ) -> dict[str, Path]:
        return publish_kinase_workflow(
            result,
            output.outdir,
            output_format=output.output_format,
        )

    def publish_signalome(
        self,
        result: SignalomeWorkflowResult,
        output: OutputTarget,
    ) -> dict[str, Path]:
        return publish_signalome_workflow(
            result,
            output.outdir,
            output_format=output.output_format,
        )
