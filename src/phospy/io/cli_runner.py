"""CLI application service for typed command execution."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.requests import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.results import KinaseWorkflowResult, SignalomeWorkflowResult
from phospy.api.workflows import KinaseWorkflow, SignalomeWorkflow
from phospy.io.cli_commands import (
    CliCommand,
    DatasetBuildCommand,
    KinaseCommand,
    SignalomeCommand,
)
from phospy.io.cli_publishing import CliOutputPublisher, WorkflowOutputPublisher
from phospy.io.cli_request_factory import (
    build_kinase_workflow_request,
    build_signalome_workflow_request,
)
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset


class DatasetBuilderLike(Protocol):
    """Execution contract for dataset build collaborator."""

    def run(self, request: DatasetBuildRequest) -> AnalysisReadyPhosphoDataset: ...


class KinaseWorkflowLike(Protocol):
    """Execution contract for kinase workflow collaborator."""

    def run(self, request: KinaseWorkflowRequest) -> KinaseWorkflowResult: ...


class SignalomeWorkflowLike(Protocol):
    """Execution contract for signalome workflow collaborator."""

    def run(
        self,
        request: SignalomeWorkflowRequest,
    ) -> SignalomeWorkflowResult: ...


@dataclass(frozen=True, slots=True)
class CommandRunResult:
    """Published command outcome."""

    command_name: str
    written: dict[str, Path]


class CliCommandRunner:
    """Run typed CLI commands and publish their outputs."""

    def __init__(
        self,
        *,
        dataset_builder: DatasetBuilderLike | None = None,
        kinase_workflow: KinaseWorkflowLike | None = None,
        signalome_workflow: SignalomeWorkflowLike | None = None,
        publisher: CliOutputPublisher | None = None,
    ) -> None:
        self._dataset_builder = (
            AnalysisReadyDatasetBuilder()
            if dataset_builder is None
            else dataset_builder
        )
        self._kinase_workflow = (
            KinaseWorkflow() if kinase_workflow is None else kinase_workflow
        )
        self._signalome_workflow = (
            SignalomeWorkflow() if signalome_workflow is None else signalome_workflow
        )
        self._publisher = WorkflowOutputPublisher() if publisher is None else publisher

    def run(self, command: CliCommand) -> CommandRunResult:
        if isinstance(command, DatasetBuildCommand):
            dataset = self._dataset_builder.run(command.request)
            written = self._publisher.publish_dataset(dataset, command.output)
            return CommandRunResult(command_name="dataset-build", written=written)
        if isinstance(command, KinaseCommand):
            kinase_result = self._run_kinase(command)
            written = self._publisher.publish_kinase(kinase_result, command.output)
            return CommandRunResult(command_name="kinase", written=written)
        if isinstance(command, SignalomeCommand):
            kinase_result = self._run_kinase(command.kinase)
            request = build_signalome_workflow_request(command, kinase_result)
            signalome_result = self._signalome_workflow.run(request)
            written = self._publisher.publish_signalome(
                signalome_result, command.output
            )
            return CommandRunResult(command_name="signalome", written=written)
        raise RuntimeError(f"unsupported command type: {type(command)!r}")

    def _run_kinase(self, command: KinaseCommand) -> KinaseWorkflowResult:
        dataset = self._dataset_builder.run(command.dataset_request)
        request = build_kinase_workflow_request(command, dataset)
        return self._kinase_workflow.run(request)
