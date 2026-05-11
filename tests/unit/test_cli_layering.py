from __future__ import annotations

from argparse import Namespace
from pathlib import Path

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.api.requests import DatasetBuildRequest, KinaseWorkflowRequest
from phospy.errors.input import PhosPyInputError
from phospy.io.cli import main
from phospy.io.cli_commands import (
    DatasetBuildCommand,
    KinaseCommand,
    OutputTarget,
    SignalomeCommand,
)
from phospy.io.cli_errors import CLI_EXIT_INTERNAL_ERROR, CLI_EXIT_USER_ERROR
from phospy.io.cli_parser import build_parser
from phospy.io.cli_publishing import WorkflowOutputPublisher
from phospy.io.cli_request_factory import build_command
from phospy.io.cli_runner import CliCommandRunner
from phospy.references.models import Organism, ReferencePreset
from tests.support.rewrite_fixture_data import load_rat_l6_phospho, site_metadata_for


def test_cli_parser_parses_signalome_defaults_without_running_workflows() -> None:
    args = build_parser().parse_args(
        [
            "signalome",
            "--phospho",
            "phospho.csv",
            "--site-metadata",
            "site_metadata.csv",
        ]
    )

    assert args.command == "signalome"
    assert args.score_preconditioning_policy == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
    )


def test_build_command_converts_parsed_args_to_typed_signalome_command() -> None:
    args = build_parser().parse_args(
        [
            "signalome",
            "--phospho",
            "phospho.csv",
            "--site-metadata",
            "site_metadata.csv",
            "--organism",
            "rat",
            "--reference",
            "auto",
            "--prediction-top-k",
            "7",
            "--skip-activity",
            "--network-correlation-threshold",
            "0.6",
        ]
    )

    command = build_command(args)
    assert isinstance(command, SignalomeCommand)
    assert isinstance(command.kinase, KinaseCommand)
    assert command.kinase.dataset_request.organism == Organism.RAT
    assert command.kinase.references == ReferencePreset.AUTO
    assert command.kinase.prediction_config.top_k == 7
    assert command.kinase.activity_config is None
    assert command.signalome_config.output.network_correlation_threshold == 0.6


def test_cli_command_runner_executes_signalome_with_fake_collaborators() -> None:
    dataset_request = DatasetBuildRequest(
        phospho=Path("phospho.csv"),
        site_metadata=Path("site_metadata.csv"),
    )
    output = OutputTarget(outdir=Path("out"), output_format="csv")
    command = SignalomeCommand(
        kinase=KinaseCommand(
            dataset_request=dataset_request,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(),
            prediction_config=KinasePredictionConfig(),
            activity_config=KinaseActivityConfig(),
            output=output,
        ),
        signalome_config=SignalomeConfig(),
        output=output,
    )

    dataset_obj = object()
    kinase_result = object()
    signalome_result = object()

    class FakeDatasetBuilder:
        def __init__(self) -> None:
            self.requests: list[DatasetBuildRequest] = []

        def run(self, request: DatasetBuildRequest) -> object:
            self.requests.append(request)
            return dataset_obj

    class FakeKinaseWorkflow:
        def __init__(self) -> None:
            self.requests: list[KinaseWorkflowRequest] = []

        def run(self, request: KinaseWorkflowRequest) -> object:
            self.requests.append(request)
            return kinase_result

    class FakeSignalomeWorkflow:
        def __init__(self) -> None:
            self.requests: list[object] = []

        def run(self, request: object) -> object:
            self.requests.append(request)
            return signalome_result

    class FakePublisher:
        def __init__(self) -> None:
            self.published: list[tuple[str, object]] = []

        def publish_dataset(
            self, dataset: object, output: OutputTarget
        ) -> dict[str, Path]:
            self.published.append(("dataset", dataset))
            return {"dataset.path": output.outdir / "dataset.csv"}

        def publish_kinase(
            self, result: object, output: OutputTarget
        ) -> dict[str, Path]:
            self.published.append(("kinase", result))
            return {"kinase.path": output.outdir / "kinase.csv"}

        def publish_signalome(
            self,
            result: object,
            output: OutputTarget,
        ) -> dict[str, Path]:
            self.published.append(("signalome", result))
            return {"signalome.path": output.outdir / "signalome.csv"}

    dataset_builder = FakeDatasetBuilder()
    kinase_workflow = FakeKinaseWorkflow()
    signalome_workflow = FakeSignalomeWorkflow()
    publisher = FakePublisher()
    runner = CliCommandRunner(
        dataset_builder=dataset_builder,
        kinase_workflow=kinase_workflow,
        signalome_workflow=signalome_workflow,
        publisher=publisher,
    )

    result = runner.run(command)

    assert result.command_name == "signalome"
    assert result.written == {"signalome.path": Path("out/signalome.csv")}
    assert dataset_builder.requests == [dataset_request]
    assert len(kinase_workflow.requests) == 1
    assert kinase_workflow.requests[0].dataset is dataset_obj
    assert signalome_workflow.requests[0].kinase_result is kinase_result
    assert publisher.published == [("signalome", signalome_result)]


def test_workflow_output_publisher_writes_dataset_outputs(tmp_path: Path) -> None:
    phospho = load_rat_l6_phospho().head(16).copy(deep=True)
    site_metadata = site_metadata_for(phospho)
    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
        )
    )

    written = WorkflowOutputPublisher().publish_dataset(
        dataset,
        OutputTarget(outdir=tmp_path / "out", output_format="csv"),
    )

    assert "dataset.phospho" in written
    assert "dataset.site_metadata" in written
    assert "dataset.manifest" in written
    assert written["dataset.phospho"].exists()
    assert written["dataset.site_metadata"].exists()
    assert written["dataset.manifest"].exists()


def test_cli_main_presents_known_user_errors(
    monkeypatch,
    capsys,
) -> None:
    class FakeParser:
        def parse_args(self, _argv: list[str] | None) -> Namespace:
            return Namespace(command="dataset-build")

    def _raise_known_error(_args: Namespace) -> DatasetBuildCommand:
        raise PhosPyInputError("bad input")

    monkeypatch.setattr("phospy.io.cli.build_parser", lambda: FakeParser())
    monkeypatch.setattr("phospy.io.cli.build_command", _raise_known_error)

    exit_code = main([])
    captured = capsys.readouterr()

    assert exit_code == CLI_EXIT_USER_ERROR
    assert "PhosPyInputError: bad input" in captured.err


def test_cli_main_presents_unknown_errors(
    monkeypatch,
    capsys,
) -> None:
    class FakeParser:
        def parse_args(self, _argv: list[str] | None) -> Namespace:
            return Namespace(command="dataset-build")

    def _raise_unknown_error(_args: Namespace) -> DatasetBuildCommand:
        raise RuntimeError("boom")

    monkeypatch.setattr("phospy.io.cli.build_parser", lambda: FakeParser())
    monkeypatch.setattr("phospy.io.cli.build_command", _raise_unknown_error)

    exit_code = main([])
    captured = capsys.readouterr()

    assert exit_code == CLI_EXIT_INTERNAL_ERROR
    assert "UnhandledError: boom" in captured.err


def test_cli_main_delegates_typed_commands_without_argparse_leak(
    monkeypatch,
    capsys,
) -> None:
    class FakeParser:
        def parse_args(self, _argv: list[str] | None) -> Namespace:
            return Namespace(command="dataset-build")

    command = DatasetBuildCommand(
        request=DatasetBuildRequest(
            phospho=Path("phospho.csv"),
            site_metadata=Path("site_metadata.csv"),
        ),
        output=OutputTarget(outdir=Path("out"), output_format="csv"),
    )

    class FakeRunner:
        def run(self, incoming_command: object):
            assert incoming_command is command
            assert not isinstance(incoming_command, Namespace)
            return Namespace(
                command_name="dataset-build",
                written={"dataset.phospho": Path("out/dataset/phospho.csv")},
            )

    monkeypatch.setattr("phospy.io.cli.build_parser", lambda: FakeParser())
    monkeypatch.setattr("phospy.io.cli.build_command", lambda _args: command)
    monkeypatch.setattr("phospy.io.cli.CliCommandRunner", lambda: FakeRunner())

    exit_code = main([])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "phospy dataset-build completed." in captured.out
