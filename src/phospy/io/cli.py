"""Rewrite-lane CLI parsing and command handlers."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.api.requests import (
    DatasetBuildRequest,
    SignalomeWorkflowRequest,
    SimpleKinaseWorkflowRequest,
)
from phospy.api.workflows import KinaseWorkflow, SignalomeWorkflow
from phospy.errors import PhosPyError
from phospy.io.adapters import (
    organism_from_value,
    reference_preset_from_value,
)
from phospy.io.publishing import (
    publish_dataset,
    publish_signalome_workflow,
    publish_simple_kinase_workflow,
)
from phospy.transformations.models import MatrixTransformationState, TransformationState

CLI_EXIT_SUCCESS = 0
CLI_EXIT_INTERNAL_ERROR = 1
CLI_EXIT_USER_ERROR = 2


def main(argv: Sequence[str] | None = None) -> int:
    """Run the rewrite-lane CLI."""

    parser = build_parser()
    try:
        args = parser.parse_args(list(argv) if argv is not None else None)
    except SystemExit as exc:
        return int(exc.code)
    try:
        if args.command == "dataset-build":
            _run_dataset_build(args)
            return CLI_EXIT_SUCCESS
        if args.command == "simple-kinase":
            _run_simple_kinase(args)
            return CLI_EXIT_SUCCESS
        if args.command == "signalome":
            _run_signalome(args)
            return CLI_EXIT_SUCCESS
        raise RuntimeError(f"unknown command: {args.command}")
    except PhosPyError as exc:
        print(f"{exc.__class__.__name__}: {exc}", file=sys.stderr)
        return CLI_EXIT_USER_ERROR
    except Exception as exc:  # pragma: no cover - defensive CLI boundary
        print(f"UnhandledError: {exc}", file=sys.stderr)
        return CLI_EXIT_INTERNAL_ERROR


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for the supported rewrite lane."""

    parser = argparse.ArgumentParser(
        prog="phospy",
        description=(
            "PhosPy rewrite CLI. Supported commands: dataset-build, simple-kinase, signalome."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    dataset_build = subparsers.add_parser(
        "dataset-build",
        help="Build an analysis-ready dataset from input files.",
    )
    _add_dataset_input_arguments(dataset_build)
    _add_output_arguments(dataset_build)

    simple_kinase = subparsers.add_parser(
        "simple-kinase",
        help="Run simple kinase workflow from input files.",
    )
    _add_dataset_input_arguments(simple_kinase)
    _add_output_arguments(simple_kinase)
    _add_simple_kinase_runtime_arguments(simple_kinase)

    signalome = subparsers.add_parser(
        "signalome",
        help="Run dataset -> simple kinase -> signalome workflow from input files.",
    )
    _add_dataset_input_arguments(signalome)
    _add_output_arguments(signalome)
    _add_simple_kinase_runtime_arguments(signalome)
    signalome.add_argument(
        "--signalome-cutoff",
        type=float,
        default=0.5,
        help="Signalome cutoff threshold.",
    )
    return parser


def _add_dataset_input_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--phospho",
        type=Path,
        required=True,
        help="Path to phospho matrix (.csv, .tsv, or .parquet).",
    )
    parser.add_argument(
        "--site-metadata",
        type=Path,
        required=True,
        help="Path to site metadata table (.csv, .tsv, or .parquet).",
    )
    parser.add_argument(
        "--sample-metadata",
        type=Path,
        default=None,
        help="Optional sample metadata table.",
    )
    parser.add_argument(
        "--total",
        type=Path,
        default=None,
        help="Optional total proteome table.",
    )
    parser.add_argument(
        "--organism",
        default=None,
        choices=["human", "mouse", "rat"],
        help="Optional dataset organism.",
    )
    parser.add_argument(
        "--transformation-state",
        default=None,
        choices=["linear", "log2"],
        help=(
            "Explicit transformation state for quantitative inputs. "
            "Required unless state is established by a supported transformer."
        ),
    )


def _add_output_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("phospy-output"),
        help="Root output directory.",
    )
    parser.add_argument(
        "--output-format",
        default="csv",
        choices=["csv", "tsv", "parquet"],
        help="File format for written outputs.",
    )


def _add_simple_kinase_runtime_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--reference",
        default="auto",
        choices=["auto", "human", "mouse", "rat"],
        help="Reference preset for workflow execution.",
    )
    parser.add_argument(
        "--scoring-min-substrates",
        type=int,
        default=1,
        help="Minimum substrates per kinase for scoring.",
    )
    parser.add_argument(
        "--prediction-top-k",
        type=int,
        default=30,
        help="Top-k predicted substrate sites per kinase.",
    )
    parser.add_argument(
        "--prediction-ensemble-size",
        type=int,
        default=10,
        help="Number of kinases included in prediction matrix output.",
    )
    parser.add_argument(
        "--skip-activity",
        action="store_true",
        help="Disable activity-stage output.",
    )
    parser.add_argument(
        "--activity-threshold",
        type=float,
        default=0.6,
        help="Activity threshold when activity stage is enabled.",
    )


def _run_dataset_build(args: argparse.Namespace) -> None:
    dataset = _build_dataset_from_args(args)
    written = publish_dataset(dataset, args.outdir, output_format=args.output_format)
    _print_written_summary("dataset-build", written)


def _run_simple_kinase(args: argparse.Namespace) -> None:
    result = _run_simple_kinase_workflow_from_args(args)
    written = publish_simple_kinase_workflow(
        result,
        args.outdir,
        output_format=args.output_format,
    )
    _print_written_summary("simple-kinase", written)


def _run_signalome(args: argparse.Namespace) -> None:
    kinase_result = _run_simple_kinase_workflow_from_args(args)
    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(signalome_cutoff=args.signalome_cutoff),
        )
    )
    written = publish_signalome_workflow(
        signalome_result,
        args.outdir,
        output_format=args.output_format,
    )
    _print_written_summary("signalome", written)


def _dataset_build_request_from_args(args: argparse.Namespace) -> DatasetBuildRequest:
    return DatasetBuildRequest(
        phospho=args.phospho,
        site_metadata=args.site_metadata,
        sample_metadata=args.sample_metadata,
        total=args.total,
        organism=organism_from_value(args.organism),
        transformation_state=_transformation_state_from_value(
            args.transformation_state,
            has_total_matrix=args.total is not None,
        ),
    )


def _transformation_state_from_value(
    value: str | None,
    *,
    has_total_matrix: bool,
) -> TransformationState | None:
    if value is None:
        return None

    phospho_state: MatrixTransformationState
    total_state: MatrixTransformationState | None
    if value == "linear":
        phospho_state = MatrixTransformationState.linear(
            established_by="phospy.io.cli.user_declared"
        )
        total_state = (
            MatrixTransformationState.linear(
                established_by="phospy.io.cli.user_declared"
            )
            if has_total_matrix
            else None
        )
    else:
        phospho_state = MatrixTransformationState.log2(
            established_by="phospy.io.cli.user_declared"
        )
        total_state = (
            MatrixTransformationState.log2(established_by="phospy.io.cli.user_declared")
            if has_total_matrix
            else None
        )
    return TransformationState(
        phospho=phospho_state,
        total=total_state,
    )


def _run_simple_kinase_workflow_from_args(args: argparse.Namespace):
    dataset = _build_dataset_from_args(args)
    reference = reference_preset_from_value(args.reference)
    activity_config = (
        None
        if args.skip_activity
        else KinaseActivityConfig(enabled=True, threshold=args.activity_threshold)
    )
    request = SimpleKinaseWorkflowRequest(
        dataset=dataset,
        references=reference,
        scoring_config=KinaseScoringConfig(
            min_substrates=args.scoring_min_substrates,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=args.prediction_top_k,
            ensemble_size=args.prediction_ensemble_size,
        ),
        activity_config=activity_config,
    )
    return KinaseWorkflow().run(request)


def _build_dataset_from_args(args: argparse.Namespace):
    return AnalysisReadyDatasetBuilder().run(_dataset_build_request_from_args(args))


def _print_written_summary(command: str, written: dict[str, Path]) -> None:
    print(f"phospy {command} completed.")
    for key in sorted(written):
        print(f"{key}: {written[key]}")
