"""CLI parsing and command handlers."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from phospy.api.builders import AnalysisReadyDatasetBuilder
from phospy.api.configs import (
    KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES,
    KINASE_ACTIVITY_DEFAULT_THRESHOLD,
    KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES,
    KINASE_ADAPTIVE_POLICY_STABLE,
    KINASE_PREDICTION_DEFAULT_ITERATIONS,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeConfig,
)
from phospy.api.requests import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.workflows import KinaseWorkflow, SignalomeWorkflow
from phospy.errors import PhosPyError
from phospy.errors.input import PhosPyInputError
from phospy.io.publishers.workflows import (
    publish_dataset,
    publish_kinase_workflow,
    publish_signalome_workflow,
)
from phospy.references.models import Organism, ReferencePreset

CLI_EXIT_SUCCESS = 0
CLI_EXIT_INTERNAL_ERROR = 1
CLI_EXIT_USER_ERROR = 2


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI."""

    parser = build_parser()
    try:
        args = parser.parse_args(list(argv) if argv is not None else None)
    except SystemExit as exc:
        return int(exc.code)
    try:
        if args.command == "dataset-build":
            _run_dataset_build(args)
            return CLI_EXIT_SUCCESS
        if args.command == "kinase":
            _run_kinase(args)
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
    """Build the CLI parser for supported public commands."""

    parser = argparse.ArgumentParser(
        prog="phospy",
        description=(
            "PhosPy CLI. Supported commands: dataset-build, kinase, signalome."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    dataset_build = subparsers.add_parser(
        "dataset-build",
        help="Build an analysis-ready dataset from input files.",
    )
    _add_dataset_input_arguments(dataset_build)
    _add_output_arguments(dataset_build)

    kinase = subparsers.add_parser(
        "kinase",
        help="Run kinase workflow from input files.",
    )
    _add_dataset_input_arguments(kinase)
    _add_output_arguments(kinase)
    _add_kinase_runtime_arguments(kinase)

    signalome = subparsers.add_parser(
        "signalome",
        help="Run dataset -> kinase -> signalome workflow from input files.",
    )
    _add_dataset_input_arguments(signalome)
    _add_output_arguments(signalome)
    _add_kinase_runtime_arguments(signalome)
    signalome.add_argument(
        "--substrate-support-cutoff",
        type=float,
        default=0.5,
        help=("Prediction support cutoff for selecting kinase-supported substrates."),
    )
    signalome.add_argument(
        "--network-correlation-threshold",
        type=float,
        default=0.5,
        help="Absolute correlation threshold for kinase network edge inclusion.",
    )
    signalome.add_argument(
        "--assignment-policy",
        default=SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
        choices=["cutoff_binary", "weighted_top"],
        help=(
            "Signalome assignment policy: cutoff-based binary support or "
            "weighted top-kinase fractional support."
        ),
    )
    return parser


def _add_dataset_input_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--phospho",
        type=Path,
        required=True,
        help="Path to phospho matrix (.csv, .tsv, .txt as tab-separated, or .parquet).",
    )
    parser.add_argument(
        "--site-metadata",
        type=Path,
        required=True,
        help=(
            "Path to site metadata table "
            "(.csv, .tsv, .txt as tab-separated, or .parquet)."
        ),
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


def _add_kinase_runtime_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--reference",
        default="auto",
        choices=["auto", "human", "mouse", "rat"],
        help="Reference preset for workflow execution.",
    )
    parser.add_argument(
        "--scoring-min-substrates",
        type=int,
        default=KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
        help=(
            "Minimum quantified substrates per kinase for scoring "
            f"(must be >= {KINASE_SCORING_MIN_SUBSTRATES_FLOOR})."
        ),
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
        help=(
            "Mode-dependent ensemble_size: deterministic lane uses it as selected "
            "kinase cap; adaptive lane uses it as number of ensemble executions."
        ),
    )
    parser.add_argument(
        "--prediction-mode",
        default=KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
        choices=["deterministic_ranking", "adaptive_ensemble"],
        help="Prediction lane: deterministic ranking or adaptive ensemble science.",
    )
    parser.add_argument(
        "--prediction-adaptive-policy",
        default=KINASE_ADAPTIVE_POLICY_STABLE,
        choices=["stable", "r_parity"],
        help="Adaptive sampling policy when --prediction-mode=adaptive_ensemble.",
    )
    parser.add_argument(
        "--prediction-n-iterations",
        type=int,
        default=KINASE_PREDICTION_DEFAULT_ITERATIONS,
        help="Adaptive sampling iterations per ensemble when adaptive mode is used.",
    )
    parser.add_argument(
        "--prediction-random-state",
        type=int,
        default=None,
        help="Optional adaptive prediction random state.",
    )
    parser.add_argument(
        "--skip-activity",
        action="store_true",
        help="Disable activity-stage output.",
    )
    parser.add_argument(
        "--activity-threshold",
        type=float,
        default=KINASE_ACTIVITY_DEFAULT_THRESHOLD,
        help="Activity threshold when activity stage is enabled.",
    )
    parser.add_argument(
        "--activity-min-substrates",
        type=int,
        default=KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES,
        help="Minimum selected substrates per kinase for activity outputs.",
    )
    parser.add_argument(
        "--activity-top-n-substrates",
        type=int,
        default=KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES,
        help="Top-N predicted substrates per kinase used in weighted activity.",
    )


def _run_dataset_build(args: argparse.Namespace) -> None:
    dataset = _build_dataset_from_args(args)
    written = publish_dataset(dataset, args.outdir, output_format=args.output_format)
    _print_written_summary("dataset-build", written)


def _run_kinase(args: argparse.Namespace) -> None:
    result = _run_kinase_workflow_from_args(args)
    written = publish_kinase_workflow(
        result,
        args.outdir,
        output_format=args.output_format,
    )
    _print_written_summary("kinase", written)


def _run_signalome(args: argparse.Namespace) -> None:
    kinase_result = _run_kinase_workflow_from_args(args)
    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(
                substrate_support_cutoff=args.substrate_support_cutoff,
                network_correlation_threshold=args.network_correlation_threshold,
                assignment_policy=args.assignment_policy,
            ),
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
        organism=_organism_from_value(args.organism),
    )


def _run_kinase_workflow_from_args(args: argparse.Namespace):
    dataset = _build_dataset_from_args(args)
    reference = _reference_preset_from_value(args.reference)
    activity_config = (
        None
        if args.skip_activity
        else KinaseActivityConfig(
            enabled=True,
            threshold=args.activity_threshold,
            min_substrates=args.activity_min_substrates,
            top_n_substrates=args.activity_top_n_substrates,
        )
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=reference,
        scoring_config=KinaseScoringConfig(
            min_substrates=args.scoring_min_substrates,
        ),
        prediction_config=KinasePredictionConfig(
            top_k=args.prediction_top_k,
            ensemble_size=args.prediction_ensemble_size,
            mode=args.prediction_mode,
            adaptive_policy=args.prediction_adaptive_policy,
            n_iterations=args.prediction_n_iterations,
            random_state=args.prediction_random_state,
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


def _organism_from_value(value: Organism | str | None) -> Organism | None:
    if value is None:
        return None
    if isinstance(value, Organism):
        return value
    if not isinstance(value, str):
        raise PhosPyInputError(
            "unsupported organism value type. expected Organism, str, or None"
        )
    normalized = value.strip().lower()
    for organism in Organism:
        if organism.value == normalized:
            return organism
    supported = ", ".join(member.value for member in Organism)
    raise PhosPyInputError(
        f"unsupported organism '{value}'. supported organisms: {supported}"
    )


def _reference_preset_from_value(value: ReferencePreset | str) -> ReferencePreset:
    if isinstance(value, ReferencePreset):
        return value
    if not isinstance(value, str):
        raise PhosPyInputError(
            "unsupported reference preset value type. expected ReferencePreset or str"
        )
    normalized = value.strip().lower()
    for preset in ReferencePreset:
        if preset.value == normalized:
            return preset
    supported = ", ".join(member.value for member in ReferencePreset)
    raise PhosPyInputError(
        f"unsupported reference preset '{value}'. supported presets: {supported}"
    )
