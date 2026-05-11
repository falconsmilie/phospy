"""Factory functions for converting parsed CLI arguments into typed commands."""

from __future__ import annotations

import argparse

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
    SignalomeClusteringConfig,
    SignalomeConfig,
    SignalomeOutputConfig,
    SignalomePerformanceConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
)
from phospy.api.requests import (
    DatasetBuildRequest,
    KinaseWorkflowRequest,
    SignalomeWorkflowRequest,
)
from phospy.api.results import KinaseWorkflowResult
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.input import PhosPyInputError
from phospy.io.cli_commands import (
    CliCommand,
    DatasetBuildCommand,
    KinaseCommand,
    OutputTarget,
    SignalomeCommand,
)
from phospy.references.models import Organism, ReferencePreset


def build_command(args: argparse.Namespace) -> CliCommand:
    """Build a typed command from parsed CLI arguments."""

    command = args.command
    if command == "dataset-build":
        return DatasetBuildCommand(
            request=_dataset_build_request_from_args(args),
            output=_output_target_from_args(args),
        )
    if command == "kinase":
        return _kinase_command_from_args(args)
    if command == "signalome":
        kinase_command = _kinase_command_from_args(args)
        return SignalomeCommand(
            kinase=kinase_command,
            signalome_config=SignalomeConfig(
                scientific=SignalomeScientificConfig(
                    substrate_support_cutoff=args.substrate_support_cutoff,
                    assignment_policy=args.assignment_policy,  # type: ignore[arg-type]
                ),
                clustering=SignalomeClusteringConfig(
                    candidate_scoring_policy=(  # type: ignore[arg-type]
                        args.candidate_scoring_policy
                    ),
                    clustering_engine=args.clustering_engine,  # type: ignore[arg-type]
                ),
                validation=SignalomeValidationConfig(
                    score_preconditioning_policy=(  # type: ignore[arg-type]
                        args.score_preconditioning_policy
                    ),
                ),
                output=SignalomeOutputConfig(
                    network_correlation_threshold=args.network_correlation_threshold,
                    network_policy=args.network_policy,  # type: ignore[arg-type]
                ),
                performance=SignalomePerformanceConfig(
                    max_exact_tree_sites=args.max_exact_tree_sites,
                    max_full_candidate_scoring_sites=args.max_full_candidate_scoring_sites,
                ),
            ),
            output=_output_target_from_args(args),
        )
    raise RuntimeError(f"unknown command: {command}")


def build_kinase_workflow_request(
    command: KinaseCommand,
    dataset: AnalysisReadyPhosphoDataset,
) -> KinaseWorkflowRequest:
    """Create a workflow request from a typed kinase command and dataset."""

    return KinaseWorkflowRequest(
        dataset=dataset,
        references=command.references,
        scoring_config=command.scoring_config,
        prediction_config=command.prediction_config,
        activity_config=command.activity_config,
    )


def build_signalome_workflow_request(
    command: SignalomeCommand,
    kinase_result: KinaseWorkflowResult,
) -> SignalomeWorkflowRequest:
    """Create a workflow request from a typed signalome command and kinase result."""

    return SignalomeWorkflowRequest(
        kinase_result=kinase_result,
        config=command.signalome_config,
    )


def _dataset_build_request_from_args(args: argparse.Namespace) -> DatasetBuildRequest:
    return DatasetBuildRequest(
        phospho=args.phospho,
        site_metadata=args.site_metadata,
        sample_metadata=args.sample_metadata,
        total=args.total,
        organism=_organism_from_value(args.organism),
    )


def _output_target_from_args(args: argparse.Namespace) -> OutputTarget:
    return OutputTarget(
        outdir=args.outdir,
        output_format=args.output_format,
    )


def _kinase_command_from_args(args: argparse.Namespace) -> KinaseCommand:
    activity_config = (
        None
        if args.skip_activity
        else KinaseActivityConfig(
            enabled=True,
            method=args.activity_method,
            threshold=args.activity_threshold,
            min_substrates=args.activity_min_substrates,
            top_n_substrates=args.activity_top_n_substrates,
            ksea_min_substrates=args.activity_ksea_min_substrates,
            ksea_evidence_threshold=args.activity_ksea_evidence_threshold,
            ksea_p_value_method=args.activity_ksea_p_value_method,
            ksea_adjust_p_values=not args.activity_ksea_no_adjust_p_values,
        )
    )
    return KinaseCommand(
        dataset_request=_dataset_build_request_from_args(args),
        references=_reference_preset_from_value(args.reference),
        scoring_config=KinaseScoringConfig(min_substrates=args.scoring_min_substrates),
        prediction_config=KinasePredictionConfig(
            top_k=args.prediction_top_k,
            deterministic_max_selected_kinases=(
                args.prediction_deterministic_max_selected_kinases
            ),
            adaptive_ensemble_runs=args.prediction_adaptive_ensemble_runs,
            mode=args.prediction_mode,
            adaptive_policy=args.prediction_adaptive_policy,
            n_iterations=args.prediction_n_iterations,
            random_state=args.prediction_random_state,
        ),
        activity_config=activity_config,
        output=_output_target_from_args(args),
    )


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
