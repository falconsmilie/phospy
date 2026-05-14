"""CLI parser construction for supported public commands."""

from __future__ import annotations

import argparse
from pathlib import Path

from phospy.api.configs import (
    KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES,
    KINASE_ACTIVITY_DEFAULT_THRESHOLD,
    KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES,
    KINASE_ACTIVITY_KSEA_DEFAULT_MIN_SUBSTRATES,
    KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION,
    KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
    KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
    KINASE_ADAPTIVE_POLICY_STABLE,
    KINASE_PREDICTION_DEFAULT_ADAPTIVE_ENSEMBLE_RUNS,
    KINASE_PREDICTION_DEFAULT_DETERMINISTIC_MAX_SELECTED_KINASES,
    KINASE_PREDICTION_DEFAULT_ITERATIONS,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
    SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
    SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
    SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
    SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD,
    SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY,
    SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT,
    SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
)


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for supported public commands."""

    parser = argparse.ArgumentParser(
        prog="phospy",
        description=(
            "PhosPy CLI. Supported commands: dataset-build, kinase, signalome. "
            "This CLI intentionally covers the file-first supported lane; use "
            "phospy.api for advanced preprocessing and full request/config control."
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
        help="Correlation threshold used by --network-policy for edge inclusion.",
    )
    signalome.add_argument(
        "--network-policy",
        default=SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
        choices=[
            SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY,
            SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD,
            SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
        ],
        help=(
            "Kinase network policy: positive-only thresholding, absolute-threshold "
            "unsigned edges, or signed absolute-threshold edges."
        ),
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
    signalome.add_argument(
        "--score-preconditioning-policy",
        default=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
        choices=[
            SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
            SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
        ],
        help=(
            "Downstream score preconditioning policy: fail when any all-missing "
            "row would be dropped (default), or opt into allow/report."
        ),
    )
    signalome.add_argument(
        "--clustering-engine",
        default=SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
        choices=[
            SIGNALOME_CLUSTERING_ENGINE_EXACT_PYTHON,
            SIGNALOME_CLUSTERING_ENGINE_SCIPY_HIERARCHICAL,
        ],
        help=(
            "Signalome clustering backend implementation. "
            "Use scipy_hierarchical for production defaults; "
            "exact_python for reference/debug parity checks."
        ),
    )
    signalome.add_argument(
        "--candidate-scoring-policy",
        default=SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
        choices=[
            SIGNALOME_CANDIDATE_SCORING_POLICY_FULL,
            SIGNALOME_CANDIDATE_SCORING_POLICY_SAMPLED,
        ],
        help=(
            "Signalome module-selection candidate scoring policy: full "
            "correlations or sampled within-cluster estimates."
        ),
    )
    signalome.add_argument(
        "--max-exact-tree-sites",
        type=int,
        default=SIGNALOME_MAX_EXACT_TREE_SITES_DEFAULT,
        help=(
            "Hard execution guard for exact tree construction. Signalome "
            "fails when interpreted site count exceeds this limit."
        ),
    )
    signalome.add_argument(
        "--max-full-candidate-scoring-sites",
        type=int,
        default=SIGNALOME_MAX_FULL_CANDIDATE_SCORING_SITES_DEFAULT,
        help=(
            "Hard execution guard for candidate_scoring_policy=full. Full "
            "candidate scoring fails above this site count."
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
    parser.add_argument(
        "--input-intensity-scale",
        default=None,
        choices=["linear", "log2"],
        help=(
            "Explicit declared input intensity scale. Required for pass-through "
            "identity preprocessing lanes to establish scientific scale state."
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
        "--prediction-deterministic-max-selected-kinases",
        type=int,
        default=KINASE_PREDICTION_DEFAULT_DETERMINISTIC_MAX_SELECTED_KINASES,
        help=("Maximum number of kinases retained in deterministic prediction mode."),
    )
    parser.add_argument(
        "--prediction-adaptive-ensemble-runs",
        type=int,
        default=KINASE_PREDICTION_DEFAULT_ADAPTIVE_ENSEMBLE_RUNS,
        help="Number of ensemble runs in adaptive prediction mode.",
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
        help=(
            "Adaptive prediction random state. Required when "
            "--prediction-mode=adaptive_ensemble."
        ),
    )
    parser.add_argument(
        "--skip-activity",
        action="store_true",
        help="Disable activity-stage output.",
    )
    parser.add_argument(
        "--activity-method",
        default=KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
        choices=[
            KINASE_ACTIVITY_METHOD_SIMPLIFIED_WEIGHTED_SUBSTRATE_ACTIVITY,
            KINASE_ACTIVITY_METHOD_KSEA_ZSCORE,
        ],
        help="Activity method: weighted heuristic or KSEA-style z-score.",
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
    parser.add_argument(
        "--activity-ksea-min-substrates",
        type=int,
        default=KINASE_ACTIVITY_KSEA_DEFAULT_MIN_SUBSTRATES,
        help="Minimum substrate count per kinase/condition for KSEA scoring.",
    )
    parser.add_argument(
        "--activity-ksea-evidence-threshold",
        type=float,
        default=None,
        help=(
            "KSEA evidence threshold for kinase-substrate membership. "
            "Defaults to --activity-threshold when omitted."
        ),
    )
    parser.add_argument(
        "--activity-ksea-p-value-method",
        default=KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION,
        choices=[KINASE_ACTIVITY_KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION],
        help="KSEA p-value method.",
    )
    parser.add_argument(
        "--activity-ksea-no-adjust-p-values",
        action="store_true",
        help="Disable Benjamini-Hochberg q-value adjustment for KSEA.",
    )
