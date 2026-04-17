from __future__ import annotations

import argparse

from .api.contracts import DatasetLoadOptions, KinaseActivityConfig
from .errors import PhospyValidationError
from .internal.defaults import (
    DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES,
    DEFAULT_KINASE_ACTIVITY_THRESHOLD,
    DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES,
    DEFAULT_LOCALIZATION_THRESHOLD,
    DEFAULT_MAX_UNMATCHED_FRACTION,
    DEFAULT_MIN_OBSERVED_VALUES,
    DEFAULT_PHOSPHO_SENTINEL,
    DEFAULT_TOTAL_SENTINEL,
)
from .internal.pipeline import PipelineRunner
from .preprocessing import CorePreprocessingConfig

CLI_EXIT_SUCCESS = 0
CLI_EXIT_USER_ERROR = 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the core phospy pipeline.")
    parser.add_argument("--total", required=True, help="Path to total proteome TSV.")
    parser.add_argument("--phospho", required=True, help="Path to phosphoproteome TSV.")
    parser.add_argument("--outdir", required=True, help="Directory for outputs.")
    parser.add_argument(
        "--pred-mat", help="Optional predMat CSV for kinase activity and KSEA."
    )
    parser.add_argument(
        "--phospho-encoding",
        help=(
            "Optional text encoding for the phosphoproteome table. "
            "Defaults to utf-8 when omitted."
        ),
    )
    parser.add_argument(
        "--localization-threshold",
        type=float,
        default=DEFAULT_LOCALIZATION_THRESHOLD,
        help="Minimum localisation probability to retain a phosphosite.",
    )
    parser.add_argument(
        "--min-observed",
        type=int,
        default=DEFAULT_MIN_OBSERVED_VALUES,
        help="Minimum number of observed values required per row.",
    )
    parser.add_argument(
        "--total-sentinel",
        type=float,
        default=DEFAULT_TOTAL_SENTINEL,
        help="Sentinel value to treat as missing in the total proteome table.",
    )
    parser.add_argument(
        "--phospho-sentinel",
        type=float,
        default=DEFAULT_PHOSPHO_SENTINEL,
        help="Sentinel value to treat as missing in the phosphoproteome table.",
    )
    parser.add_argument(
        "--kinase-activity-threshold",
        type=float,
        default=DEFAULT_KINASE_ACTIVITY_THRESHOLD,
        help="Score threshold used for downstream kinase activity summaries.",
    )
    parser.add_argument(
        "--kinase-activity-min-substrates",
        type=int,
        default=DEFAULT_KINASE_ACTIVITY_MIN_SUBSTRATES,
        help="Minimum substrate count used for downstream kinase activity summaries.",
    )
    parser.add_argument(
        "--kinase-activity-top-n-substrates",
        type=int,
        default=DEFAULT_KINASE_ACTIVITY_TOP_N_SUBSTRATES,
        help="Top-N substrates used for weighted downstream kinase activity summaries.",
    )
    parser.add_argument(
        "--max-unmatched-fraction",
        type=float,
        default=DEFAULT_MAX_UNMATCHED_FRACTION,
        help=(
            "Maximum allowed fraction of phosphosite rows that may be dropped "
            "during protein correction because no matching protein row exists."
        ),
    )
    return parser


def build_cli_configs(
    args: argparse.Namespace,
) -> tuple[DatasetLoadOptions, CorePreprocessingConfig, KinaseActivityConfig]:
    """Build typed pipeline config objects from parsed CLI arguments."""

    dataset_options = DatasetLoadOptions(
        phospho_encoding=args.phospho_encoding,
    )
    preprocessing_config = CorePreprocessingConfig(
        localization_threshold=args.localization_threshold,
        min_observed=args.min_observed,
        total_sentinel=args.total_sentinel,
        phospho_sentinel=args.phospho_sentinel,
        max_unmatched_fraction=args.max_unmatched_fraction,
    )
    activity_config = KinaseActivityConfig(
        threshold=args.kinase_activity_threshold,
        min_substrates=args.kinase_activity_min_substrates,
        top_n_substrates=args.kinase_activity_top_n_substrates,
    )
    return dataset_options, preprocessing_config, activity_config


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        dataset_options, preprocessing_config, activity_config = build_cli_configs(args)
        pipeline = PipelineRunner.from_files(
            total_path=args.total,
            phospho_path=args.phospho,
            pred_mat_path=args.pred_mat,
            dataset_options=dataset_options,
            preprocessing_config=preprocessing_config,
            activity_config=activity_config,
        )
        pipeline.run(outdir=args.outdir)
    except PhospyValidationError as error:
        parser.exit(
            status=CLI_EXIT_USER_ERROR,
            message=f"{parser.prog}: error: {error}\n",
        )

    return None


if __name__ == "__main__":
    main()
