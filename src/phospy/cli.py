from __future__ import annotations

import argparse

from .pipeline import PhosRPipeline
from .validation.errors import RequestValidationError
from .validation.requests import CorePipelineRequest


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
        default=0.75,
        help="Minimum localisation probability to retain a phosphosite.",
    )
    parser.add_argument(
        "--min-observed",
        type=int,
        default=4,
        help="Minimum number of observed values required per row.",
    )
    parser.add_argument(
        "--total-sentinel",
        type=float,
        default=10.0,
        help="Sentinel value to treat as missing in the total proteome table.",
    )
    parser.add_argument(
        "--phospho-sentinel",
        type=float,
        default=12.0,
        help="Sentinel value to treat as missing in the phosphoproteome table.",
    )
    parser.add_argument(
        "--max-unmatched-fraction",
        type=float,
        default=0.0,
        help=(
            "Maximum allowed fraction of phosphosite rows that may be dropped "
            "during protein correction because no matching protein row exists."
        ),
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        request = CorePipelineRequest.validate_request(
            total_path=args.total,
            phospho_path=args.phospho,
            pred_mat_path=args.pred_mat,
            phospho_encoding=args.phospho_encoding,
            localization_threshold=args.localization_threshold,
            min_observed=args.min_observed,
            total_sentinel=args.total_sentinel,
            phospho_sentinel=args.phospho_sentinel,
            max_unmatched_fraction=args.max_unmatched_fraction,
        )
    except RequestValidationError as error:
        parser.error(str(error))

    pipeline = PhosRPipeline.from_request(request)
    pipeline.run(outdir=args.outdir)


if __name__ == "__main__":
    main()
