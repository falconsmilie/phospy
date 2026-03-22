from __future__ import annotations

import argparse

from .pipeline import run_core_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the core phosrpy pipeline.")
    parser.add_argument("--total", required=True, help="Path to total proteome TSV.")
    parser.add_argument("--phospho", required=True, help="Path to phosphoproteome TSV.")
    parser.add_argument("--outdir", required=True, help="Directory for outputs.")
    parser.add_argument(
        "--pred-mat", help="Optional predMat CSV for kinase activity and KSEA."
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    run_core_pipeline(
        total_path=args.total,
        phospho_path=args.phospho,
        outdir=args.outdir,
        pred_mat_path=args.pred_mat,
    )


if __name__ == "__main__":
    main()
