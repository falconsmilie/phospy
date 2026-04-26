#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
DEFAULT_OUTPUT_DIR = ROOT / "tests" / "fixtures" / "public_workflow_reference"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate rewrite-owned public predMat benchmark references "
            "(stable and r_parity adaptive-policy lanes)."
        )
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where predmat_rewrite_*.csv files are written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if str(SRC) not in sys.path:
        sys.path.insert(0, str(SRC))
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    from tests.support.public_predmat_parity_metrics import _run_public_predmat_lane

    stable = _run_public_predmat_lane(
        adaptive_policy="stable",
        reverse_reference_order=False,
    )
    r_parity = _run_public_predmat_lane(
        adaptive_policy="r_parity",
        reverse_reference_order=False,
    )

    stable.to_csv(output_dir / "predmat_rewrite_stable.csv")
    r_parity.to_csv(output_dir / "predmat_rewrite_r_parity.csv")


if __name__ == "__main__":
    main()
