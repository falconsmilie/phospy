from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import platform
import sys
from pathlib import Path

from tests.support.l6_prediction_parity_metrics import TableParityMetrics

L6_PREDICTION_EVIDENCE_PATH_ENV = "PHOSPY_L6_PREDICTION_EVIDENCE_PATH"
L6_PREDICTION_LANE_ENV = "PHOSPY_L6_PREDICTION_LANE"


def record_l6_prediction_portability_evidence(
    *,
    test_node_id: str,
    junit_properties: list[tuple[str, object]],
    parity: TableParityMetrics,
    mean_abs_diff_ceiling: float,
    max_abs_diff_ceiling: float,
    pearson_floor: float,
    spearman_floor: float,
    decision_evidence: dict[str, object],
) -> Path | None:
    """Attach L6 measurements to JUnit and optionally persist structured evidence."""
    lane = os.environ.get(L6_PREDICTION_LANE_ENV, "local")
    runtime = _runtime_identity()
    measurements = {
        "mean_abs_diff": parity.mean_abs_diff,
        "max_abs_diff": parity.max_abs_diff,
        "per_kinase_mean_pearson_corr": parity.mean_pearson_corr,
        "per_kinase_mean_spearman_corr": parity.mean_spearman_corr,
        "flattened_global_pearson_corr": parity.flattened_pearson_corr,
        "flattened_global_spearman_corr": parity.flattened_spearman_corr,
    }
    bounds = {
        "mean_abs_diff_ceiling": mean_abs_diff_ceiling,
        "max_abs_diff_ceiling": max_abs_diff_ceiling,
        "mean_pearson_corr_floor": pearson_floor,
        "mean_spearman_corr_floor": spearman_floor,
    }

    observed_properties = {
        "l6_lane": lane,
        "l6_prediction_matrix_mean_abs_diff": parity.mean_abs_diff,
        "l6_prediction_matrix_max_abs_diff": parity.max_abs_diff,
        "l6_prediction_matrix_mean_pearson_corr": parity.mean_pearson_corr,
        "l6_prediction_matrix_mean_spearman_corr": parity.mean_spearman_corr,
        "l6_platform": runtime["platform"]["platform"],
        "l6_python": runtime["python"]["version"],
        "l6_numpy": runtime["dependencies"]["numpy"],
        "l6_pandas": runtime["dependencies"]["pandas"],
        "l6_scipy": runtime["dependencies"]["scipy"],
        "l6_scikit_learn": runtime["dependencies"]["scikit-learn"],
    }
    junit_properties.extend(observed_properties.items())

    evidence_path_value = os.environ.get(L6_PREDICTION_EVIDENCE_PATH_ENV)
    if evidence_path_value is None:
        return None

    evidence_path = Path(evidence_path_value)
    evidence_path.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "schema": "phospy-l6-prediction-portability-v1",
        "test_node_id": test_node_id,
        "lane": lane,
        "source": {
            "commit": os.environ.get("GITHUB_SHA"),
            "repository": os.environ.get("GITHUB_REPOSITORY"),
            "run_id": os.environ.get("GITHUB_RUN_ID"),
            "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        },
        "runtime": runtime,
        "surface": {
            "name": "prediction_matrix",
            "observed_shape": list(parity.observed_shape),
            "expected_shape": list(parity.expected_shape),
            "shared_site_count": parity.shared_site_count,
            "shared_kinase_count": parity.shared_kinase_count,
        },
        "measurements": measurements,
        "correlation_metric_definitions": {
            "per_kinase_mean": (
                "arithmetic mean of 27 independently calculated per-kinase "
                "site-vector correlations; this is the maintained L6 gate metric"
            ),
            "flattened_global": (
                "one correlation over all aligned finite matrix cells; this is "
                "the Wave 10 diagnostic comparison metric"
            ),
        },
        "acceptance_bounds": bounds,
        "decision_contracts": decision_evidence,
    }
    evidence_path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return evidence_path


def compare_l6_prediction_portability_evidence(
    *,
    left_path: Path,
    right_path: Path,
    output_path: Path,
) -> dict[str, object]:
    """Assert deterministic scientific evidence from two fresh processes."""

    left = json.loads(left_path.read_text(encoding="utf-8"))
    right = json.loads(right_path.read_text(encoding="utf-8"))
    compared_sections = (
        "surface",
        "measurements",
        "correlation_metric_definitions",
        "acceptance_bounds",
        "decision_contracts",
    )
    differences = [
        section for section in compared_sections if left[section] != right[section]
    ]
    comparison = {
        "schema": "phospy-l6-fresh-process-repeatability-v1",
        "left": str(left_path),
        "right": str(right_path),
        "compared_sections": list(compared_sections),
        "differing_sections": differences,
        "exact_repeatability": not differences,
        "measurements": left["measurements"],
        "decision_contracts": left["decision_contracts"],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(comparison, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    if differences:
        raise AssertionError(
            "L6 fresh-process evidence differs in sections: " + ", ".join(differences)
        )
    return comparison


def _runtime_identity() -> dict[str, dict[str, str]]:
    return {
        "python": {
            "version": platform.python_version(),
            "full_version": sys.version,
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "platform": platform.platform(),
        },
        "dependencies": {
            distribution: importlib.metadata.version(distribution)
            for distribution in (
                "numpy",
                "pandas",
                "scipy",
                "scikit-learn",
            )
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("--left", required=True, type=Path)
    compare_parser.add_argument("--right", required=True, type=Path)
    compare_parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "compare":
        comparison = compare_l6_prediction_portability_evidence(
            left_path=args.left,
            right_path=args.right,
            output_path=args.output,
        )
        print("l6_repeatability=" + json.dumps(comparison, sort_keys=True))


if __name__ == "__main__":
    main()


__all__ = [
    "L6_PREDICTION_EVIDENCE_PATH_ENV",
    "L6_PREDICTION_LANE_ENV",
    "compare_l6_prediction_portability_evidence",
    "record_l6_prediction_portability_evidence",
]
