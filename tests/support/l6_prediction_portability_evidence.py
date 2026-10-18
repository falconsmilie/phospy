from __future__ import annotations

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
) -> Path | None:
    """Attach L6 measurements to JUnit and optionally persist structured evidence."""
    lane = os.environ.get(L6_PREDICTION_LANE_ENV, "local")
    runtime = _runtime_identity()
    measurements = {
        "mean_abs_diff": parity.mean_abs_diff,
        "max_abs_diff": parity.max_abs_diff,
        "mean_pearson_corr": parity.mean_pearson_corr,
        "mean_spearman_corr": parity.mean_spearman_corr,
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
        "acceptance_bounds": bounds,
    }
    evidence_path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return evidence_path


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


__all__ = [
    "L6_PREDICTION_EVIDENCE_PATH_ENV",
    "L6_PREDICTION_LANE_ENV",
    "record_l6_prediction_portability_evidence",
]
