from __future__ import annotations

import json
from pathlib import Path

from tests.support.l6_prediction_parity_metrics import TableParityMetrics
from tests.support.l6_prediction_portability_evidence import (
    L6_PREDICTION_EVIDENCE_PATH_ENV,
    L6_PREDICTION_LANE_ENV,
    record_l6_prediction_portability_evidence,
)

TEST_NODE_ID = (
    "tests/parity/test_l6_prediction_parity.py::"
    "test_l6_prediction_matrix_numeric_parity_matches_promoted_reference_surfaces"
)


def _parity_metrics() -> TableParityMetrics:
    return TableParityMetrics(
        observed_shape=(600, 27),
        expected_shape=(600, 28),
        shared_site_count=600,
        shared_kinase_count=27,
        mean_abs_diff=1.25e-5,
        max_abs_diff=8.5e-4,
        mean_pearson_corr=0.9999998,
        mean_spearman_corr=0.9999997,
    )


def test_l6_portability_evidence_records_junit_properties_without_output_path(
    monkeypatch,
) -> None:
    junit_properties: list[tuple[str, object]] = []
    monkeypatch.delenv(L6_PREDICTION_EVIDENCE_PATH_ENV, raising=False)
    monkeypatch.setenv(L6_PREDICTION_LANE_ENV, "local-unit")

    written = record_l6_prediction_portability_evidence(
        test_node_id=TEST_NODE_ID,
        junit_properties=junit_properties,
        parity=_parity_metrics(),
        mean_abs_diff_ceiling=5e-5,
        max_abs_diff_ceiling=3e-3,
        pearson_floor=0.999999,
        spearman_floor=0.999999,
    )

    assert written is None
    properties = dict(junit_properties)
    assert properties["l6_lane"] == "local-unit"
    assert properties["l6_prediction_matrix_mean_abs_diff"] == 1.25e-5
    assert properties["l6_prediction_matrix_max_abs_diff"] == 8.5e-4
    assert properties["l6_python"]
    assert properties["l6_numpy"]


def test_l6_portability_evidence_writes_structured_runtime_record(
    monkeypatch,
    tmp_path: Path,
) -> None:
    evidence_path = tmp_path / "reports" / "l6-portability.json"
    junit_properties: list[tuple[str, object]] = []
    monkeypatch.setenv(L6_PREDICTION_EVIDENCE_PATH_ENV, str(evidence_path))
    monkeypatch.setenv(L6_PREDICTION_LANE_ENV, "Windows-py3.12")
    monkeypatch.setenv("GITHUB_SHA", "a" * 40)

    written = record_l6_prediction_portability_evidence(
        test_node_id=TEST_NODE_ID,
        junit_properties=junit_properties,
        parity=_parity_metrics(),
        mean_abs_diff_ceiling=5e-5,
        max_abs_diff_ceiling=3e-3,
        pearson_floor=0.999999,
        spearman_floor=0.999999,
    )

    assert written == evidence_path
    record = json.loads(evidence_path.read_text(encoding="utf-8"))
    assert record["schema"] == "phospy-l6-prediction-portability-v1"
    assert record["test_node_id"] == TEST_NODE_ID
    assert record["lane"] == "Windows-py3.12"
    assert record["source"]["commit"] == "a" * 40
    assert record["surface"] == {
        "expected_shape": [600, 28],
        "name": "prediction_matrix",
        "observed_shape": [600, 27],
        "shared_kinase_count": 27,
        "shared_site_count": 600,
    }
    assert record["measurements"]["mean_abs_diff"] == 1.25e-5
    assert record["measurements"]["max_abs_diff"] == 8.5e-4
    assert record["acceptance_bounds"]["mean_abs_diff_ceiling"] == 5e-5
    assert record["acceptance_bounds"]["max_abs_diff_ceiling"] == 3e-3
    assert set(record["runtime"]["dependencies"]) == {
        "numpy",
        "pandas",
        "scikit-learn",
        "scipy",
    }
    assert record["runtime"]["python"]["version"]
    assert record["runtime"]["platform"]["platform"]
