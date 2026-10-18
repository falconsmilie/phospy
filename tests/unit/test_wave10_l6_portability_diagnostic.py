from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from tools.diagnostics.capture_wave10_l6_portability import (
    _align_prediction_matrices,
    _checkpoint_differences,
    array_payload,
    matrix_metrics,
    verify_manifest,
    write_manifest,
)


def test_array_payload_preserves_binary64_values_and_hashes() -> None:
    values = np.asarray([[0.1, -0.0], [np.pi, 2.0**-52]], dtype=np.float64)

    payload = array_payload(values)
    values_hex = cast(list[str], payload["values_hex"])
    shape = cast(list[int], payload["shape"])

    reconstructed = np.asarray(
        [float.fromhex(value) for value in values_hex], dtype="<f8"
    ).reshape(shape)
    assert reconstructed.tobytes() == values.astype("<f8").tobytes()
    assert len(cast(str, payload["binary64_bytes_sha256"])) == 64
    assert len(cast(str, payload["canonical_sha256"])) == 64


def test_matrix_metrics_aligns_ids_and_reports_largest_difference() -> None:
    expected = pd.DataFrame({"A": [0.1, 0.2], "B": [0.3, 0.4]}, index=["s1", "s2"])
    observed = pd.DataFrame({"B": [0.5, 0.3], "A": [0.2, 0.1]}, index=["s2", "s1"])

    expected_aligned, observed_aligned, difference = _align_prediction_matrices(
        observed, expected
    )
    metrics = matrix_metrics(observed, expected)

    assert expected_aligned.index.tolist() == ["s1", "s2"]
    assert observed_aligned.columns.tolist() == ["A", "B"]
    assert np.isclose(float(cast(Any, difference.at["s2", "B"])), 0.1)
    largest = cast(list[dict[str, object]], metrics["largest_100_absolute_differences"])
    assert largest[0]["row_id"] == "s2"
    assert largest[0]["kinase"] == "B"


def test_checkpoint_comparison_finds_first_divergence(tmp_path: Path) -> None:
    left = tmp_path / "left"
    right = tmp_path / "right"
    left.mkdir()
    right.mkdir()
    common = {"kind": "input", "kinase": "K", "ensemble": 0, "iteration": 0}
    (left / "checkpoints.json").write_text(
        json.dumps({"checkpoints": [{**common, "payload_sha256": "a"}]})
    )
    (right / "checkpoints.json").write_text(
        json.dumps({"checkpoints": [{**common, "payload_sha256": "b"}]})
    )

    result = _checkpoint_differences(left, right)

    assert result["difference_count"] == 1
    first = cast(dict[str, object], result["first_divergence"])
    assert first["position"] == 0


def test_manifest_is_self_excluding(tmp_path: Path) -> None:
    (tmp_path / "evidence.json").write_text("{}\n", encoding="utf-8", newline="\n")

    write_manifest(tmp_path)
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))

    assert manifest["manifest_policy"].startswith("self-excluding")
    assert [item["path"] for item in manifest["files"]] == ["evidence.json"]
    assert verify_manifest(tmp_path) == manifest["files"]


def test_manifest_verification_rejects_changed_file(tmp_path: Path) -> None:
    evidence = tmp_path / "evidence.json"
    evidence.write_text("{}\n", encoding="utf-8", newline="\n")
    write_manifest(tmp_path)
    evidence.write_text('{"changed":true}\n', encoding="utf-8", newline="\n")

    with np.testing.assert_raises(AssertionError):
        verify_manifest(tmp_path)
