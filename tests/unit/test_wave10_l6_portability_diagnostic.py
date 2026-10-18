from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from tools.diagnostics.capture_wave10_l6_portability import (
    _align_prediction_matrices,
    _checkpoint_differences,
    _public_environment,
    _ranking_differences,
    _sequence_difference,
    array_payload,
    matrix_metrics,
    verify_manifest,
    write_manifest,
    write_redacted_manifest,
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


def test_public_environment_excludes_runner_and_private_paths() -> None:
    environment = {
        "source": {
            "commit": "a" * 40,
            "tree": "b" * 40,
            "parents": ["c" * 40],
            "branch": "diagnostics/wave10-l6-portability",
            "github": {
                "github_repository": "owner/repo",
                "github_run_id": "123",
                "runner_name": "private-runner-name",
            },
        },
        "operating_system": {
            "system": "Linux",
            "release": "6",
            "version": "version",
            "platform": "Linux-6-x86_64",
            "architecture": "x86_64",
            "processor": "processor identifier",
        },
        "python": {
            "version": "3.12.0",
            "implementation": "CPython",
            "full_version": "private build details",
            "executable": "/private/path/python",
        },
        "dependencies": {"numpy": "2.0.0"},
        "threadpools": [
            {
                "internal_api": "openblas",
                "user_api": "blas",
                "prefix": "libopenblas",
                "version": "0.3.0",
                "architecture": "Haswell",
                "threading_layer": "pthreads",
                "num_threads": 2,
                "filepath": "/private/path/libopenblas.so",
            }
        ],
        "numpy_runtime": {
            "cpu_features": {"AVX": True, "AVX512_SKX": False, "SSE": True},
            "configuration": {
                "Build Dependencies": {
                    "blas": {
                        "name": "accelerate",
                        "found": True,
                        "version": "unknown",
                        "detection method": "system",
                        "lib directory": "/private/path/lib",
                    }
                }
            },
        },
    }

    public = _public_environment(environment)
    serialized = json.dumps(public, sort_keys=True)

    assert "private" not in serialized
    assert "runner_name" not in serialized
    assert "filepath" not in serialized
    assert "lib directory" not in serialized
    assert "SSE" not in serialized
    assert public["numpy_cpu_features"] == {
        "AVX": True,
        "AVX2": False,
        "AVX512F": False,
        "AVX512CD": False,
        "AVX512_SKX": False,
        "FMA3": False,
    }
    assert public["numpy_blas_lapack"] == {
        "blas": {
            "name": "accelerate",
            "found": True,
            "version": "unknown",
            "detection method": "system",
            "openblas configuration": None,
        }
    }


def test_ranking_differences_reports_membership_and_order() -> None:
    left = {
        "candidate_tie_widths": ["0x0.0p+0"],
        "by_kinase": {
            "K": {
                "top_1": ["a"],
                "top_5": ["a", "b"],
                "top_10": ["a", "b"],
                "tie_groups": {},
            }
        },
    }
    right = {
        "candidate_tie_widths": ["0x0.0p+0"],
        "by_kinase": {
            "K": {
                "top_1": ["b"],
                "top_5": ["b", "a"],
                "top_10": ["b", "a"],
                "tie_groups": {},
            }
        },
    }

    result = _ranking_differences(left, right)

    assert result["top_order_or_membership_change_count"] == 3
    changes = cast(list[dict[str, object]], result["top_order_or_membership_changes"])
    assert changes[0]["membership_equal"] is False
    assert changes[1]["membership_equal"] is True


def test_sequence_difference_supports_structured_assignments() -> None:
    left = [{"kinase": "A", "score": 1.0}, {"kinase": "B", "score": 0.5}]
    right = [{"kinase": "B", "score": 0.5}, {"kinase": "A", "score": 1.0}]

    result = _sequence_difference(left, right)

    assert result["ordering_equal"] is False
    assert result["included_only_right"] == []
    assert result["excluded_from_right"] == []
    assert result["first_order_difference"] == {
        "position": 0,
        "left": left[0],
        "right": right[0],
    }


def test_redacted_manifest_uses_self_excluding_policy(tmp_path: Path) -> None:
    (tmp_path / "summary.md").write_text("summary\n", encoding="utf-8", newline="\n")

    write_redacted_manifest(tmp_path)
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))

    assert manifest["manifest_policy"].startswith("self-excluding")
    assert verify_manifest(tmp_path) == manifest["files"]
