from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
import tomllib
from copy import deepcopy
from pathlib import Path
from typing import Any, cast

import pytest

from tools.architecture.capture_control_site_regression import (
    _invalid_quantitative_state_payload,
)

pytestmark = [pytest.mark.release_gate, pytest.mark.reproducibility]

EXPECTED_NORMALIZED_CASES_SHA256 = (
    "c6d6bc1e87190bfd972dfed0ce74ca0b5ca595b3993b9cba5ccc8db6cb7aa037"
)
EXPECTED_NORMALIZED_PAYLOAD_SHA256 = (
    "4615d7fc68bc84c030fb52edd7e021107386dc50e15899618f446418eeeab386"
)
EXPECTED_LEGACY_CONTROL_DISCOVERY_IDENTITY = "sha256-stable-json-v1:1d4091e2ba3aaae65d54b25b3688cf53a28a9a96fcf553de5f8805f4d69004f3"
PYTHON_VERSION_PROVENANCE_PATH = (
    "cases",
    "correction_from_restored_controls",
    "provenance",
    "python_version",
)
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CAPTURE_TOOL = (
    REPOSITORY_ROOT / "tools" / "architecture" / "capture_control_site_regression.py"
)


def _canonical_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _fresh_process_comparison_payload() -> dict[str, object]:
    """Capture with source-tree distribution metadata isolated from site packages."""

    environment = os.environ.copy()
    project = cast(
        dict[str, object],
        tomllib.loads((REPOSITORY_ROOT / "pyproject.toml").read_text("utf-8"))[
            "project"
        ],
    )
    assert project["name"] == "phospy"
    source_version = str(project["version"])
    with tempfile.TemporaryDirectory(prefix="phospy-control-site-regression-") as root:
        metadata_root = Path(root)
        dist_info = metadata_root / f"phospy-{source_version}.dist-info"
        dist_info.mkdir()
        (dist_info / "METADATA").write_bytes(
            (
                f"Metadata-Version: 2.1\nName: phospy\nVersion: {source_version}\n"
            ).encode()
        )
        environment["PYTHONPATH"] = os.pathsep.join(
            (
                str(metadata_root),
                str(REPOSITORY_ROOT / "src"),
                str(REPOSITORY_ROOT),
            )
        )
        completed = subprocess.run(
            (sys.executable, str(CAPTURE_TOOL)),
            cwd=REPOSITORY_ROOT,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )
    assert not completed.stderr
    payload = json.loads(completed.stdout)
    assert isinstance(payload, dict)
    typed_payload = cast(dict[str, object], payload)
    assert (
        _payload_value(
            typed_payload,
            (
                "cases",
                "correction_from_restored_controls",
                "provenance",
                "phospy_version",
            ),
        )
        == source_version
    )
    return typed_payload


def _payload_value(payload: dict[str, object], path: tuple[str, ...]) -> object:
    value: object = payload
    for key in path:
        value = cast(dict[str, object], value)[key]
    return value


def _set_payload_value(
    payload: dict[str, object], path: tuple[str, ...], value: object
) -> None:
    parent: dict[str, object] = payload
    for key in path[:-1]:
        parent = cast(dict[str, object], parent[key])
    parent[path[-1]] = value


def _assert_regression_payload_match(
    payload: dict[str, object], *, actual_python_version: str
) -> None:
    """Check every leaf while normalizing the Python runtime version."""

    actual_provenance = _payload_value(payload, PYTHON_VERSION_PROVENANCE_PATH)
    assert actual_provenance == actual_python_version
    actual_parts = actual_python_version.split(".")
    assert len(actual_parts) == 3
    assert tuple(int(part) for part in actual_parts[:2]) == sys.version_info[:2]
    cases = cast(dict[str, object], payload["cases"])
    assert payload["cases_sha256"] == _canonical_sha256(cases)

    qualified = cast(dict[str, object], deepcopy(payload))
    _set_payload_value(
        qualified,
        PYTHON_VERSION_PROVENANCE_PATH,
        "runtime-python",
    )
    qualified_cases = cast(dict[str, object], qualified["cases"])
    qualified["cases_sha256"] = _canonical_sha256(qualified_cases)
    assert qualified["cases_sha256"] == EXPECTED_NORMALIZED_CASES_SHA256
    assert _canonical_sha256(qualified) == EXPECTED_NORMALIZED_PAYLOAD_SHA256


def test_control_site_payload_matches_regression_baseline() -> None:
    payload = _fresh_process_comparison_payload()
    _assert_regression_payload_match(
        payload,
        actual_python_version=platform.python_version(),
    )
    cases = cast(dict[str, Any], payload["cases"])
    legacy_controls = cast(dict[str, Any], cases["legacy_controls"])
    legacy_controls_restored = cast(dict[str, Any], cases["legacy_controls_restored"])
    legacy_result = cast(dict[str, Any], cases["legacy_restored"])
    source_metadata = cast(dict[str, Any], legacy_controls["source_metadata"])

    assert legacy_controls_restored == legacy_controls
    assert source_metadata["source_version"] == "1.0.0"
    assert (
        source_metadata["sps_discovery_identity"] == legacy_result["discovery_identity"]
    )
    assert (
        source_metadata["sps_discovery_identity"]
        == EXPECTED_LEGACY_CONTROL_DISCOVERY_IDENTITY
    )


def test_invalid_quantitative_state_is_rejected_before_ranking_without_mutation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_if_ranked(*args: object, **kwargs: object) -> None:
        raise AssertionError("ranking executed for invalid quantitative state")

    monkeypatch.setattr(
        "phospy.preprocessing.control_sites.discovery_execution."
        "SpsDiscoveryExecutor.run",
        fail_if_ranked,
    )
    payload = _invalid_quantitative_state_payload()
    validation = cast(dict[str, Any], payload["validation_error"])
    issues = cast(list[dict[str, Any]], validation["issues"])

    assert [issue["code"] for issue in issues] == [
        "missing_scale_establishment_evidence",
        "missing_baseline_establishment_evidence",
    ]
    assert payload["matrix_after"] == payload["matrix_before"]
    assert payload["quantitative_state_after"] == (payload["quantitative_state_before"])
    assert payload["quantitative_state_identity_preserved"] is True
    assert payload["rejected_before_ranking"] is True


def test_control_site_regression_permits_only_python_patch_difference() -> None:
    payload = _fresh_process_comparison_payload()
    current = platform.python_version()
    major, minor, patch = (int(part) for part in current.split("."))
    simulated_runtime = f"{major}.{minor}.{patch + 1}"
    _set_payload_value(payload, PYTHON_VERSION_PROVENANCE_PATH, simulated_runtime)
    payload["cases_sha256"] = _canonical_sha256(payload["cases"])

    _assert_regression_payload_match(
        payload,
        actual_python_version=simulated_runtime,
    )


@pytest.mark.parametrize("change", ["scientific", "dependency", "major_minor"])
def test_control_site_regression_rejects_unapproved_payload_changes(
    change: str,
) -> None:
    payload = _fresh_process_comparison_payload()
    actual_python_version = platform.python_version()
    if change == "scientific":
        cases = cast(dict[str, object], payload["cases"])
        warnings = cast(list[str], cases["suspicious_warnings"])
        warnings.append("unapproved scientific change")
    elif change == "dependency":
        provenance = cast(
            dict[str, object],
            cast(dict[str, object], payload["cases"])[
                "correction_from_restored_controls"
            ],
        )["provenance"]
        dependency_versions = cast(
            dict[str, object],
            cast(dict[str, object], provenance)["dependency_versions"],
        )
        dependency_versions["numpy"] = "999.0.0"
    else:
        major, minor, patch = (int(part) for part in actual_python_version.split("."))
        actual_python_version = f"{major}.{minor + 1}.{patch}"
        _set_payload_value(
            payload,
            PYTHON_VERSION_PROVENANCE_PATH,
            actual_python_version,
        )
    payload["cases_sha256"] = _canonical_sha256(payload["cases"])

    with pytest.raises(AssertionError):
        _assert_regression_payload_match(
            payload,
            actual_python_version=actual_python_version,
        )
