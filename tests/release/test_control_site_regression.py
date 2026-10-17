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
    "75d4d0d112fc7597e1b61f9b6d4d8018d3a61dcbbfd031215fc3ca5266ce2084"
)
EXPECTED_NORMALIZED_PAYLOAD_SHA256 = (
    "a11ce6b4dd26739ea2e0aef4dac94758b8ff92aa7b409bba2b3d5af70811ea23"
)
EXPECTED_LEGACY_CONTROL_DISCOVERY_IDENTITY = "sha256-stable-json-v1:1d4091e2ba3aaae65d54b25b3688cf53a28a9a96fcf553de5f8805f4d69004f3"
EXPECTED_CORRECTED_MATRIX_VALUES = (
    (2.0, 3.0, 2.0, 3.0),
    (1.5, 3.5, 1.5, 3.5),
    (5.625, 6.125, 5.375, 6.375),
    (8.0, 7.0, 9.0, 6.0),
)
EXPECTED_CORRECTION_DIAGNOSTICS = {
    ("batch_associated_variance", "after", "max_r_squared"): 0.0,
    ("batch_associated_variance", "after", "mean_r_squared"): 0.0,
    ("batch_associated_variance", "after", "median_r_squared"): 0.0,
    ("batch_associated_variance", "before", "max_r_squared"): 1.0,
    ("batch_associated_variance", "before", "mean_r_squared"): 0.740384615384615,
    ("batch_associated_variance", "before", "median_r_squared"): 0.9807692307692302,
    ("batch_associated_variance", "delta_mean_r_squared"): -0.740384615384615,
    ("batch_associated_variance", "delta_median_r_squared"): -0.9807692307692302,
    ("max_abs_adjustment",): 1.0,
    ("mean_abs_adjustment",): 0.53125,
}
EXPECTED_COMPLETE_RANKING_NUMERICS = (
    (0.375, 0.375, 0.328125, 0.015625),
    (0.5, 0.75, 0.5, 0.75, 0.875, 0.375, 0.125, 0.125),
    (0.2, 0.1, 0.2, 0.1, 0.1, 0.2, 2.0, 3.0),
)
EXPECTED_PARTIAL_RANKING_NUMERICS = (
    (0.6759846280642754, 0.46060675273531304, 0.018100631647363866, 0.375),
    (
        0.5,
        0.75,
        0.8333333333333334,
        0.875,
        0.375,
        0.5,
        0.125,
        0.125,
        0.16666666666666666,
        0.5,
        0.75,
    ),
    (0.2, 0.1, 0.2, 0.1, 0.2, 0.3, 2.0, 3.0, 4.0, 0.2, 0.1),
)
RANKING_NUMERIC_BASELINES = {
    "complete": EXPECTED_COMPLETE_RANKING_NUMERICS,
    "partial": EXPECTED_PARTIAL_RANKING_NUMERICS,
    "permuted": EXPECTED_PARTIAL_RANKING_NUMERICS,
    "current_restored": EXPECTED_PARTIAL_RANKING_NUMERICS,
    "legacy_input": (
        EXPECTED_COMPLETE_RANKING_NUMERICS[0],
        (),
        EXPECTED_COMPLETE_RANKING_NUMERICS[2],
    ),
    "legacy_restored": EXPECTED_COMPLETE_RANKING_NUMERICS,
}
CORRECTION_CASE_PATH = ("cases", "correction_from_restored_controls")
PLATFORM_SENSITIVE_NUMBER = "platform-sensitive-number"
PLATFORM_SENSITIVE_HASH = "platform-sensitive-hash"
PLATFORM_SENSITIVE_IDENTITY = "platform-sensitive-identity"
PYTHON_VERSION_PROVENANCE_PATH = (
    "cases",
    "correction_from_restored_controls",
    "provenance",
    "python_version",
)
DEPENDENCY_VERSIONS_PATH = (
    "cases",
    "correction_from_restored_controls",
    "provenance",
    "dependency_versions",
)
CANONICAL_DEPENDENCY_VERSIONS = {
    "numpy": "2.2.6",
    "pandas": "2.2.3",
    "scipy": "1.15.2",
    "scikit-learn": "1.6.1",
}
APPROVED_DEPENDENCY_VERSION_PROFILES = (
    CANONICAL_DEPENDENCY_VERSIONS,
    {
        "numpy": "1.24.0",
        "pandas": "2.0.0",
        "scipy": "1.10.0",
        "scikit-learn": "1.6.0",
    },
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


def _normalize_platform_sensitive_payload(
    value: Any, path: tuple[str | int, ...] = ()
) -> Any:
    """Project backend-sensitive derived values to stable sentinels."""

    if isinstance(value, float):
        return PLATFORM_SENSITIVE_NUMBER
    if isinstance(value, list):
        return [
            _normalize_platform_sensitive_payload(item, (*path, index))
            for index, item in enumerate(value)
        ]
    if isinstance(value, dict):
        if path == DEPENDENCY_VERSIONS_PATH and any(
            value == profile for profile in APPROVED_DEPENDENCY_VERSION_PROFILES
        ):
            return dict(CANONICAL_DEPENDENCY_VERSIONS)
        normalized = {}
        for key, item in value.items():
            if (
                key.endswith("_identity")
                and isinstance(item, str)
                and item.startswith("sha256-stable-json-v1:")
            ):
                normalized[key] = PLATFORM_SENSITIVE_IDENTITY
            else:
                normalized[key] = _normalize_platform_sensitive_payload(
                    item, (*path, key)
                )
        fingerprint_fields = {
            "exact_hash_algorithm",
            "exact_hash_value",
            "tolerance_hash_algorithm",
            "tolerance_hash_value",
        }
        if fingerprint_fields <= normalized.keys():
            normalized["exact_hash_value"] = PLATFORM_SENSITIVE_HASH
            normalized["tolerance_hash_value"] = PLATFORM_SENSITIVE_HASH
        return normalized
    return value


def _assert_correction_numerics(payload: dict[str, object]) -> None:
    """Compare scientific correction outputs independently of BLAS byte identity."""

    corrected_matrix = cast(
        dict[str, object],
        _payload_value(payload, (*CORRECTION_CASE_PATH, "corrected_matrix")),
    )
    actual_rows = cast(list[list[float]], corrected_matrix["values"])
    assert [len(row) for row in actual_rows] == [4, 4, 4, 4]
    actual_values = [value for row in actual_rows for value in row]
    expected_values = [
        value for row in EXPECTED_CORRECTED_MATRIX_VALUES for value in row
    ]
    assert actual_values == pytest.approx(expected_values, rel=1e-7, abs=1e-8)

    executor = cast(
        dict[str, object],
        _payload_value(payload, (*CORRECTION_CASE_PATH, "diagnostics", "executor")),
    )
    for path, expected in EXPECTED_CORRECTION_DIAGNOSTICS.items():
        value: object = executor
        for key in path:
            value = cast(dict[str, object], value)[key]
        assert isinstance(value, float)
        assert value == pytest.approx(expected, rel=1e-7, abs=1e-8)
    singular_values = cast(list[float], executor["singular_values"])
    assert singular_values == pytest.approx([2.23606797749979], rel=1e-7, abs=1e-8)


def _assert_ranking_numerics(payload: dict[str, object]) -> None:
    """Compare SPS ranking outputs independently of platform math libraries."""

    cases = cast(dict[str, object], payload["cases"])
    for case_name, expected in RANKING_NUMERIC_BASELINES.items():
        case = cast(dict[str, object], cases[case_name])
        ranking = cast(list[dict[str, object]], case["site_ranking"])
        consensus = [
            cast(float, record["consensus_stability_score"]) for record in ranking
        ]
        statistics = [
            cast(dict[str, object], statistic)
            for record in ranking
            for statistic in cast(list[object], record["dataset_statistics"])
        ]
        quantiles = [
            cast(float, statistic["rank_quantile"])
            for statistic in statistics
            if "rank_quantile" in statistic
        ]
        stability = [
            cast(float, statistic["stability_score"]) for statistic in statistics
        ]
        assert consensus == pytest.approx(expected[0], rel=1e-7, abs=1e-8)
        assert quantiles == pytest.approx(expected[1], rel=1e-7, abs=1e-8)
        assert stability == pytest.approx(expected[2], rel=1e-7, abs=1e-8)


def _assert_regression_payload_match(
    payload: dict[str, object], *, actual_python_version: str
) -> None:
    """Check stable structure exactly and correction numerics within tolerance."""

    actual_provenance = _payload_value(payload, PYTHON_VERSION_PROVENANCE_PATH)
    assert actual_provenance == actual_python_version
    actual_parts = actual_python_version.split(".")
    assert len(actual_parts) == 3
    assert tuple(int(part) for part in actual_parts[:2]) == sys.version_info[:2]
    dependency_versions = _payload_value(payload, DEPENDENCY_VERSIONS_PATH)
    assert dependency_versions in APPROVED_DEPENDENCY_VERSION_PROFILES
    cases = cast(dict[str, object], payload["cases"])
    assert payload["cases_sha256"] == _canonical_sha256(cases)
    _assert_ranking_numerics(payload)
    _assert_correction_numerics(payload)

    qualified = cast(
        dict[str, object],
        _normalize_platform_sensitive_payload(deepcopy(payload)),
    )
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


def test_control_site_regression_permits_svd_factor_orientation_difference() -> None:
    payload = _fresh_process_comparison_payload()
    fingerprint = cast(
        dict[str, object],
        _payload_value(
            payload,
            (
                "cases",
                "correction_from_restored_controls",
                "provenance",
                "resolved_parameters",
                "executor",
                "estimated_unwanted_factors_fingerprint",
            ),
        ),
    )
    fingerprint["exact_hash_value"] = "opposite-svd-orientation-exact"
    fingerprint["tolerance_hash_value"] = "opposite-svd-orientation-tolerance"
    payload["cases_sha256"] = _canonical_sha256(payload["cases"])

    _assert_regression_payload_match(
        payload,
        actual_python_version=platform.python_version(),
    )


@pytest.mark.parametrize(
    "change",
    ["scientific", "dependency", "ranking", "corrected_matrix", "major_minor"],
)
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
    elif change == "corrected_matrix":
        corrected_matrix = cast(
            dict[str, object],
            _payload_value(
                payload,
                (
                    "cases",
                    "correction_from_restored_controls",
                    "corrected_matrix",
                ),
            ),
        )
        values = cast(list[list[float]], corrected_matrix["values"])
        values[0][0] += 0.1
    elif change == "ranking":
        ranking = cast(
            list[dict[str, object]],
            cast(
                dict[str, object], cast(dict[str, object], payload["cases"])["partial"]
            )["site_ranking"],
        )
        ranking[0]["consensus_stability_score"] = (
            cast(float, ranking[0]["consensus_stability_score"]) + 0.1
        )
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
