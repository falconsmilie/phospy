from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
from threadpoolctl import threadpool_info

import tools.architecture.capture_batch_correction_regression as comparison
from tools.architecture.capture_batch_correction_regression import (
    build_comparison_payload,
)

pytestmark = [pytest.mark.release_gate, pytest.mark.reproducibility]

EXPECTED_CASES_SHA256 = (
    "8426feab0e287fad266ed7a09746e5edf5dc2f7687752539c8604de6bfe24326"
)
EXPECTED_PROVENANCE_COMPUTE_ENVIRONMENT = {
    "platform_system": "Windows",
    "platform_machine": "AMD64",
    "python_distribution": "Anaconda, Inc.",
    "openblas_coretype_env": "SkylakeX",
    "openblas_architectures": ("SkylakeX",),
}
NUMERICAL_RTOL = 1e-10
NUMERICAL_ATOL = 1e-12
NUMERICAL_EXPECTATIONS_PATH = (
    Path(__file__).parents[1]
    / "fixtures"
    / "architecture"
    / "batch_correction_regression_frames.json"
)
PROVENANCE_EXPECTATIONS_PATH = (
    Path(__file__).parents[1]
    / "fixtures"
    / "architecture"
    / "batch_correction_provenance_regression.json"
)
EXPECTED_WORKFLOW_CAVEATS = {
    "sps_ruv_style": [],
    "sps_ruv_style_permuted": [],
    "ruv_iii_style": [],
}


def test_regression_harness_uses_exact_provenance_owners() -> None:
    from phospy.preparation_records.batch_correction import (
        BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
        BATCH_CORRECTION_STATUS_DISABLED,
        BatchCorrectionDiagnostics,
        BatchCorrectionPolicy,
        BatchCorrectionReport,
    )
    from phospy.preprocessing.provenance.models import (
        PreprocessingStageProvenance as CurrentPreprocessingStageProvenance,
    )
    from phospy.provenance.models import EnvironmentProvenance
    from phospy.run_lineage.models import RunProvenance as CurrentRunProvenance
    from phospy.run_lineage.serialization import (
        from_payload as current_from_payload,
    )
    from phospy.run_lineage.serialization import to_payload as current_to_payload

    assert comparison.EnvironmentProvenance is EnvironmentProvenance
    assert (
        comparison.PreprocessingStageProvenance is CurrentPreprocessingStageProvenance
    )
    assert comparison.RunProvenance is CurrentRunProvenance
    assert comparison.from_payload is current_from_payload
    assert comparison.to_payload is current_to_payload
    assert (
        comparison.BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE
        == BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE
    )
    assert (
        comparison.BATCH_CORRECTION_STATUS_DISABLED == BATCH_CORRECTION_STATUS_DISABLED
    )
    assert comparison.BatchCorrectionDiagnostics is BatchCorrectionDiagnostics
    assert comparison.BatchCorrectionPolicy is BatchCorrectionPolicy
    assert comparison.BatchCorrectionReport is BatchCorrectionReport


def test_batch_correction_matches_regression_invariants() -> None:
    payload = build_comparison_payload()

    assert payload["cases_sha256"] == EXPECTED_CASES_SHA256


def test_batch_correction_matches_regression_frames() -> None:
    payload = build_comparison_payload()
    cases = cast(dict[str, Any], payload["cases"])
    expected = json.loads(NUMERICAL_EXPECTATIONS_PATH.read_text(encoding="utf-8"))
    actual = _numerical_frame_values(cases)

    assert tuple(actual) == tuple(expected)
    for path, expected_values in expected.items():
        np.testing.assert_allclose(
            np.asarray(actual[path], dtype=float),
            np.asarray(expected_values, dtype=float),
            rtol=NUMERICAL_RTOL,
            atol=NUMERICAL_ATOL,
            equal_nan=True,
            err_msg=f"batch-correction regression frame changed: {path}",
        )


def test_batch_correction_preserves_ordered_caveats() -> None:
    payload = build_comparison_payload()
    cases = cast(dict[str, Any], payload["cases"])

    assert {
        case_name: cases[case_name]["caveats"]
        for case_name in EXPECTED_WORKFLOW_CAVEATS
    } == EXPECTED_WORKFLOW_CAVEATS


def test_batch_correction_provenance_bytes_match_pinned_regression() -> None:
    expected = json.loads(PROVENANCE_EXPECTATIONS_PATH.read_text(encoding="utf-8"))
    runtime_versions = {
        "python": platform.python_version(),
        "numpy": importlib.metadata.version("numpy"),
        "pandas": importlib.metadata.version("pandas"),
        "scipy": importlib.metadata.version("scipy"),
        "scikit-learn": importlib.metadata.version("scikit-learn"),
    }
    openblas_architectures = tuple(
        sorted(
            {
                str(record.get("architecture"))
                for record in threadpool_info()
                if record.get("internal_api") == "openblas"
            }
        )
    )
    compute_environment = {
        "platform_system": platform.system(),
        "platform_machine": platform.machine(),
        "python_distribution": (
            "Anaconda, Inc." if "packaged by Anaconda, Inc." in sys.version else None
        ),
        "openblas_coretype_env": os.environ.get("OPENBLAS_CORETYPE"),
        "openblas_architectures": openblas_architectures,
    }
    if (
        runtime_versions != expected["runtime"]
        or compute_environment != EXPECTED_PROVENANCE_COMPUTE_ENVIRONMENT
    ):
        pytest.skip(
            "exact provenance bytes require the documented pinned "
            "lane: "
            f"expected_runtime={expected['runtime']!r}, "
            f"actual_runtime={runtime_versions!r}, "
            "expected_compute_environment="
            f"{EXPECTED_PROVENANCE_COMPUTE_ENVIRONMENT!r}, "
            f"actual_compute_environment={compute_environment!r}"
        )

    payload = build_comparison_payload()
    cases = cast(dict[str, Any], payload["cases"])
    for case_name, case_expected in expected["direct_payloads"].items():
        case = cases[case_name]
        provenance = case["direct_provenance"]
        provenance_bytes = case["direct_provenance_bytes"].encode("utf-8")
        assert len(provenance_bytes) == case_expected["byte_length"]
        assert _sha256(provenance_bytes) == case_expected["sha256"]
        assert tuple(provenance) == tuple(expected["direct_payload_ordered_fields"])
        assert all(
            field not in provenance
            for field in expected["direct_payload_omitted_optional_fields"]
        )
        assert provenance["schema_version"] == case_expected["schema_version"]
        assert provenance["requested_method"] == case_expected["requested_method"]
        assert (
            provenance["control_site_source"]["source_type"]
            == case_expected["control_source_type"]
        )
        assert provenance["phospy_version"] == case_expected["phospy_version"]
        assert provenance["python_version"] == case_expected["python_version"]
        assert _fingerprint_values(provenance) == {
            key: case_expected[key]
            for key in (
                "input_exact_hash",
                "input_tolerance_hash",
                "output_exact_hash",
                "output_tolerance_hash",
            )
        }

    nested = cases["nested_provenance"]
    nested_expected = expected["nested_run_provenance"]
    nested_payload = nested["payload"]
    nested_bytes = nested["payload_bytes"].encode("utf-8")
    assert len(nested_bytes) == nested_expected["byte_length"]
    assert _sha256(nested_bytes) == nested_expected["sha256"]
    assert tuple(nested_payload) == tuple(nested_expected["ordered_fields"])
    assert nested_payload["workflow_name"] == nested_expected["workflow_name"]
    assert (
        nested_payload["environment"]["package_name"]
        == nested_expected["environment_package_name"]
    )
    assert (
        nested_payload["environment"]["package_version"]
        == nested_expected["environment_package_version"]
    )
    stage = nested_payload["preprocessing_stages"][0]
    assert stage["stage"] == nested_expected["stage"]
    assert stage["operation"] == nested_expected["operation"]
    assert (
        stage["batch_correction_provenance"]["requested_method"]
        == nested_expected["nested_requested_method"]
    )
    assert (
        nested_payload["input_tables"][0]["exact_hash_value"]
        == (nested_expected["input_exact_hash"])
    )
    assert (
        nested_payload["output_tables"][0]["exact_hash_value"]
        == (nested_expected["output_exact_hash"])
    )
    assert stage["input_hash"] == nested_expected["stage_input_hash"]
    assert stage["output_hash"] == nested_expected["stage_output_hash"]


def test_matrix_missingness_and_persistence_contracts_are_complete() -> None:
    payload = build_comparison_payload()
    cases = cast(dict[str, Any], payload["cases"])

    assert cases["disabled_report"]["status"] == "disabled"
    assert cases["linear_residualisation"]["input_non_mutation"] is True
    assert (
        cases["linear_residualisation"]["defensive_dataframe_access_preserved"] is True
    )
    for case_name in (
        "sps_ruv_style",
        "sps_ruv_style_permuted",
        "ruv_iii_style",
    ):
        case = cases[case_name]
        assert case["input_non_mutation"] is True
        assert case["defensive_dataframe_access_preserved"] is True
        assert case["report"]["status"] == "applied"
        assert case["caveats"] == EXPECTED_WORKFLOW_CAVEATS[case_name]
        assert case["direct_provenance"]["schema_version"] == 1
        assert case["direct_provenance_bytes"] == json.dumps(
            case["direct_provenance"],
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
        )

    missingness = cases["missingness_preparation"]
    assert missingness["actual_missing"]["temporary_completion_applied"] is True
    assert missingness["actual_missing"]["input_non_mutation"] is True
    assert missingness["upstream_imputed"]["temporary_completion_applied"] is False
    assert missingness["actual_missing"]["restored"]["values"][1][1] is None

    nested = cases["nested_provenance"]
    assert nested["payload"] == nested["round_trip"]
    assert nested["runtime_identity_preserved"] is True
    assert nested["payload_bytes"] == json.dumps(
        nested["payload"],
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    )

    errors = cases["supported_errors"]
    assert errors["ungoverned_missing_cell"]["type"] == "PhosPyInputError"
    assert errors["malformed_historical_payload"]["type"] == "PhosPyInputError"


def _numerical_frame_values(
    value: object, path: tuple[str, ...] = ()
) -> dict[str, object]:
    frames: dict[str, object] = {}
    if isinstance(value, dict):
        dtypes = value.get("dtypes")
        values = value.get("values")
        if (
            isinstance(dtypes, list)
            and dtypes
            and all(
                isinstance(dtype, str)
                and (dtype.startswith("float") or dtype.startswith("int"))
                for dtype in dtypes
            )
            and isinstance(values, list)
        ):
            frames["/".join(path)] = values
            return frames
        for key, item in value.items():
            frames.update(_numerical_frame_values(item, (*path, str(key))))
    elif isinstance(value, list):
        for position, item in enumerate(value):
            frames.update(_numerical_frame_values(item, (*path, str(position))))
    return frames


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _fingerprint_values(provenance: dict[str, Any]) -> dict[str, str]:
    input_fingerprint = provenance["input_matrix_fingerprint"]
    output_fingerprint = provenance["output_matrix_fingerprint"]
    return {
        "input_exact_hash": input_fingerprint["exact_hash_value"],
        "input_tolerance_hash": input_fingerprint["tolerance_hash_value"],
        "output_exact_hash": output_fingerprint["exact_hash_value"],
        "output_tolerance_hash": output_fingerprint["tolerance_hash_value"],
    }
