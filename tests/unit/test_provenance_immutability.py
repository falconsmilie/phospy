from __future__ import annotations

from typing import cast

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.models import (
    ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    EnvironmentProvenance,
    JsonValue,
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.serialization import from_payload, to_payload


def _fingerprint() -> TableFingerprint:
    return TableFingerprint(
        name="dataset.phospho",
        rows=1,
        columns=1,
        index_name="site_id",
        column_names=("sample_a",),
        dtypes=("float64",),
        exact_hash_algorithm="sha256-stable-json-v1",
        exact_hash_value="a" * 64,
        tolerance_hash_algorithm="sha256-float-round-8dp-v1",
        tolerance_hash_value="b" * 64,
        index_structure={"axis": {"labels": ["site_a"]}},
        column_index_structure={"axis": {"labels": ["sample_a"]}},
    )


def _run_with_parameters(parameters: dict[str, object]) -> RunProvenance:
    fingerprint = _fingerprint()
    stage = PreprocessingStageProvenance(
        stage="missing_data",
        operation="forbid",
        parameters=parameters,
        input_shape=(1, 1),
        output_shape=(1, 1),
        input_hash=fingerprint.tolerance_hash_value,
        output_hash=fingerprint.tolerance_hash_value,
        phospho_input_hash=fingerprint.tolerance_hash_value,
        phospho_output_hash=fingerprint.tolerance_hash_value,
        dropped_row_ids=(),
        dropped_row_count=0,
        schema_version=PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
        consumed_input_tables=(fingerprint,),
        produced_output_tables=(fingerprint,),
    )
    return RunProvenance(
        environment=EnvironmentProvenance(
            package_name="phospy",
            package_version="0",
            python_version="3.13",
            dependency_versions={"pandas": "2.0"},
            schema_version=ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
            platform={"system": "test"},
            blas_lapack={"backend": {"libraries": ["openblas"]}},
        ),
        input_tables=(fingerprint,),
        preprocessing_stages=(stage,),
        reference=None,
        workflow_name="immutability_test",
        workflow_parameters={"parameters": parameters},
        random_state=None,
        random_seed_policy=None,
        output_tables=(fingerprint,),
    )


def test_provenance_constructor_recursively_freezes_source_input() -> None:
    source = {
        "policy": "test",
        "nested": {"items": ["a", {"score": 1.0}]},
    }

    provenance = _run_with_parameters(source)
    parameters = provenance.workflow_parameters["parameters"]
    assert parameters["nested"]["items"] == ("a", {"score": 1.0})

    with pytest.raises(TypeError):
        provenance.workflow_parameters["new"] = "value"  # type: ignore[index]
    with pytest.raises(TypeError):
        parameters["nested"]["extra"] = "value"  # type: ignore[index]
    with pytest.raises((AttributeError, TypeError)):
        parameters["nested"]["items"].append("b")  # type: ignore[attr-defined]
    with pytest.raises(TypeError):
        parameters["nested"]["items"][1]["score"] = 2.0  # type: ignore[index]

    source["nested"]["items"].append("source-only")
    source["nested"]["items"][1]["score"] = 9.0

    assert parameters["nested"]["items"] == ("a", {"score": 1.0})


def test_provenance_serialization_returns_fresh_mutable_payloads() -> None:
    provenance = _run_with_parameters({"nested": {"items": ["a", {"score": 1.0}]}})

    payload = to_payload(provenance)
    payload_parameters = payload["workflow_parameters"]["parameters"]
    payload_parameters["nested"]["items"][1]["score"] = 9.0
    payload_parameters["nested"]["items"].append("payload-only")

    fresh_payload = to_payload(provenance)
    fresh_parameters = fresh_payload["workflow_parameters"]["parameters"]

    assert fresh_parameters["nested"]["items"] == ["a", {"score": 1.0}]
    assert provenance.workflow_parameters["parameters"]["nested"]["items"] == (
        "a",
        {"score": 1.0},
    )


def test_provenance_rejects_unsupported_and_non_finite_json_values() -> None:
    with pytest.raises(PhosPyInputError, match="JSON-compatible"):
        _run_with_parameters({"bad": object()})

    with pytest.raises(PhosPyInputError, match="finite JSON number"):
        _run_with_parameters({"bad": float("nan")})

    with pytest.raises(PhosPyInputError, match="finite JSON number"):
        _run_with_parameters({"bad": float("inf")})


def test_provenance_round_trip_payload_and_hash_are_stable() -> None:
    provenance = _run_with_parameters({"nested": {"items": ["a", {"score": 1.0}]}})

    payload = to_payload(provenance)
    restored = from_payload(payload)
    restored_payload = to_payload(restored)

    assert restored_payload == payload
    assert hash_json_payload(cast(JsonValue, restored_payload)) == hash_json_payload(
        cast(JsonValue, payload)
    )
