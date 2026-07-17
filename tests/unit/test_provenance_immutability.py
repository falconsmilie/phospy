from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import cast

import pytest

from phospy.errors.input import PhosPyInputError
from phospy.provenance.derived_quantitative import (
    DerivedQuantitativeDataProvenance,
    DerivedSampleMapping,
)
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.immutability import FrozenJsonMapping
from phospy.provenance.models import (
    ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2,
    PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3,
    EnvironmentProvenance,
    JsonValue,
    PreprocessingStageProvenance,
    ReproducibilityCaveat,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.scientific_policy_models import (
    ScientificPolicyId,
    ScientificPolicyRecord,
)
from phospy.provenance.serialization import from_payload, to_payload


class _DuplicateKeyMapping(Mapping[str, object]):
    def __iter__(self) -> Iterator[str]:
        return iter(("duplicate", "duplicate"))

    def __len__(self) -> int:
        return 2

    def __getitem__(self, key: str) -> object:
        if key == "duplicate":
            return "value"
        raise KeyError(key)


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


def _run_with_parameters(parameters: Mapping[object, object]) -> RunProvenance:
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
    assert isinstance(parameters, FrozenJsonMapping)
    assert not isinstance(parameters, dict)
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


def test_base_class_mutation_bypass_is_unavailable() -> None:
    provenance = _run_with_parameters({"nested": {"items": ["a"]}})
    parameters = provenance.workflow_parameters["parameters"]
    nested = parameters["nested"]
    items = nested["items"]

    assert not isinstance(parameters, dict)
    assert isinstance(items, tuple)

    before_hash = hash_json_payload(cast(JsonValue, to_payload(provenance)))

    with pytest.raises(TypeError):
        dict.__setitem__(parameters, "new", "value")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        list.append(items, "b")  # type: ignore[arg-type]

    after_hash = hash_json_payload(cast(JsonValue, to_payload(provenance)))
    assert after_hash == before_hash
    assert provenance.workflow_parameters["parameters"]["nested"]["items"] == ("a",)


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


def test_provenance_rejects_non_string_and_colliding_json_object_keys() -> None:
    with pytest.raises(PhosPyInputError, match="keys must be strings"):
        _run_with_parameters({1: "numeric-key"})

    with pytest.raises(PhosPyInputError, match="duplicate JSON object key"):
        _run_with_parameters(_DuplicateKeyMapping())


def test_provenance_deserialization_rejects_non_string_root_keys() -> None:
    payload = to_payload(_run_with_parameters({"policy": "test"}))
    payload[1] = "numeric-key"

    with pytest.raises(PhosPyInputError, match="keys must be strings"):
        from_payload(cast(Mapping[str, object], payload))


def test_hashing_rejects_non_string_json_object_keys() -> None:
    with pytest.raises(PhosPyInputError, match="keys must be strings"):
        hash_json_payload(cast(JsonValue, {1: "numeric-key"}))

    with pytest.raises(PhosPyInputError, match="JSON-compatible"):
        hash_json_payload(cast(JsonValue, {"bad": object()}))

    with pytest.raises(PhosPyInputError, match="finite JSON number"):
        hash_json_payload(cast(JsonValue, {"bad": float("nan")}))


def test_provenance_model_families_use_same_immutable_json_container() -> None:
    nested_source = {"nested": {"items": ["a", {"score": 1.0}]}}
    fingerprint = _fingerprint()
    caveat = ReproducibilityCaveat(
        code="test_caveat",
        severity="warning",
        message="Test caveat.",
        details=nested_source,
    )
    policy = ScientificPolicyRecord(
        id=ScientificPolicyId.KINASE_PROFILE_SCORING,
        name="Test policy",
        version="1",
        description="Test policy.",
        parameters=nested_source,
        assumptions=("test",),
    )
    lineage = DerivedQuantitativeDataProvenance(
        derivation_type="technical_replicate_aggregation",
        parent_dataset_type="analysis_ready",
        derived_dataset_type="analysis_ready",
        parent_dataset_fingerprints=(fingerprint,),
        derived_dataset_fingerprints=(fingerprint,),
        sample_mapping=(
            DerivedSampleMapping(
                output_sample_id="sample_a",
                input_sample_ids=("sample_a",),
                condition="treated",
                biological_replicate_id="bio_1",
            ),
        ),
        aggregation_method="mean",
        input_intensity_scale="log2",
        output_intensity_scale="log2",
        quantitative_meaning="phosphosite_abundance",
        missingness_policy=nested_source,
        matrices_transformed={"phospho": True},
        implementation="test",
        implementation_version="1",
        parameters=nested_source,
    )

    for mapping in (
        caveat.details,
        policy.parameters,
        lineage.missingness_policy,
        lineage.matrices_transformed,
        lineage.parameters,
    ):
        assert isinstance(mapping, FrozenJsonMapping)
        assert not isinstance(mapping, dict)

    nested_source["nested"]["items"].append("source-only")
    assert caveat.details["nested"]["items"] == ("a", {"score": 1.0})
    assert policy.parameters["nested"]["items"] == ("a", {"score": 1.0})
    assert lineage.parameters["nested"]["items"] == ("a", {"score": 1.0})


def test_table_fingerprint_rejects_invalid_shape_and_hash_state() -> None:
    valid = {
        "name": "dataset.phospho",
        "rows": 1,
        "columns": 1,
        "index_name": "site_id",
        "column_names": ("sample_a",),
        "dtypes": ("float64",),
        "exact_hash_algorithm": "sha256-stable-json-v1",
        "exact_hash_value": "a" * 64,
        "tolerance_hash_algorithm": "sha256-float-round-8dp-v1",
        "tolerance_hash_value": "b" * 64,
    }

    with pytest.raises(PhosPyInputError, match="rows"):
        TableFingerprint(**{**valid, "rows": -1})
    with pytest.raises(PhosPyInputError, match="column_names length"):
        TableFingerprint(**{**valid, "columns": 2})
    with pytest.raises(PhosPyInputError, match="exact_hash_value"):
        TableFingerprint(**{**valid, "exact_hash_value": ""})


def test_provenance_round_trip_payload_and_hash_are_stable() -> None:
    provenance = _run_with_parameters({"nested": {"items": ["a", {"score": 1.0}]}})

    payload = to_payload(provenance)
    restored = from_payload(payload)
    restored_payload = to_payload(restored)

    assert restored_payload == payload
    assert hash_json_payload(cast(JsonValue, restored_payload)) == hash_json_payload(
        cast(JsonValue, payload)
    )
