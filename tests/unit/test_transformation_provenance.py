from __future__ import annotations

from collections.abc import Iterator, Mapping
from typing import cast

import numpy as np
import pandas as pd
import pytest

from phospy.errors.transformations import InvalidTransformationStateError
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.models import JsonValue
from phospy.science.datasets.builders.transformation_resolver import (
    DatasetIntensityScaleResolver,
)
from phospy.science.transformations._authority import (
    dataset_quantitative_meaning_transition_authority,
    dataset_resolver_establishment_authority,
)
from phospy.science.transformations.models import (
    QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE,
    IntensityScaleEstablishmentMode,
    IntensityScaleEstablishmentProvenance,
    IntensityScaleEstablishmentSource,
    IntensityScaleEvidenceLevel,
    IntensityScaleState,
    MatrixIntensityScaleState,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
    QuantitativeMeaningTransitionProvenance,
)
from phospy.science.transformations.transformers import IdentityTransformer


def _phospho() -> pd.DataFrame:
    return pd.DataFrame(
        {"sample_a": [1.0, 2.0]},
        index=pd.Index(["MAPK14;Y182;", "GSK3B;S9;"], name="site_id"),
    )


class _DuplicateKeyMapping(Mapping[str, object]):
    def __iter__(self) -> Iterator[str]:
        return iter(("duplicate", "duplicate"))

    def __len__(self) -> int:
        return 2

    def __getitem__(self, key: str) -> object:
        if key == "duplicate":
            return "value"
        raise KeyError(key)


def _scale_provenance(
    parameters: Mapping[object, object],
) -> IntensityScaleEstablishmentProvenance:
    return IntensityScaleEstablishmentProvenance(
        scale="log2",
        mode=IntensityScaleEstablishmentMode.TRANSFORMED,
        source=IntensityScaleEstablishmentSource.TRANSFORMED_BY_PHOSPY,
        evidence_level=IntensityScaleEvidenceLevel.OBSERVED_TRANSFORMATION,
        transformer_name="tests.transformer",
        parameters=cast(Mapping[str, object], parameters),
        trace_id="trace-1",
        diagnostic_warnings=("review scale evidence",),
    )


def _scale_provenance_from_payload(
    payload: Mapping[str, object],
) -> IntensityScaleEstablishmentProvenance:
    return IntensityScaleEstablishmentProvenance(
        scale=cast(str, payload["scale"]),
        mode=cast(str, payload["establishment_mode"]),
        source=cast(str, payload["establishment_source"]),
        evidence_level=cast(str, payload["evidence_level"]),
        transformer_name=cast(str | None, payload["transformer_name"]),
        input_declaration_source=cast(
            str | None,
            payload["input_declaration_source"],
        ),
        parameters=cast(Mapping[str, object], payload["parameters"]),
        trace_id=cast(str | None, payload["trace_id"]),
        diagnostic_warnings=tuple(cast(list[str], payload["diagnostic_warnings"])),
    )


def _meaning_provenance(
    parameters: Mapping[object, object],
) -> QuantitativeMeaningTransitionProvenance:
    return QuantitativeMeaningTransitionProvenance(
        source_quantity=None,
        target_quantity=QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        operation_id=QUANTITATIVE_MEANING_OPERATION_SCALE_CONTRACT_INFERENCE,
        producer_id="tests.transformer",
        evidence_mode=QuantitativeMeaningEvidenceMode.INFERRED_FROM_SCALE_CONTRACT,
        parameters=cast(Mapping[str, object], parameters),
        diagnostic_caveat_codes=("quantitative_meaning_test",),
    )


def test_establishment_provenance_payload_exposes_mode_and_scale() -> None:
    declared_linear = IntensityScaleState(
        phospho=MatrixIntensityScaleState.linear(established_by="trusted.input"),
        total=None,
    )
    state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=_phospho(),
            total=None,
            declared_input_scale_state=declared_linear,
            input_declaration_source="tests.unit",
        )
        .intensity_scale_state
    )
    provenance = state.establishment_provenance

    assert provenance is not None
    payload = provenance.to_payload()
    assert payload["scale"] == "linear"
    assert (
        payload["establishment_mode"] == IntensityScaleEstablishmentMode.DECLARED.value
    )
    assert (
        payload["establishment_source"]
        == IntensityScaleEstablishmentSource.DECLARED_BY_USER.value
    )
    assert payload["diagnostic_warnings"] == []


def test_with_quantitative_meaning_is_blocked() -> None:
    declared_linear = IntensityScaleState(
        phospho=MatrixIntensityScaleState.linear(established_by="trusted.input"),
        total=None,
    )
    state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=_phospho(),
            total=None,
            declared_input_scale_state=declared_linear,
            input_declaration_source="tests.unit",
        )
        .intensity_scale_state
    )

    with pytest.raises(InvalidTransformationStateError, match="no longer supported"):
        state.with_quantitative_meaning(QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE)


def test_quantitative_meaning_transition_has_separate_provenance() -> None:
    declared_linear = IntensityScaleState(
        phospho=MatrixIntensityScaleState.linear(established_by="trusted.input"),
        total=None,
    )
    state = (
        DatasetIntensityScaleResolver(transformer=IdentityTransformer())
        .run(
            phospho=_phospho(),
            total=None,
            declared_input_scale_state=declared_linear,
            input_declaration_source="tests.unit",
        )
        .intensity_scale_state
    )
    provenance = _meaning_provenance({"scale_kind": "linear"})

    updated = state.transition_quantitative_meaning(
        target_quantity=QuantitativeMeaning.PHOSPHOSITE_ABUNDANCE,
        provenance=provenance,
        authority=dataset_quantitative_meaning_transition_authority(),
    )

    assert updated.establishment_mode is IntensityScaleEstablishmentMode.DECLARED
    assert updated.establishment_provenance == state.establishment_provenance
    assert updated.quantitative_meaning_provenance == provenance
    assert updated.quantitative_meaning_provenance != updated.establishment_provenance


def test_intensity_scale_establishment_parameters_are_recursively_immutable() -> None:
    source: dict[object, object] = {
        "nested_list": ["phospho", {"thresholds": [0.25, 0.75]}],
        "nested_mapping": {"labels": ["linear"], "details": {"version": 1}},
    }

    provenance = _scale_provenance(source)

    cast(list[object], source["nested_list"]).append("source-only")
    cast(
        list[object],
        cast(dict[str, object], source["nested_mapping"])["labels"],
    ).append("source-only")
    cast(
        dict[str, object],
        cast(dict[str, object], source["nested_mapping"])["details"],
    )["version"] = 2

    payload = provenance.to_payload()
    assert payload["parameters"] == {
        "nested_list": ["phospho", {"thresholds": [0.25, 0.75]}],
        "nested_mapping": {"labels": ["linear"], "details": {"version": 1}},
    }


def test_quantitative_meaning_parameters_are_recursively_immutable() -> None:
    source: dict[object, object] = {
        "nested_list": ["meaning", {"thresholds": [0.25, 0.75]}],
        "nested_mapping": {"labels": ["abundance"], "details": {"version": 1}},
    }

    provenance = _meaning_provenance(source)

    cast(list[object], source["nested_list"]).append("source-only")
    cast(
        list[object],
        cast(dict[str, object], source["nested_mapping"])["labels"],
    ).append("source-only")
    cast(
        dict[str, object],
        cast(dict[str, object], source["nested_mapping"])["details"],
    )["version"] = 2

    payload = provenance.to_payload()
    assert payload["parameters"] == {
        "nested_list": ["meaning", {"thresholds": [0.25, 0.75]}],
        "nested_mapping": {"labels": ["abundance"], "details": {"version": 1}},
    }


def test_intensity_scale_payload_is_fresh_and_hash_stable() -> None:
    provenance = _scale_provenance(
        {
            "nested_list": ["phospho", {"thresholds": [0.25, 0.75]}],
            "nested_mapping": {"labels": ["linear"], "details": {"version": 1}},
        }
    )

    first_payload = provenance.to_payload()
    first_hash = hash_json_payload(cast(JsonValue, first_payload))
    first_parameters = cast(dict[str, object], first_payload["parameters"])
    cast(list[object], first_parameters["nested_list"]).append("payload-only")
    cast(
        list[object],
        cast(dict[str, object], first_parameters["nested_mapping"])["labels"],
    ).append("payload-only")

    second_payload = provenance.to_payload()
    restored = _scale_provenance_from_payload(second_payload)

    assert second_payload["parameters"] == {
        "nested_list": ["phospho", {"thresholds": [0.25, 0.75]}],
        "nested_mapping": {"labels": ["linear"], "details": {"version": 1}},
    }
    assert hash_json_payload(cast(JsonValue, second_payload)) == first_hash
    assert restored.to_payload() == second_payload
    assert hash_json_payload(
        cast(JsonValue, restored.to_payload())
    ) == hash_json_payload(cast(JsonValue, second_payload))


def test_quantitative_meaning_payload_is_fresh_and_hash_stable() -> None:
    provenance = _meaning_provenance(
        {
            "nested_list": ["meaning", {"thresholds": [0.25, 0.75]}],
            "nested_mapping": {"labels": ["abundance"], "details": {"version": 1}},
        }
    )

    first_payload = provenance.to_payload()
    first_hash = hash_json_payload(cast(JsonValue, first_payload))
    first_parameters = cast(dict[str, object], first_payload["parameters"])
    cast(list[object], first_parameters["nested_list"]).append("payload-only")
    cast(
        list[object],
        cast(dict[str, object], first_parameters["nested_mapping"])["labels"],
    ).append("payload-only")

    second_payload = provenance.to_payload()
    restored = QuantitativeMeaningTransitionProvenance.from_payload(second_payload)

    assert second_payload["parameters"] == {
        "nested_list": ["meaning", {"thresholds": [0.25, 0.75]}],
        "nested_mapping": {"labels": ["abundance"], "details": {"version": 1}},
    }
    assert hash_json_payload(cast(JsonValue, second_payload)) == first_hash
    assert restored.to_payload() == second_payload
    assert hash_json_payload(
        cast(JsonValue, restored.to_payload())
    ) == hash_json_payload(cast(JsonValue, second_payload))


def test_intensity_scale_establishment_rejects_invalid_json_state() -> None:
    with pytest.raises(InvalidTransformationStateError, match="keys must be strings"):
        _scale_provenance({1: "numeric-key"})

    with pytest.raises(
        InvalidTransformationStateError,
        match="duplicate JSON object key",
    ):
        _scale_provenance(_DuplicateKeyMapping())

    with pytest.raises(
        InvalidTransformationStateError,
        match="duplicate JSON object key",
    ):
        establishable_state = IntensityScaleState.raw()
        establishable_state.with_establishment(
            established_via="tests.establishment",
            authority=dataset_resolver_establishment_authority(),
            establishment_mode=IntensityScaleEstablishmentMode.DERIVED,
            evidence_level=IntensityScaleEvidenceLevel.UNKNOWN,
            transformer_name=None,
            input_declaration_source=None,
            parameters=_DuplicateKeyMapping(),
            trace_id=None,
            diagnostic_warnings=(),
        )

    invalid_values = (
        {"bad": float("nan")},
        {"bad": float("inf")},
        {"bad": {"unsupported"}},
        {"bad": np.array([1.0])},
    )
    for parameters in invalid_values:
        with pytest.raises(InvalidTransformationStateError):
            _scale_provenance(parameters)
