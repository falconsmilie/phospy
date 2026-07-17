"""Typed provenance for derived quantitative dataset objects."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import (
    DEFAULT_STABLE_JSON_HASH_ALGORITHM,
    hash_json_payload,
)
from phospy.provenance.immutability import (
    freeze_json_mapping,
    freeze_json_value,
    thaw_json_mapping,
)
from phospy.provenance.models import (
    EnvironmentProvenance,
    JsonValue,
    ReferenceContextProtocol,
    RunProvenance,
    TableFingerprint,
)
from phospy.provenance.serialization import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)

DERIVED_QUANTITATIVE_DATA_PROVENANCE_SCHEMA_VERSION_V1 = 1
DERIVED_QUANTITATIVE_DATA_FINGERPRINT_ALGORITHM = DEFAULT_STABLE_JSON_HASH_ALGORITHM
TECHNICAL_REPLICATE_AGGREGATION_DERIVATION_TYPE = "technical_replicate_aggregation"
TECHNICAL_REPLICATE_AGGREGATOR_IMPLEMENTATION = (
    "phospy.workflows.differential.replicates.TechnicalReplicateAggregator"
)


@dataclass(frozen=True, slots=True)
class DerivedSampleMapping:
    """Source-to-derived sample mapping for one derived sample."""

    output_sample_id: str
    input_sample_ids: tuple[str, ...]
    condition: str
    biological_replicate_id: str
    technical_replicate_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "output_sample_id",
            _required_text(
                self.output_sample_id,
                field_name="derived_sample_mapping.output_sample_id",
            ),
        )
        input_sample_ids = _required_text_tuple(
            self.input_sample_ids,
            field_name="derived_sample_mapping.input_sample_ids",
        )
        if not input_sample_ids:
            raise PhosPyInputError(
                "derived_sample_mapping.input_sample_ids must not be empty"
            )
        object.__setattr__(self, "input_sample_ids", input_sample_ids)
        object.__setattr__(
            self,
            "condition",
            _required_text(
                self.condition, field_name="derived_sample_mapping.condition"
            ),
        )
        object.__setattr__(
            self,
            "biological_replicate_id",
            _required_text(
                self.biological_replicate_id,
                field_name="derived_sample_mapping.biological_replicate_id",
            ),
        )
        object.__setattr__(
            self,
            "technical_replicate_ids",
            _text_tuple(
                self.technical_replicate_ids,
                field_name="derived_sample_mapping.technical_replicate_ids",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-safe sample lineage payload."""

        return {
            "output_sample_id": self.output_sample_id,
            "input_sample_ids": list(self.input_sample_ids),
            "source_sample_ids": list(self.input_sample_ids),
            "condition": self.condition,
            "biological_replicate_id": self.biological_replicate_id,
            "technical_replicate_ids": list(self.technical_replicate_ids),
            "n_source_samples": int(len(self.input_sample_ids)),
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> DerivedSampleMapping:
        """Deserialize one sample mapping from a decoded payload."""

        payload = _require_mapping(payload, field_name="derived_sample_mapping")
        return cls(
            output_sample_id=_require_payload_str(
                payload.get("output_sample_id"),
                field_name="derived_sample_mapping.output_sample_id",
            ),
            input_sample_ids=tuple(
                _require_payload_str(
                    item,
                    field_name="derived_sample_mapping.input_sample_ids[]",
                )
                for item in _require_sequence(
                    payload.get("input_sample_ids"),
                    field_name="derived_sample_mapping.input_sample_ids",
                )
            ),
            condition=_require_payload_str(
                payload.get("condition"),
                field_name="derived_sample_mapping.condition",
            ),
            biological_replicate_id=_require_payload_str(
                payload.get("biological_replicate_id"),
                field_name="derived_sample_mapping.biological_replicate_id",
            ),
            technical_replicate_ids=tuple(
                _require_payload_str(
                    item,
                    field_name="derived_sample_mapping.technical_replicate_ids[]",
                )
                for item in _require_sequence(
                    payload.get("technical_replicate_ids", []),
                    field_name="derived_sample_mapping.technical_replicate_ids",
                )
            ),
        )


@dataclass(frozen=True, slots=True)
class DerivedQuantitativeDataProvenance:
    """Lineage for a quantitative dataset derived from another dataset."""

    derivation_type: str
    parent_dataset_type: str
    derived_dataset_type: str
    parent_dataset_fingerprints: tuple[TableFingerprint, ...]
    derived_dataset_fingerprints: tuple[TableFingerprint, ...]
    sample_mapping: tuple[DerivedSampleMapping, ...]
    aggregation_method: str
    input_intensity_scale: str
    output_intensity_scale: str
    quantitative_meaning: str
    missingness_policy: Mapping[str, JsonValue]
    matrices_transformed: Mapping[str, bool]
    implementation: str
    implementation_version: str
    schema_version: int = DERIVED_QUANTITATIVE_DATA_PROVENANCE_SCHEMA_VERSION_V1
    parameters: Mapping[str, JsonValue] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if (
            int(self.schema_version)
            != DERIVED_QUANTITATIVE_DATA_PROVENANCE_SCHEMA_VERSION_V1
        ):
            raise PhosPyInputError(
                "derived_quantitative_data.schema_version must be "
                f"{DERIVED_QUANTITATIVE_DATA_PROVENANCE_SCHEMA_VERSION_V1}"
            )
        object.__setattr__(
            self,
            "derivation_type",
            _required_text(
                self.derivation_type,
                field_name="derived_quantitative_data.derivation_type",
            ),
        )
        object.__setattr__(
            self,
            "parent_dataset_type",
            _required_text(
                self.parent_dataset_type,
                field_name="derived_quantitative_data.parent_dataset_type",
            ),
        )
        object.__setattr__(
            self,
            "derived_dataset_type",
            _required_text(
                self.derived_dataset_type,
                field_name="derived_quantitative_data.derived_dataset_type",
            ),
        )
        object.__setattr__(
            self,
            "parent_dataset_fingerprints",
            _table_fingerprint_tuple(
                self.parent_dataset_fingerprints,
                field_name="derived_quantitative_data.parent_dataset_fingerprints",
            ),
        )
        derived_fingerprints = _table_fingerprint_tuple(
            self.derived_dataset_fingerprints,
            field_name="derived_quantitative_data.derived_dataset_fingerprints",
        )
        object.__setattr__(self, "derived_dataset_fingerprints", derived_fingerprints)
        sample_mapping = _sample_mapping_tuple(self.sample_mapping)
        _validate_mapping_matches_derived_phospho(
            sample_mapping=sample_mapping,
            fingerprints=derived_fingerprints,
        )
        object.__setattr__(self, "sample_mapping", sample_mapping)
        object.__setattr__(
            self,
            "aggregation_method",
            _required_text(
                self.aggregation_method,
                field_name="derived_quantitative_data.aggregation_method",
            ),
        )
        object.__setattr__(
            self,
            "input_intensity_scale",
            _required_text(
                self.input_intensity_scale,
                field_name="derived_quantitative_data.input_intensity_scale",
            ),
        )
        object.__setattr__(
            self,
            "output_intensity_scale",
            _required_text(
                self.output_intensity_scale,
                field_name="derived_quantitative_data.output_intensity_scale",
            ),
        )
        object.__setattr__(
            self,
            "quantitative_meaning",
            _required_text(
                self.quantitative_meaning,
                field_name="derived_quantitative_data.quantitative_meaning",
            ),
        )
        object.__setattr__(
            self,
            "missingness_policy",
            _json_mapping(
                self.missingness_policy,
                field_name="derived_quantitative_data.missingness_policy",
            ),
        )
        object.__setattr__(
            self,
            "matrices_transformed",
            _bool_mapping(
                self.matrices_transformed,
                field_name="derived_quantitative_data.matrices_transformed",
            ),
        )
        object.__setattr__(
            self,
            "implementation",
            _required_text(
                self.implementation,
                field_name="derived_quantitative_data.implementation",
            ),
        )
        object.__setattr__(
            self,
            "implementation_version",
            _required_text(
                self.implementation_version,
                field_name="derived_quantitative_data.implementation_version",
            ),
        )
        object.__setattr__(
            self,
            "parameters",
            _json_mapping(
                self.parameters,
                field_name="derived_quantitative_data.parameters",
            ),
        )

    @property
    def parent_dataset_fingerprint_value(self) -> str:
        """Return a stable digest over all recorded parent table fingerprints."""

        return _fingerprint_collection_digest(self.parent_dataset_fingerprints)

    @property
    def derived_dataset_fingerprint_value(self) -> str:
        """Return a stable digest over all recorded derived table fingerprints."""

        return _fingerprint_collection_digest(self.derived_dataset_fingerprints)

    @property
    def lineage_hash_value(self) -> str:
        """Return a stable digest over the complete lineage payload."""

        return hash_json_payload(cast(JsonValue, self._payload_without_lineage_hash()))

    @property
    def sample_count(self) -> int:
        """Return the number of derived samples recorded in lineage."""

        return int(len(self.sample_mapping))

    def sample_groups(self) -> tuple[tuple[str, tuple[str, ...]], ...]:
        """Return replay-ready output-to-input sample groups."""

        return tuple(
            (item.output_sample_id, item.input_sample_ids)
            for item in self.sample_mapping
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-safe derived-lineage payload."""

        payload = self._payload_without_lineage_hash()
        payload["lineage_hash_algorithm"] = DEFAULT_STABLE_JSON_HASH_ALGORITHM
        payload["lineage_hash_value"] = self.lineage_hash_value
        return payload

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> DerivedQuantitativeDataProvenance:
        """Deserialize derived-lineage provenance from a decoded payload."""

        payload = _require_mapping(payload, field_name="derived_quantitative_data")
        schema_version = _require_int(
            payload.get("schema_version"),
            field_name="derived_quantitative_data.schema_version",
        )
        aggregation = _require_mapping(
            payload.get("aggregation"),
            field_name="derived_quantitative_data.aggregation",
        )
        parent_summary = _require_mapping(
            payload.get("parent_dataset_fingerprint"),
            field_name="derived_quantitative_data.parent_dataset_fingerprint",
        )
        derived_summary = _require_mapping(
            payload.get("derived_dataset_fingerprint"),
            field_name="derived_quantitative_data.derived_dataset_fingerprint",
        )
        restored = cls(
            schema_version=schema_version,
            derivation_type=_require_payload_str(
                payload.get("derivation_type"),
                field_name="derived_quantitative_data.derivation_type",
            ),
            parent_dataset_type=_require_payload_str(
                payload.get("parent_dataset_type"),
                field_name="derived_quantitative_data.parent_dataset_type",
            ),
            derived_dataset_type=_require_payload_str(
                payload.get("derived_dataset_type"),
                field_name="derived_quantitative_data.derived_dataset_type",
            ),
            parent_dataset_fingerprints=_fingerprints_from_summary(
                parent_summary,
                field_name="derived_quantitative_data.parent_dataset_fingerprint",
            ),
            derived_dataset_fingerprints=_fingerprints_from_summary(
                derived_summary,
                field_name="derived_quantitative_data.derived_dataset_fingerprint",
            ),
            sample_mapping=tuple(
                DerivedSampleMapping.from_payload(
                    _require_mapping(
                        item,
                        field_name=(
                            f"derived_quantitative_data.sample_mapping[{position}]"
                        ),
                    )
                )
                for position, item in enumerate(
                    _require_sequence(
                        payload.get("sample_mapping"),
                        field_name="derived_quantitative_data.sample_mapping",
                    )
                )
            ),
            aggregation_method=_require_payload_str(
                aggregation.get("method"),
                field_name="derived_quantitative_data.aggregation.method",
            ),
            input_intensity_scale=_require_payload_str(
                aggregation.get("input_intensity_scale"),
                field_name=(
                    "derived_quantitative_data.aggregation.input_intensity_scale"
                ),
            ),
            output_intensity_scale=_require_payload_str(
                aggregation.get("output_intensity_scale"),
                field_name=(
                    "derived_quantitative_data.aggregation.output_intensity_scale"
                ),
            ),
            quantitative_meaning=_require_payload_str(
                aggregation.get("quantitative_meaning"),
                field_name="derived_quantitative_data.aggregation.quantitative_meaning",
            ),
            missingness_policy=_json_mapping(
                _require_mapping(
                    payload.get("missingness_policy"),
                    field_name="derived_quantitative_data.missingness_policy",
                ),
                field_name="derived_quantitative_data.missingness_policy",
            ),
            matrices_transformed=_bool_mapping(
                _require_mapping(
                    payload.get("matrices_transformed"),
                    field_name="derived_quantitative_data.matrices_transformed",
                ),
                field_name="derived_quantitative_data.matrices_transformed",
            ),
            implementation=_require_payload_str(
                payload.get("implementation"),
                field_name="derived_quantitative_data.implementation",
            ),
            implementation_version=_require_payload_str(
                payload.get("implementation_version"),
                field_name="derived_quantitative_data.implementation_version",
            ),
            parameters=_json_mapping(
                _require_mapping(
                    payload.get("parameters", {}),
                    field_name="derived_quantitative_data.parameters",
                ),
                field_name="derived_quantitative_data.parameters",
            ),
        )
        expected_hash = _optional_payload_str(
            payload.get("lineage_hash_value"),
            field_name="derived_quantitative_data.lineage_hash_value",
        )
        if expected_hash is not None and expected_hash != restored.lineage_hash_value:
            raise PhosPyInputError(
                "derived_quantitative_data.lineage_hash_value does not match "
                "the decoded lineage payload"
            )
        return restored

    def _payload_without_lineage_hash(self) -> dict[str, object]:
        return {
            "schema_version": int(self.schema_version),
            "derivation_type": self.derivation_type,
            "parent_dataset_type": self.parent_dataset_type,
            "derived_dataset_type": self.derived_dataset_type,
            "implementation": self.implementation,
            "implementation_version": self.implementation_version,
            "aggregation": {
                "method": self.aggregation_method,
                "input_intensity_scale": self.input_intensity_scale,
                "output_intensity_scale": self.output_intensity_scale,
                "quantitative_meaning": self.quantitative_meaning,
            },
            "missingness_policy": thaw_json_mapping(
                self.missingness_policy,
                field_name="derived_quantitative_data.missingness_policy",
            ),
            "matrices_transformed": thaw_json_mapping(
                self.matrices_transformed,
                field_name="derived_quantitative_data.matrices_transformed",
            ),
            "parameters": thaw_json_mapping(
                self.parameters,
                field_name="derived_quantitative_data.parameters",
            ),
            "parent_dataset_fingerprint": _fingerprint_summary_payload(
                self.parent_dataset_fingerprints
            ),
            "derived_dataset_fingerprint": _fingerprint_summary_payload(
                self.derived_dataset_fingerprints
            ),
            "sample_mapping": [item.to_payload() for item in self.sample_mapping],
        }


def build_derived_quantitative_run_provenance(
    *,
    lineage: DerivedQuantitativeDataProvenance,
    environment: EnvironmentProvenance | None = None,
    reference_context: ReferenceContextProtocol | None = None,
) -> RunProvenance:
    """Build truthful run provenance for a derived quantitative dataset."""

    resolved_environment = (
        collect_environment_provenance() if environment is None else environment
    )
    return RunProvenance(
        environment=resolved_environment,
        input_tables=lineage.parent_dataset_fingerprints,
        preprocessing_stages=(),
        reference=None,
        workflow_name=lineage.derivation_type,
        workflow_parameters={
            "derived_quantitative_data": lineage.to_payload(),
            "construction": {
                "method": lineage.implementation,
                "dataset_type": lineage.derived_dataset_type,
                "parent_dataset_type": lineage.parent_dataset_type,
                "model_constructor": (
                    "DerivedAnalysisReadyPhosphoDataset._from_owned_derived_tables"
                ),
                "source_dataset_provenance_reused": False,
                "source_dataset_preprocessing_report_reused": False,
            },
        },
        random_state=None,
        random_seed_policy=None,
        output_tables=lineage.derived_dataset_fingerprints,
        scientific_policies=(),
        reference_context=reference_context,
    )


def _fingerprint_summary_payload(
    fingerprints: tuple[TableFingerprint, ...],
) -> dict[str, object]:
    return {
        "algorithm": DERIVED_QUANTITATIVE_DATA_FINGERPRINT_ALGORITHM,
        "value": _fingerprint_collection_digest(fingerprints),
        "tables": [table_fingerprint_to_payload(item) for item in fingerprints],
    }


def _fingerprint_collection_digest(
    fingerprints: tuple[TableFingerprint, ...],
) -> str:
    return hash_json_payload(
        cast(
            JsonValue,
            [table_fingerprint_to_payload(item) for item in fingerprints],
        )
    )


def _fingerprints_from_summary(
    payload: Mapping[str, object],
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    algorithm = _require_payload_str(
        payload.get("algorithm"),
        field_name=f"{field_name}.algorithm",
    )
    if algorithm != DERIVED_QUANTITATIVE_DATA_FINGERPRINT_ALGORITHM:
        raise PhosPyInputError(
            f"{field_name}.algorithm must be "
            f"{DERIVED_QUANTITATIVE_DATA_FINGERPRINT_ALGORITHM!r}"
        )
    fingerprints = tuple(
        table_fingerprint_from_payload(
            _require_mapping(item, field_name=f"{field_name}.tables[{position}]")
        )
        for position, item in enumerate(
            _require_sequence(payload.get("tables"), field_name=f"{field_name}.tables")
        )
    )
    expected_value = _require_payload_str(
        payload.get("value"),
        field_name=f"{field_name}.value",
    )
    observed_value = _fingerprint_collection_digest(fingerprints)
    if expected_value != observed_value:
        raise PhosPyInputError(f"{field_name}.value does not match table fingerprints")
    return _table_fingerprint_tuple(fingerprints, field_name=f"{field_name}.tables")


def _validate_mapping_matches_derived_phospho(
    *,
    sample_mapping: tuple[DerivedSampleMapping, ...],
    fingerprints: tuple[TableFingerprint, ...],
) -> None:
    phospho_fingerprint = _fingerprint_by_name(fingerprints, "dataset.phospho")
    if phospho_fingerprint is None:
        raise PhosPyInputError(
            "derived_quantitative_data.derived_dataset_fingerprints must include "
            "dataset.phospho"
        )
    expected_columns = tuple(item.output_sample_id for item in sample_mapping)
    if phospho_fingerprint.columns != len(expected_columns):
        raise PhosPyInputError(
            "derived_quantitative_data.sample_mapping length must match "
            "derived dataset.phospho column count"
        )
    if phospho_fingerprint.column_names != expected_columns:
        raise PhosPyInputError(
            "derived_quantitative_data.sample_mapping output_sample_id values "
            "must match derived dataset.phospho columns"
        )


def _fingerprint_by_name(
    fingerprints: tuple[TableFingerprint, ...],
    name: str,
) -> TableFingerprint | None:
    for fingerprint in fingerprints:
        if fingerprint.name == name:
            return fingerprint
    return None


def _table_fingerprint_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    fingerprints = tuple(values)
    if not fingerprints:
        raise PhosPyInputError(f"{field_name} must not be empty")
    for fingerprint in fingerprints:
        if not isinstance(fingerprint, TableFingerprint):
            raise PhosPyInputError(
                f"{field_name} must contain only TableFingerprint values"
            )
    return fingerprints


def _sample_mapping_tuple(
    values: object,
) -> tuple[DerivedSampleMapping, ...]:
    field_name = "derived_quantitative_data.sample_mapping"
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    sample_mapping = tuple(values)
    if not sample_mapping:
        raise PhosPyInputError(f"{field_name} must not be empty")
    for item in sample_mapping:
        if not isinstance(item, DerivedSampleMapping):
            raise PhosPyInputError(
                f"{field_name} must contain only DerivedSampleMapping values"
            )
    return sample_mapping


def _bool_mapping(values: object, *, field_name: str) -> Mapping[str, bool]:
    if not isinstance(values, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    result: dict[str, bool] = {}
    for key, value in values.items():
        if not isinstance(key, str):
            raise PhosPyInputError(
                f"{field_name} JSON object keys must be strings; "
                f"got {type(key).__name__}"
            )
        if key in result:
            raise PhosPyInputError(
                f"{field_name} contains duplicate JSON object key {key!r}"
            )
        if not isinstance(value, bool):
            raise PhosPyInputError(f"{field_name}.{key!r} must be a bool")
        result[key] = bool(value)
    return cast(
        Mapping[str, bool],
        freeze_json_mapping(result, field_name=field_name),
    )


def _json_mapping(values: object, *, field_name: str) -> Mapping[str, JsonValue]:
    if not isinstance(values, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    return cast(
        Mapping[str, JsonValue],
        freeze_json_mapping(values, field_name=field_name),
    )


def _to_json_value(value: object) -> JsonValue:
    return cast(
        JsonValue,
        freeze_json_value(value, field_name="derived_quantitative_data"),
    )


def _required_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def _text_tuple(values: object, *, field_name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    return tuple(
        _required_text(value, field_name=f"{field_name}[]") for value in values
    )


def _required_text_tuple(values: object, *, field_name: str) -> tuple[str, ...]:
    return _text_tuple(values, field_name=field_name)


def _require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if isinstance(value, Mapping):
        result: dict[str, object] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise PhosPyInputError(
                    f"{field_name} JSON object keys must be strings; "
                    f"got {type(key).__name__}"
                )
            if key in result:
                raise PhosPyInputError(
                    f"{field_name} contains duplicate JSON object key {key!r}"
                )
            result[key] = item
        return result
    raise PhosPyInputError(f"{field_name} must be an object")


def _require_sequence(value: object, *, field_name: str) -> list[object]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(value)
    raise PhosPyInputError(f"{field_name} must be an array")


def _require_payload_str(value: object, *, field_name: str) -> str:
    return _required_text(value, field_name=field_name)


def _optional_payload_str(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_payload_str(value, field_name=field_name)


def _require_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(f"{field_name} must be an int")
    return int(value)


__all__ = [
    "DERIVED_QUANTITATIVE_DATA_FINGERPRINT_ALGORITHM",
    "DERIVED_QUANTITATIVE_DATA_PROVENANCE_SCHEMA_VERSION_V1",
    "TECHNICAL_REPLICATE_AGGREGATION_DERIVATION_TYPE",
    "TECHNICAL_REPLICATE_AGGREGATOR_IMPLEMENTATION",
    "DerivedQuantitativeDataProvenance",
    "DerivedSampleMapping",
    "build_derived_quantitative_run_provenance",
]
