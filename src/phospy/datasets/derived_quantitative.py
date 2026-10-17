"""Derived quantitative lineage and Dataset model records."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import cast

import pandas as pd

from phospy.datasets._construction.service import (
    _AnalysisReadyDatasetConstructionService,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.datasets.processing_state import DatasetProcessingState
from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import DatasetValidationError
from phospy.provenance.hashing import (
    DEFAULT_STABLE_JSON_HASH_ALGORITHM,
    fingerprint_optional_table_strict,
    hash_json_payload,
)
from phospy.provenance.immutability import (
    freeze_json_mapping,
    freeze_json_value,
    thaw_json_mapping,
    thaw_json_value,
)
from phospy.provenance.models import (
    JsonValue,
    TableFingerprint,
)
from phospy.provenance.serialization import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)
from phospy.quantitative.models import IntensityScaleState
from phospy.run_lineage.models import RunProvenance
from phospy.taxonomy.organisms import Organism

DERIVED_QUANTITATIVE_DATA_PROVENANCE_SCHEMA_VERSION_V1 = 1
DERIVED_QUANTITATIVE_DATA_FINGERPRINT_ALGORITHM = DEFAULT_STABLE_JSON_HASH_ALGORITHM
TECHNICAL_REPLICATE_AGGREGATION_DERIVATION_TYPE = "technical_replicate_aggregation"
TECHNICAL_REPLICATE_AGGREGATOR_IMPLEMENTATION = (
    "phospy.differential.planning.replicates.TechnicalReplicateAggregator"
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
    "CertifiedDerivedQuantitativeParentState",
    "DerivedAnalysisReadyPhosphoDataset",
    "DerivedQuantitativeDataProvenance",
    "DerivedSampleMapping",
]
_CERTIFIED_PARENT_STATE_AUTHORITY = object()
_MATRIX_TRANSFORMATION_TABLES: Mapping[str, str] = {
    "phospho": "dataset.phospho",
    "total_protein": "dataset.total",
    "sample_metadata": "dataset.sample_metadata",
    "imputation_observation_mask": "dataset.imputation_observation_mask",
    "comparisons": "dataset.comparisons",
}
_TECHNICAL_REPLICATE_AGGREGATION = "technical_replicate_aggregation"
_TECHNICAL_REPLICATE_AGGREGATION_METHODS = frozenset(("mean", "median"))


@dataclass(frozen=True, slots=True, init=False)
class CertifiedDerivedQuantitativeParentState:
    """Certified parent state used to validate derived quantitative lineage."""

    parent_dataset_type: str
    parent_dataset_fingerprints: tuple[TableFingerprint, ...]
    phospho_sample_ids: tuple[str, ...]
    intensity_scale: str
    quantitative_meaning: str

    def __init__(
        self,
        *,
        parent_dataset_type: str,
        parent_dataset_fingerprints: tuple[TableFingerprint, ...],
        phospho_sample_ids: tuple[str, ...],
        intensity_scale: str,
        quantitative_meaning: str,
        _authority: object,
    ) -> None:
        if _authority is not _CERTIFIED_PARENT_STATE_AUTHORITY:
            raise PhosPyInputError(
                "CertifiedDerivedQuantitativeParentState must be created from an "
                "actual parent dataset or a certified internal snapshot"
            )
        object.__setattr__(
            self,
            "parent_dataset_type",
            _dataset_required_text(
                parent_dataset_type,
                field_name="derived_quantitative_parent_state.parent_dataset_type",
            ),
        )
        object.__setattr__(
            self,
            "parent_dataset_fingerprints",
            _require_table_fingerprint_tuple(
                parent_dataset_fingerprints,
                field_name=(
                    "derived_quantitative_parent_state.parent_dataset_fingerprints"
                ),
            ),
        )
        object.__setattr__(
            self,
            "phospho_sample_ids",
            _dataset_required_text_tuple(
                phospho_sample_ids,
                field_name="derived_quantitative_parent_state.phospho_sample_ids",
            ),
        )
        object.__setattr__(
            self,
            "intensity_scale",
            _dataset_required_text(
                intensity_scale,
                field_name="derived_quantitative_parent_state.intensity_scale",
            ),
        )
        object.__setattr__(
            self,
            "quantitative_meaning",
            _dataset_required_text(
                quantitative_meaning,
                field_name="derived_quantitative_parent_state.quantitative_meaning",
            ),
        )

    @classmethod
    def from_dataset(
        cls,
        dataset: AnalysisReadyPhosphoDataset,
    ) -> CertifiedDerivedQuantitativeParentState:
        """Capture certified state from an already-validated parent dataset."""

        if not isinstance(dataset, AnalysisReadyPhosphoDataset):
            raise PhosPyInputError(
                "derived quantitative parent_state must be captured from an "
                "AnalysisReadyPhosphoDataset"
            )
        phospho = dataset.phospho
        return cls._from_owned_parent_tables(
            parent_dataset_type=type(dataset).__name__,
            phospho=phospho,
            site_metadata=dataset.site_metadata,
            sample_metadata=dataset.sample_metadata,
            total=dataset.total,
            comparisons=dataset.comparisons,
            imputation_observation_mask=dataset.imputation_observed_mask_dataframe(),
            intensity_scale_state=dataset.intensity_scale_state,
        )

    @classmethod
    def _from_owned_parent_tables(
        cls,
        *,
        parent_dataset_type: str,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
    ) -> CertifiedDerivedQuantitativeParentState:
        """Package-private certified snapshot path for internally owned tables."""

        return cls(
            parent_dataset_type=parent_dataset_type,
            parent_dataset_fingerprints=_fingerprints_for_table_state(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                comparisons=comparisons,
                imputation_observation_mask=imputation_observation_mask,
            ),
            phospho_sample_ids=tuple(str(label) for label in phospho.columns.tolist()),
            intensity_scale=str(intensity_scale_state.label),
            quantitative_meaning=_quantitative_meaning_for_intensity_state(
                intensity_scale_state
            ),
            _authority=_CERTIFIED_PARENT_STATE_AUTHORITY,
        )


class DerivedAnalysisReadyPhosphoDataset(AnalysisReadyPhosphoDataset):
    """Analysis-ready dataset whose quantitative matrices were derived from a parent."""

    __slots__ = ("derived_lineage", "derived_parent_state")

    def __init__(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        parent_state: CertifiedDerivedQuantitativeParentState,
        derived_lineage: DerivedQuantitativeDataProvenance,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        allow_opaque_site_values: bool = False,
    ) -> None:
        self._init_from_derived_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            parent_state=parent_state,
            derived_lineage=derived_lineage,
            provenance=provenance,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=False,
        )

    def _init_from_derived_tables(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        parent_state: CertifiedDerivedQuantitativeParentState,
        derived_lineage: DerivedQuantitativeDataProvenance,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        allow_opaque_site_values: bool = False,
        _assume_owned: bool,
    ) -> None:
        if not isinstance(derived_lineage, DerivedQuantitativeDataProvenance):
            raise PhosPyInputError(
                "derived dataset requires DerivedQuantitativeDataProvenance"
            )
        if not isinstance(parent_state, CertifiedDerivedQuantitativeParentState):
            raise PhosPyInputError(
                "derived dataset requires CertifiedDerivedQuantitativeParentState"
            )
        if not isinstance(provenance, RunProvenance):
            raise PhosPyInputError("derived dataset requires RunProvenance")
        validated_tables = (
            _AnalysisReadyDatasetConstructionService().from_provenanced_tables(
                phospho=phospho,
                site_metadata=site_metadata,
                intensity_scale_state=intensity_scale_state,
                processing_state=processing_state,
                sample_metadata=sample_metadata,
                total=total,
                comparisons=comparisons,
                imputation_observation_mask=imputation_observation_mask,
                organism=organism,
                preprocessing_report=None,
                protein_aware_preparation=None,
                provenance=provenance,
                allow_opaque_site_values=allow_opaque_site_values,
                assume_owned=_assume_owned,
            )
        )
        self._init_from_validated_tables(validated_tables)
        _require_lineage_matches_owned_state(
            dataset=self,
            parent_state=parent_state,
            derived_lineage=derived_lineage,
            provenance=provenance,
        )
        object.__setattr__(self, "derived_lineage", derived_lineage)
        object.__setattr__(self, "derived_parent_state", parent_state)

    @classmethod
    def _from_owned_derived_tables(
        cls,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        intensity_scale_state: IntensityScaleState,
        processing_state: DatasetProcessingState,
        parent_state: CertifiedDerivedQuantitativeParentState,
        derived_lineage: DerivedQuantitativeDataProvenance,
        provenance: RunProvenance,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        comparisons: pd.DataFrame | None = None,
        imputation_observation_mask: pd.DataFrame | None = None,
        organism: Organism | None = None,
        allow_opaque_site_values: bool = False,
    ) -> DerivedAnalysisReadyPhosphoDataset:
        """Package-private owned transfer for internally produced derived tables."""

        dataset = cls.__new__(cls)
        dataset._init_from_derived_tables(
            phospho=phospho,
            site_metadata=site_metadata,
            intensity_scale_state=intensity_scale_state,
            processing_state=processing_state,
            parent_state=parent_state,
            sample_metadata=sample_metadata,
            total=total,
            comparisons=comparisons,
            imputation_observation_mask=imputation_observation_mask,
            organism=organism,
            provenance=provenance,
            derived_lineage=derived_lineage,
            allow_opaque_site_values=allow_opaque_site_values,
            _assume_owned=True,
        )
        return dataset


def _require_lineage_matches_owned_state(
    *,
    dataset: DerivedAnalysisReadyPhosphoDataset,
    parent_state: CertifiedDerivedQuantitativeParentState,
    derived_lineage: DerivedQuantitativeDataProvenance,
    provenance: RunProvenance,
) -> None:
    actual_fingerprints = _fingerprints_for_owned_derived_dataset(dataset)
    _require_sample_mapping_matches_actual_phospho(
        dataset=dataset,
        derived_lineage=derived_lineage,
    )
    _require_parent_identity_matches_certified_state(
        parent_state=parent_state,
        derived_lineage=derived_lineage,
    )
    _require_sample_mapping_matches_certified_parent(
        parent_state=parent_state,
        derived_lineage=derived_lineage,
    )
    _require_quantitative_state_matches_certified_boundaries(
        parent_state=parent_state,
        dataset=dataset,
        derived_lineage=derived_lineage,
    )
    _require_fingerprint_sets_match(
        expected=derived_lineage.derived_dataset_fingerprints,
        actual=actual_fingerprints,
        field_name="derived_quantitative_data.derived_dataset_fingerprints",
        expected_source="actual derived dataset tables",
    )
    _require_matrices_transformed_match_table_state(
        parent_state=parent_state,
        actual_derived_fingerprints=actual_fingerprints,
        derived_lineage=derived_lineage,
    )
    _require_fingerprint_sets_match(
        expected=derived_lineage.parent_dataset_fingerprints,
        actual=provenance.input_tables,
        field_name="run_provenance.input_tables",
        expected_source="derived_lineage.parent_dataset_fingerprints",
    )
    _require_fingerprint_sets_match(
        expected=derived_lineage.derived_dataset_fingerprints,
        actual=provenance.output_tables,
        field_name="run_provenance.output_tables",
        expected_source="derived_lineage.derived_dataset_fingerprints",
    )
    if provenance.workflow_name != derived_lineage.derivation_type:
        raise DatasetValidationError(
            "derived dataset run_provenance.workflow_name must match "
            "derived_lineage.derivation_type"
        )
    _require_run_provenance_embeds_lineage(
        provenance=provenance,
        derived_lineage=derived_lineage,
    )


def _fingerprints_for_owned_derived_dataset(
    dataset: DerivedAnalysisReadyPhosphoDataset,
) -> tuple[TableFingerprint, ...]:
    return _fingerprints_for_table_state(
        phospho=dataset._phospho,
        site_metadata=dataset._site_metadata,
        sample_metadata=dataset._sample_metadata,
        total=dataset._total,
        comparisons=dataset._comparisons,
        imputation_observation_mask=_owned_imputation_observation_mask(dataset),
    )


def _fingerprints_for_table_state(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    sample_metadata: pd.DataFrame | None,
    total: pd.DataFrame | None,
    comparisons: pd.DataFrame | None,
    imputation_observation_mask: pd.DataFrame | None,
) -> tuple[TableFingerprint, ...]:
    entries: tuple[tuple[str, pd.DataFrame | None], ...] = (
        ("dataset.phospho", phospho),
        ("dataset.site_metadata", site_metadata),
        ("dataset.sample_metadata", sample_metadata),
        ("dataset.total", total),
        ("dataset.comparisons", comparisons),
        ("dataset.imputation_observation_mask", imputation_observation_mask),
    )
    return _collect_table_fingerprints(entries)


def _owned_imputation_observation_mask(
    dataset: DerivedAnalysisReadyPhosphoDataset,
) -> pd.DataFrame | None:
    metadata = dataset._imputation_observation_metadata
    if metadata is None:
        return None
    return metadata.observed_mask_dataframe()


def _collect_table_fingerprints(
    entries: Sequence[tuple[str, pd.DataFrame | None]],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table_strict(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)


def _require_sample_mapping_matches_actual_phospho(
    *,
    dataset: DerivedAnalysisReadyPhosphoDataset,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    expected_columns = tuple(
        item.output_sample_id for item in derived_lineage.sample_mapping
    )
    actual_columns = tuple(str(label) for label in dataset._phospho.columns.tolist())
    if expected_columns != actual_columns:
        raise DatasetValidationError(
            "derived_quantitative_data.sample_mapping output_sample_id values "
            "must match actual derived dataset.phospho columns"
        )


def _require_parent_identity_matches_certified_state(
    *,
    parent_state: CertifiedDerivedQuantitativeParentState,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    if derived_lineage.parent_dataset_type != parent_state.parent_dataset_type:
        raise DatasetValidationError(
            "derived_quantitative_data.parent_dataset_type must match certified "
            "parent dataset state"
        )
    _require_fingerprint_sets_match(
        expected=derived_lineage.parent_dataset_fingerprints,
        actual=parent_state.parent_dataset_fingerprints,
        field_name="derived_quantitative_data.parent_dataset_fingerprints",
        expected_source="certified parent dataset state",
    )


def _require_sample_mapping_matches_certified_parent(
    *,
    parent_state: CertifiedDerivedQuantitativeParentState,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    parent_sample_ids = parent_state.phospho_sample_ids
    declared_input_sample_ids = tuple(
        sample_id
        for item in derived_lineage.sample_mapping
        for sample_id in item.input_sample_ids
    )
    missing_from_parent = sorted(
        set(declared_input_sample_ids) - set(parent_sample_ids)
    )
    if missing_from_parent:
        raise DatasetValidationError(
            "derived_quantitative_data.sample_mapping input_sample_ids must exist "
            "in certified parent dataset.phospho columns; missing: "
            + ", ".join(missing_from_parent)
        )
    duplicate_input_ids = sorted(
        sample_id
        for sample_id, count in Counter(declared_input_sample_ids).items()
        if count > 1
    )
    if duplicate_input_ids:
        raise DatasetValidationError(
            "derived_quantitative_data.sample_mapping input_sample_ids must not "
            "duplicate parent dataset.phospho columns for this operation; "
            "duplicates: " + ", ".join(duplicate_input_ids)
        )
    omitted_parent_ids = [
        sample_id
        for sample_id in parent_sample_ids
        if sample_id not in set(declared_input_sample_ids)
    ]
    if omitted_parent_ids:
        raise DatasetValidationError(
            "derived_quantitative_data.sample_mapping input_sample_ids must cover "
            "every certified parent dataset.phospho column for this operation; "
            "omitted: " + ", ".join(omitted_parent_ids)
        )
    _require_operation_specific_sample_mapping(derived_lineage=derived_lineage)


def _require_operation_specific_sample_mapping(
    *,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    if derived_lineage.derivation_type != _TECHNICAL_REPLICATE_AGGREGATION:
        return
    if (
        derived_lineage.aggregation_method
        not in _TECHNICAL_REPLICATE_AGGREGATION_METHODS
    ):
        raise DatasetValidationError(
            "technical_replicate_aggregation derived_quantitative_data."
            "aggregation_method must be 'mean' or 'median'"
        )
    group_keys = tuple(
        (item.condition, item.biological_replicate_id)
        for item in derived_lineage.sample_mapping
    )
    duplicate_group_keys = sorted(
        f"{condition}:{biological_replicate_id}"
        for (condition, biological_replicate_id), count in Counter(group_keys).items()
        if count > 1
    )
    if duplicate_group_keys:
        raise DatasetValidationError(
            "technical_replicate_aggregation sample_mapping must have one output "
            "per condition+biological_replicate_id group; duplicate groups: "
            + ", ".join(duplicate_group_keys)
        )
    biological_replicate_counts = Counter(
        item.biological_replicate_id for item in derived_lineage.sample_mapping
    )
    for item in derived_lineage.sample_mapping:
        expected_output_sample_id = item.biological_replicate_id
        if biological_replicate_counts[item.biological_replicate_id] > 1:
            expected_output_sample_id = (
                f"{item.condition}__{item.biological_replicate_id}"
            )
        if item.output_sample_id != expected_output_sample_id:
            raise DatasetValidationError(
                "technical_replicate_aggregation sample_mapping output_sample_id "
                "must match the condition+biological_replicate_id grouping rule; "
                f"expected {expected_output_sample_id!r} for "
                f"condition={item.condition!r}, "
                f"biological_replicate_id={item.biological_replicate_id!r}"
            )
        duplicate_technical_ids = sorted(
            technical_id
            for technical_id, count in Counter(item.technical_replicate_ids).items()
            if count > 1
        )
        if duplicate_technical_ids:
            raise DatasetValidationError(
                "technical_replicate_aggregation sample_mapping "
                "technical_replicate_ids must be unique within each output "
                "sample; duplicates: " + ", ".join(duplicate_technical_ids)
            )
        if len(item.technical_replicate_ids) > len(item.input_sample_ids):
            raise DatasetValidationError(
                "technical_replicate_aggregation sample_mapping "
                "technical_replicate_ids cannot outnumber input_sample_ids"
            )


def _require_quantitative_state_matches_certified_boundaries(
    *,
    parent_state: CertifiedDerivedQuantitativeParentState,
    dataset: DerivedAnalysisReadyPhosphoDataset,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    output_intensity_scale = str(dataset.intensity_scale_state.label)
    output_quantitative_meaning = _quantitative_meaning_for_intensity_state(
        dataset.intensity_scale_state
    )
    if derived_lineage.input_intensity_scale != parent_state.intensity_scale:
        raise DatasetValidationError(
            "derived_quantitative_data.input_intensity_scale must match certified "
            "parent dataset intensity scale"
        )
    if derived_lineage.output_intensity_scale != output_intensity_scale:
        raise DatasetValidationError(
            "derived_quantitative_data.output_intensity_scale must match actual "
            "derived dataset intensity scale"
        )
    if derived_lineage.quantitative_meaning != parent_state.quantitative_meaning:
        raise DatasetValidationError(
            "derived_quantitative_data.quantitative_meaning must match certified "
            "parent dataset quantitative meaning"
        )
    if derived_lineage.quantitative_meaning != output_quantitative_meaning:
        raise DatasetValidationError(
            "derived_quantitative_data.quantitative_meaning must match actual "
            "derived dataset quantitative meaning"
        )


def _require_matrices_transformed_match_table_state(
    *,
    parent_state: CertifiedDerivedQuantitativeParentState,
    actual_derived_fingerprints: tuple[TableFingerprint, ...],
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    expected = _actual_matrices_transformed(
        parent_fingerprints=parent_state.parent_dataset_fingerprints,
        derived_fingerprints=actual_derived_fingerprints,
    )
    declared = dict(derived_lineage.matrices_transformed)
    missing_keys = sorted(set(expected) - set(declared))
    unexpected_keys = sorted(set(declared) - set(expected))
    if missing_keys or unexpected_keys:
        detail_parts: list[str] = []
        if missing_keys:
            detail_parts.append("missing keys: " + ", ".join(missing_keys))
        if unexpected_keys:
            detail_parts.append("unexpected keys: " + ", ".join(unexpected_keys))
        raise DatasetValidationError(
            "derived_quantitative_data.matrices_transformed keys must match "
            "certified parent and actual derived table state; "
            + "; ".join(detail_parts)
        )
    mismatched = sorted(
        key
        for key, expected_value in expected.items()
        if declared[key] != expected_value
    )
    if mismatched:
        details = ", ".join(f"{key} expected {expected[key]!r}" for key in mismatched)
        raise DatasetValidationError(
            "derived_quantitative_data.matrices_transformed must match certified "
            "parent and actual derived table state; mismatched flags: " + details
        )


def _actual_matrices_transformed(
    *,
    parent_fingerprints: tuple[TableFingerprint, ...],
    derived_fingerprints: tuple[TableFingerprint, ...],
) -> dict[str, bool]:
    return {
        key: (
            _dataset_fingerprint_by_name(parent_fingerprints, table_name)
            != _dataset_fingerprint_by_name(derived_fingerprints, table_name)
        )
        for key, table_name in _MATRIX_TRANSFORMATION_TABLES.items()
    }


def _dataset_fingerprint_by_name(
    fingerprints: tuple[TableFingerprint, ...],
    name: str,
) -> TableFingerprint | None:
    for fingerprint in fingerprints:
        if fingerprint.name == name:
            return fingerprint
    return None


def _require_fingerprint_sets_match(
    *,
    expected: tuple[TableFingerprint, ...],
    actual: tuple[TableFingerprint, ...],
    field_name: str,
    expected_source: str,
) -> None:
    expected_by_name = _fingerprint_map(expected, field_name=field_name)
    actual_by_name = _fingerprint_map(actual, field_name=field_name)
    missing = sorted(set(expected_by_name) - set(actual_by_name))
    unexpected = sorted(set(actual_by_name) - set(expected_by_name))
    if missing or unexpected:
        detail_parts: list[str] = []
        if missing:
            detail_parts.append("missing actual tables: " + ", ".join(missing))
        if unexpected:
            detail_parts.append("unexpected actual tables: " + ", ".join(unexpected))
        raise DatasetValidationError(
            f"{field_name} table set must match {expected_source}; "
            + "; ".join(detail_parts)
        )
    for name in expected_by_name:
        _require_fingerprint_matches(
            expected=expected_by_name[name],
            actual=actual_by_name[name],
            field_name=f"{field_name}.{name}",
            expected_source=expected_source,
        )


def _fingerprint_map(
    fingerprints: tuple[TableFingerprint, ...],
    *,
    field_name: str,
) -> dict[str, TableFingerprint]:
    result: dict[str, TableFingerprint] = {}
    for fingerprint in fingerprints:
        if fingerprint.name in result:
            raise DatasetValidationError(
                f"{field_name} contains duplicate table fingerprint "
                f"{fingerprint.name!r}"
            )
        result[fingerprint.name] = fingerprint
    return result


def _require_fingerprint_matches(
    *,
    expected: TableFingerprint,
    actual: TableFingerprint,
    field_name: str,
    expected_source: str,
) -> None:
    checks: tuple[tuple[str, object, object], ...] = (
        ("rows", expected.rows, actual.rows),
        ("columns", expected.columns, actual.columns),
        ("index_name", expected.index_name, actual.index_name),
        ("column_names", expected.column_names, actual.column_names),
        ("dtypes", expected.dtypes, actual.dtypes),
        ("index_structure", expected.index_structure, actual.index_structure),
        (
            "column_index_structure",
            expected.column_index_structure,
            actual.column_index_structure,
        ),
        (
            "exact_hash_algorithm",
            expected.exact_hash_algorithm,
            actual.exact_hash_algorithm,
        ),
        ("exact_hash_value", expected.exact_hash_value, actual.exact_hash_value),
        (
            "tolerance_hash_algorithm",
            expected.tolerance_hash_algorithm,
            actual.tolerance_hash_algorithm,
        ),
        (
            "tolerance_hash_value",
            expected.tolerance_hash_value,
            actual.tolerance_hash_value,
        ),
    )
    mismatched = [
        name
        for name, expected_value, actual_value in checks
        if expected_value != actual_value
    ]
    if mismatched:
        raise DatasetValidationError(
            f"{field_name} does not match {expected_source}; "
            "mismatched fields: " + ", ".join(mismatched)
        )


def _require_run_provenance_embeds_lineage(
    *,
    provenance: RunProvenance,
    derived_lineage: DerivedQuantitativeDataProvenance,
) -> None:
    payload = provenance.workflow_parameters.get("derived_quantitative_data")
    if not isinstance(payload, Mapping):
        raise DatasetValidationError(
            "derived dataset run_provenance.workflow_parameters must include "
            "typed derived_quantitative_data lineage"
        )
    try:
        restored = DerivedQuantitativeDataProvenance.from_payload(
            _thaw_mapping(payload)
        )
    except PhosPyInputError as exc:
        raise DatasetValidationError(
            "derived dataset run_provenance.workflow_parameters "
            "derived_quantitative_data must decode as typed lineage"
        ) from exc
    if restored.to_payload() != derived_lineage.to_payload():
        raise DatasetValidationError(
            "derived dataset run_provenance.workflow_parameters "
            "derived_quantitative_data must match derived_lineage"
        )


def _thaw_mapping(payload: Mapping[str, object]) -> Mapping[str, object]:
    thawed = thaw_json_value(payload, field_name="derived_dataset.lineage")
    if not isinstance(thawed, Mapping):
        raise DatasetValidationError(
            "derived dataset lineage payload must thaw to a JSON object"
        )
    return thawed


def _quantitative_meaning_for_intensity_state(
    intensity_scale_state: IntensityScaleState,
) -> str:
    quantity = intensity_scale_state.quantity
    if quantity is None:
        return "unknown"
    return quantity.value


def _require_table_fingerprint_tuple(
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


def _dataset_required_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    normalized = value.strip()
    if not normalized:
        raise PhosPyInputError(f"{field_name} must be a non-empty string")
    return normalized


def _dataset_required_text_tuple(values: object, *, field_name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    result = tuple(
        _dataset_required_text(value, field_name=f"{field_name}[]") for value in values
    )
    if not result:
        raise PhosPyInputError(f"{field_name} must not be empty")
    return result
