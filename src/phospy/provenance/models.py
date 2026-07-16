"""Typed machine-readable provenance models."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from numbers import Integral
from typing import Protocol, TypeAlias

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.immutability import (
    FrozenJsonValue,
    freeze_json_mapping,
    freeze_optional_json_mapping,
    thaw_json_mapping,
)
from phospy.provenance.reference_identifiers import (
    ReferenceIdentifierNormalisationReport,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V1 = 1
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2 = 2
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3 = 3
ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V1 = 1
ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2 = 2
BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1 = 1
TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V1 = 1
TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V2 = 2
_TRUSTED_DATASET_CONSTRUCTION_REQUIRED_DIMENSIONS = (
    "identity",
    "quantitative_meaning",
    "localisation",
    "sequence",
    "reference_context",
)
_TRUSTED_DATASET_CONSTRUCTION_MISSING_ASSERTION_NAMES = {
    "identity": "identity_user_asserted",
    "quantitative_meaning": "quantitative_meaning_user_asserted",
    "localisation": "localisation_user_asserted",
    "sequence": "sequence_user_asserted",
    "reference_context": "reference_context_user_asserted",
}
_TRUSTED_DATASET_CONSTRUCTION_EVIDENCE_KINDS = frozenset({"evidence", "waiver"})


class DeterminismKind(str, Enum):
    """Declared execution determinism for preprocessing stage provenance."""

    DETERMINISTIC = "deterministic"
    SEEDED_STOCHASTIC = "seeded_stochastic"
    EXTERNALLY_NONDETERMINISTIC = "externally_nondeterministic"


PREPROCESSING_STAGE_DETERMINISM_DETERMINISTIC = DeterminismKind.DETERMINISTIC.value
PREPROCESSING_STAGE_DETERMINISM_PURE = PREPROCESSING_STAGE_DETERMINISM_DETERMINISTIC
PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC = (
    DeterminismKind.SEEDED_STOCHASTIC.value
)
PREPROCESSING_STAGE_DETERMINISM_EXTERNALLY_NONDETERMINISTIC = (
    DeterminismKind.EXTERNALLY_NONDETERMINISTIC.value
)
PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY = (
    PREPROCESSING_STAGE_DETERMINISM_EXTERNALLY_NONDETERMINISTIC
)
PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE = (
    "preprocessing_external_nondeterminism"
)
BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS = frozenset(
    {"not_provided", "unknown", "none", "null"}
)
_REPRODUCIBILITY_CAVEAT_SEVERITIES = frozenset({"warning", "error"})

JsonPrimitive = str | int | float | bool | None
JsonValue: TypeAlias = (
    FrozenJsonValue
    | tuple["JsonValue", ...]
    | list["JsonValue"]
    | dict[str, "JsonValue"]
)


class ReferenceContextProtocol(Protocol):
    """Structural protocol for reference-context provenance values."""

    @property
    def organism(self) -> str: ...

    @property
    def protein_namespace(self) -> str: ...

    @property
    def source_name(self) -> str: ...

    @property
    def source_version(self) -> str: ...

    @property
    def proteome_version(self) -> str | None: ...

    @property
    def reference_table_sha256(self) -> str | None: ...

    @property
    def reference_context_id(self) -> str: ...

    def to_payload(self) -> Mapping[str, JsonValue]: ...


@dataclass(frozen=True, slots=True)
class ReproducibilityCaveat:
    """Machine-readable reproducibility caveat for provenance records."""

    code: str
    severity: str
    message: str
    details: Mapping[str, JsonValue] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "code",
            _required_provenance_text(
                self.code, field_name="reproducibility_caveat.code"
            ),
        )
        severity = _required_provenance_text(
            self.severity,
            field_name="reproducibility_caveat.severity",
        )
        if severity not in _REPRODUCIBILITY_CAVEAT_SEVERITIES:
            supported = ", ".join(sorted(_REPRODUCIBILITY_CAVEAT_SEVERITIES))
            raise PhosPyInputError(
                "reproducibility_caveat.severity must be one of: " + supported
            )
        object.__setattr__(self, "severity", severity)
        object.__setattr__(
            self,
            "message",
            _required_provenance_text(
                self.message,
                field_name="reproducibility_caveat.message",
            ),
        )
        object.__setattr__(
            self,
            "details",
            freeze_json_mapping(
                self.details,
                field_name="reproducibility_caveat.details",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible caveat payload."""

        return {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "details": thaw_json_mapping(
                self.details,
                field_name="reproducibility_caveat.details",
            ),
        }


def _empty_platform_provenance() -> dict[str, str]:
    return {}


def _empty_json_mapping() -> dict[str, JsonValue]:
    return {}


def _required_provenance_text(value: object, *, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        raise PhosPyInputError(f"{field_name} must be non-empty")
    return text


def _required_provenance_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise PhosPyInputError(f"{field_name} must be a bool")
    return bool(value)


def _optional_provenance_text(value: object | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _optional_provenance_float(
    value: object | None, *, field_name: str
) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhosPyInputError(f"{field_name} must be a finite number")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise PhosPyInputError(f"{field_name} must be a finite number")
    return numeric


@dataclass(frozen=True, slots=True)
class InputIntensityScaleEvidence:
    """Workflow-visible evidence for how input intensity scale was established."""

    input_intensity_scale: str
    input_intensity_scale_evidence_level: str
    input_intensity_scale_source: str
    input_intensity_scale_source_detail: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "input_intensity_scale",
            _required_provenance_text(
                self.input_intensity_scale,
                field_name="input_intensity_scale",
            ),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_evidence_level",
            _required_provenance_text(
                self.input_intensity_scale_evidence_level,
                field_name="input_intensity_scale_evidence_level",
            ),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_source",
            _required_provenance_text(
                self.input_intensity_scale_source,
                field_name="input_intensity_scale_source",
            ),
        )
        object.__setattr__(
            self,
            "input_intensity_scale_source_detail",
            _optional_provenance_text(self.input_intensity_scale_source_detail),
        )

    def to_payload(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "input_intensity_scale": self.input_intensity_scale,
            "input_intensity_scale_evidence_level": (
                self.input_intensity_scale_evidence_level
            ),
            "input_intensity_scale_source": self.input_intensity_scale_source,
        }
        if self.input_intensity_scale_source_detail is not None:
            payload["input_intensity_scale_source_detail"] = (
                self.input_intensity_scale_source_detail
            )
        return payload


@dataclass(frozen=True, slots=True)
class TrustedDatasetConstructionEvidence:
    """Typed evidence or waiver for one trusted construction assertion."""

    kind: str
    source: str | None = None
    policy: str | None = None
    threshold: float | None = None
    waiver_reason: str | None = None
    details: Mapping[str, JsonValue] = field(default_factory=dict)

    def __post_init__(self) -> None:
        kind = _required_provenance_text(
            self.kind,
            field_name="trusted_dataset_construction_evidence.kind",
        )
        if kind not in _TRUSTED_DATASET_CONSTRUCTION_EVIDENCE_KINDS:
            supported = ", ".join(sorted(_TRUSTED_DATASET_CONSTRUCTION_EVIDENCE_KINDS))
            raise PhosPyInputError(
                "trusted_dataset_construction_evidence.kind must be one of: "
                + supported
            )
        source = _optional_provenance_text(self.source)
        policy = _optional_provenance_text(self.policy)
        waiver_reason = _optional_provenance_text(self.waiver_reason)
        threshold = _optional_provenance_float(
            self.threshold,
            field_name="trusted_dataset_construction_evidence.threshold",
        )
        if kind == "evidence" and source is None:
            raise PhosPyInputError(
                "trusted_dataset_construction_evidence.source is required when "
                "kind='evidence'"
            )
        if kind == "waiver" and waiver_reason is None:
            raise PhosPyInputError(
                "trusted_dataset_construction_evidence.waiver_reason is required "
                "when kind='waiver'"
            )
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "policy", policy)
        object.__setattr__(self, "threshold", threshold)
        object.__setattr__(self, "waiver_reason", waiver_reason)
        object.__setattr__(
            self,
            "details",
            freeze_json_mapping(
                self.details,
                field_name="trusted_dataset_construction_evidence.details",
            ),
        )

    @classmethod
    def evidence(
        cls,
        *,
        source: str,
        policy: str | None = None,
        threshold: float | None = None,
        details: Mapping[str, JsonValue] | None = None,
    ) -> TrustedDatasetConstructionEvidence:
        """Create an explicit evidence record."""

        return cls(
            kind="evidence",
            source=source,
            policy=policy,
            threshold=threshold,
            details={} if details is None else details,
        )

    @classmethod
    def waiver(
        cls,
        *,
        reason: str,
        policy: str | None = None,
        details: Mapping[str, JsonValue] | None = None,
    ) -> TrustedDatasetConstructionEvidence:
        """Create an explicit waiver record."""

        return cls(
            kind="waiver",
            policy=policy,
            waiver_reason=reason,
            details={} if details is None else details,
        )

    @property
    def is_waiver(self) -> bool:
        """Return whether this assertion is an explicit waiver."""

        return self.kind == "waiver"

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible evidence payload."""

        return {
            "kind": self.kind,
            "source": self.source,
            "policy": self.policy,
            "threshold": self.threshold,
            "waiver_reason": self.waiver_reason,
            "details": thaw_json_mapping(
                self.details,
                field_name="trusted_dataset_construction_evidence.details",
            ),
        }


def _optional_trusted_construction_evidence(
    value: object | None,
    *,
    field_name: str,
) -> TrustedDatasetConstructionEvidence | None:
    if value is None:
        return None
    if not isinstance(value, TrustedDatasetConstructionEvidence):
        raise PhosPyInputError(
            f"{field_name} must be TrustedDatasetConstructionEvidence or None"
        )
    return value


@dataclass(frozen=True, slots=True)
class TrustedDatasetConstructionAssertions:
    """User assertion provenance for trusted direct dataset construction.

    Complete trusted construction metadata records typed evidence or an
    explicit waiver for identity, quantitative meaning, localisation, sequence,
    and reference context. A missing assertion bundle is reserved for legacy
    direct-construction audit markers.
    """

    identity: TrustedDatasetConstructionEvidence | None = None
    quantitative_meaning: TrustedDatasetConstructionEvidence | None = None
    localisation: TrustedDatasetConstructionEvidence | None = None
    sequence: TrustedDatasetConstructionEvidence | None = None
    reference_context: TrustedDatasetConstructionEvidence | None = None
    assertion_metadata_provided: bool = True
    asserted_by: str | None = None
    assertion_source: str | None = None
    notes: str | None = None
    schema_version: int = TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V2

    def __post_init__(self) -> None:
        identity = _optional_trusted_construction_evidence(
            self.identity,
            field_name="trusted_dataset_construction_assertions.identity",
        )
        quantitative_meaning = _optional_trusted_construction_evidence(
            self.quantitative_meaning,
            field_name="trusted_dataset_construction_assertions.quantitative_meaning",
        )
        localisation = _optional_trusted_construction_evidence(
            self.localisation,
            field_name="trusted_dataset_construction_assertions.localisation",
        )
        sequence = _optional_trusted_construction_evidence(
            self.sequence,
            field_name="trusted_dataset_construction_assertions.sequence",
        )
        reference_context = _optional_trusted_construction_evidence(
            self.reference_context,
            field_name="trusted_dataset_construction_assertions.reference_context",
        )
        assertion_metadata_provided = _required_provenance_bool(
            self.assertion_metadata_provided,
            field_name=(
                "trusted_dataset_construction_assertions.assertion_metadata_provided"
            ),
        )
        supplied_assertions = {
            "identity": identity,
            "quantitative_meaning": quantitative_meaning,
            "localisation": localisation,
            "sequence": sequence,
            "reference_context": reference_context,
        }
        if not assertion_metadata_provided and any(
            value is not None for value in supplied_assertions.values()
        ):
            raise PhosPyInputError(
                "trusted_dataset_construction_assertions cannot record user "
                "assertions when assertion_metadata_provided is False"
            )
        if assertion_metadata_provided:
            missing = tuple(
                name for name, value in supplied_assertions.items() if value is None
            )
            if missing:
                raise PhosPyInputError(
                    "trusted_dataset_construction_assertions requires typed "
                    "evidence or an explicit waiver for: " + ", ".join(missing)
                )
            _require_localisation_evidence_or_waiver(localisation)
        schema_version = _required_non_negative_row_count(
            self.schema_version,
            field_name="trusted_dataset_construction_assertions.schema_version",
        )
        if schema_version != TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V2:
            raise PhosPyInputError(
                "trusted_dataset_construction_assertions.schema_version must be "
                f"{TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V2}"
            )
        object.__setattr__(self, "identity", identity)
        object.__setattr__(self, "quantitative_meaning", quantitative_meaning)
        object.__setattr__(self, "localisation", localisation)
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "reference_context", reference_context)
        object.__setattr__(
            self,
            "assertion_metadata_provided",
            assertion_metadata_provided,
        )
        object.__setattr__(
            self,
            "asserted_by",
            _optional_provenance_text(self.asserted_by),
        )
        object.__setattr__(
            self,
            "assertion_source",
            _optional_provenance_text(self.assertion_source),
        )
        object.__setattr__(self, "notes", _optional_provenance_text(self.notes))
        object.__setattr__(self, "schema_version", schema_version)

    @property
    def identity_user_asserted(self) -> bool:
        """Return whether identity has typed evidence or a waiver."""

        return self.identity is not None

    @property
    def quantitative_meaning_user_asserted(self) -> bool:
        """Return whether quantitative meaning has typed evidence or a waiver."""

        return self.quantitative_meaning is not None

    @property
    def localisation_user_asserted(self) -> bool:
        """Return whether localisation has typed evidence or a waiver."""

        return self.localisation is not None

    @property
    def sequence_user_asserted(self) -> bool:
        """Return whether sequence has typed evidence or a waiver."""

        return self.sequence is not None

    @property
    def reference_context_user_asserted(self) -> bool:
        """Return whether reference context has typed evidence or a waiver."""

        return self.reference_context is not None

    @property
    def assertion_fingerprint(self) -> str:
        """Return a stable fingerprint of the assertion payload."""

        encoded = json.dumps(
            self._to_payload(include_fingerprint=False),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    @property
    def waived_assertions(self) -> tuple[str, ...]:
        """Return assertion dimensions satisfied by explicit waiver."""

        waived: list[str] = []
        for dimension in _TRUSTED_DATASET_CONSTRUCTION_REQUIRED_DIMENSIONS:
            record = getattr(self, dimension)
            if isinstance(record, TrustedDatasetConstructionEvidence) and (
                record.is_waiver
            ):
                waived.append(dimension)
        return tuple(waived)

    @classmethod
    def missing(cls) -> TrustedDatasetConstructionAssertions:
        """Return explicit metadata for absent trusted assertion provenance."""

        return cls(
            assertion_metadata_provided=False,
            notes="No typed trusted construction assertion metadata was supplied.",
        )

    @property
    def missing_assertions(self) -> tuple[str, ...]:
        """Return required assertion fields not recorded or waived."""

        if not self.assertion_metadata_provided:
            return tuple(
                _TRUSTED_DATASET_CONSTRUCTION_MISSING_ASSERTION_NAMES[dimension]
                for dimension in _TRUSTED_DATASET_CONSTRUCTION_REQUIRED_DIMENSIONS
            )
        missing: list[str] = []
        for dimension in _TRUSTED_DATASET_CONSTRUCTION_REQUIRED_DIMENSIONS:
            if getattr(self, dimension) is None:
                missing.append(
                    _TRUSTED_DATASET_CONSTRUCTION_MISSING_ASSERTION_NAMES[dimension]
                )
        return tuple(missing)

    @property
    def all_required_assertions_present(self) -> bool:
        """Return whether all required trusted construction assertions are present."""

        return not self.missing_assertions

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible trusted assertion payload."""

        return self._to_payload(include_fingerprint=True)

    def _to_payload(self, *, include_fingerprint: bool) -> dict[str, object]:
        payload: dict[str, object] = {
            "schema_version": int(self.schema_version),
            "assertion_metadata_provided": bool(self.assertion_metadata_provided),
            "identity": _trusted_evidence_payload(self.identity),
            "quantitative_meaning": _trusted_evidence_payload(
                self.quantitative_meaning
            ),
            "localisation": _trusted_evidence_payload(self.localisation),
            "sequence": _trusted_evidence_payload(self.sequence),
            "reference_context": _trusted_evidence_payload(self.reference_context),
            "identity_user_asserted": bool(self.identity_user_asserted),
            "quantitative_meaning_user_asserted": bool(
                self.quantitative_meaning_user_asserted
            ),
            "localisation_user_asserted": bool(self.localisation_user_asserted),
            "sequence_user_asserted": bool(self.sequence_user_asserted),
            "reference_context_user_asserted": bool(
                self.reference_context_user_asserted
            ),
            "waived_assertions": list(self.waived_assertions),
            "missing_assertions": list(self.missing_assertions),
            "asserted_by": self.asserted_by,
            "assertion_source": self.assertion_source,
            "notes": self.notes,
        }
        if include_fingerprint:
            payload["assertion_fingerprint"] = self.assertion_fingerprint
        return payload


def _trusted_evidence_payload(
    value: TrustedDatasetConstructionEvidence | None,
) -> dict[str, object] | None:
    if value is None:
        return None
    return value.to_payload()


def _require_localisation_evidence_or_waiver(
    value: TrustedDatasetConstructionEvidence | None,
) -> None:
    if value is None:
        return
    if value.kind == "waiver":
        return
    if value.policy is None:
        raise PhosPyInputError(
            "trusted_dataset_construction_assertions.localisation.policy is "
            "required when localisation is recorded as evidence"
        )
    if value.threshold is None:
        raise PhosPyInputError(
            "trusted_dataset_construction_assertions.localisation.threshold is "
            "required when localisation is recorded as evidence"
        )
    if value.threshold < 0.0 or value.threshold > 1.0:
        raise PhosPyInputError(
            "trusted_dataset_construction_assertions.localisation.threshold must "
            "be between 0 and 1"
        )


@dataclass(frozen=True, slots=True)
class RowAttritionRecord:
    """Standard count-only provenance for rows removed at one workflow stage."""

    stage: str
    input_rows: int
    output_rows: int
    removed_rows: int
    reason: str
    examples: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        input_rows = _required_non_negative_row_count(
            self.input_rows,
            field_name="row_attrition_record.input_rows",
        )
        output_rows = _required_non_negative_row_count(
            self.output_rows,
            field_name="row_attrition_record.output_rows",
        )
        removed_rows = _required_non_negative_row_count(
            self.removed_rows,
            field_name="row_attrition_record.removed_rows",
        )
        if output_rows > input_rows:
            raise PhosPyInputError(
                "row_attrition_record.output_rows must be less than or equal "
                "to input_rows"
            )
        expected_removed_rows = input_rows - output_rows
        if removed_rows != expected_removed_rows:
            raise PhosPyInputError(
                "row_attrition_record.removed_rows must equal input_rows - output_rows"
            )
        object.__setattr__(
            self,
            "stage",
            _required_provenance_text(
                self.stage,
                field_name="row_attrition_record.stage",
            ),
        )
        object.__setattr__(self, "input_rows", input_rows)
        object.__setattr__(self, "output_rows", output_rows)
        object.__setattr__(self, "removed_rows", removed_rows)
        object.__setattr__(
            self,
            "reason",
            _required_provenance_text(
                self.reason,
                field_name="row_attrition_record.reason",
            ),
        )
        object.__setattr__(
            self,
            "examples",
            _required_provenance_text_tuple(
                self.examples,
                field_name="row_attrition_record.examples",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible row-attrition record payload."""

        return {
            "stage": self.stage,
            "input_rows": int(self.input_rows),
            "output_rows": int(self.output_rows),
            "removed_rows": int(self.removed_rows),
            "reason": self.reason,
            "examples": list(self.examples),
        }


@dataclass(frozen=True, slots=True)
class RowAttritionReport:
    """Ordered row-attrition provenance across workflow stages."""

    records: tuple[RowAttritionRecord, ...]
    input_rows: int
    final_rows: int

    def __post_init__(self) -> None:
        records = _required_row_attrition_record_tuple(self.records)
        input_rows = _required_non_negative_row_count(
            self.input_rows,
            field_name="row_attrition_report.input_rows",
        )
        final_rows = _required_non_negative_row_count(
            self.final_rows,
            field_name="row_attrition_report.final_rows",
        )
        if final_rows > input_rows:
            raise PhosPyInputError(
                "row_attrition_report.final_rows must be less than or equal "
                "to input_rows"
            )
        if records:
            if records[0].input_rows != input_rows:
                raise PhosPyInputError(
                    "row_attrition_report.input_rows must match the first "
                    "record input_rows"
                )
            if records[-1].output_rows != final_rows:
                raise PhosPyInputError(
                    "row_attrition_report.final_rows must match the last "
                    "record output_rows"
                )
            for previous, current in zip(records[:-1], records[1:], strict=True):
                if current.input_rows != previous.output_rows:
                    raise PhosPyInputError(
                        "row_attrition_report records must form a continuous "
                        "row-count chain"
                    )
        elif final_rows != input_rows:
            raise PhosPyInputError(
                "row_attrition_report.final_rows must equal input_rows when "
                "records is empty"
            )
        object.__setattr__(self, "records", records)
        object.__setattr__(self, "input_rows", input_rows)
        object.__setattr__(self, "final_rows", final_rows)

    @classmethod
    def from_records(
        cls,
        records: Sequence[RowAttritionRecord],
    ) -> RowAttritionReport:
        """Create a report using the first input and last output row counts."""

        record_tuple = _required_row_attrition_record_tuple(records)
        if not record_tuple:
            raise PhosPyInputError(
                "row_attrition_report.records must contain at least one record"
            )
        return cls(
            records=record_tuple,
            input_rows=record_tuple[0].input_rows,
            final_rows=record_tuple[-1].output_rows,
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible row-attrition report payload."""

        return {
            "records": [record.to_payload() for record in self.records],
            "input_rows": int(self.input_rows),
            "final_rows": int(self.final_rows),
        }


@dataclass(frozen=True, slots=True)
class TableFingerprint:
    """Deterministic table fingerprint metadata."""

    name: str
    rows: int
    columns: int
    index_name: str | None
    column_names: tuple[str, ...]
    dtypes: tuple[str, ...]
    exact_hash_algorithm: str
    exact_hash_value: str
    tolerance_hash_algorithm: str
    tolerance_hash_value: str
    index_structure: Mapping[str, JsonValue] | None = None
    column_index_structure: Mapping[str, JsonValue] | None = None

    def __post_init__(self) -> None:
        name = _required_provenance_text(self.name, field_name="table_fingerprint.name")
        rows = _required_non_negative_row_count(
            self.rows,
            field_name="table_fingerprint.rows",
        )
        columns = _required_non_negative_row_count(
            self.columns,
            field_name="table_fingerprint.columns",
        )
        column_names = _provenance_string_tuple(
            self.column_names,
            field_name="table_fingerprint.column_names",
        )
        dtypes = _provenance_string_tuple(
            self.dtypes,
            field_name="table_fingerprint.dtypes",
        )
        if len(column_names) != columns:
            raise PhosPyInputError(
                "table_fingerprint.column_names length must match "
                "table_fingerprint.columns"
            )
        if len(dtypes) != columns:
            raise PhosPyInputError(
                "table_fingerprint.dtypes length must match table_fingerprint.columns"
            )
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "rows", rows)
        object.__setattr__(self, "columns", columns)
        object.__setattr__(
            self,
            "index_name",
            None if self.index_name is None else str(self.index_name),
        )
        object.__setattr__(self, "column_names", column_names)
        object.__setattr__(self, "dtypes", dtypes)
        object.__setattr__(
            self,
            "exact_hash_algorithm",
            _required_provenance_text(
                self.exact_hash_algorithm,
                field_name="table_fingerprint.exact_hash_algorithm",
            ),
        )
        object.__setattr__(
            self,
            "exact_hash_value",
            _required_provenance_text(
                self.exact_hash_value,
                field_name="table_fingerprint.exact_hash_value",
            ),
        )
        object.__setattr__(
            self,
            "tolerance_hash_algorithm",
            _required_provenance_text(
                self.tolerance_hash_algorithm,
                field_name="table_fingerprint.tolerance_hash_algorithm",
            ),
        )
        object.__setattr__(
            self,
            "tolerance_hash_value",
            _required_provenance_text(
                self.tolerance_hash_value,
                field_name="table_fingerprint.tolerance_hash_value",
            ),
        )
        object.__setattr__(
            self,
            "index_structure",
            freeze_optional_json_mapping(
                self.index_structure,
                field_name="table_fingerprint.index_structure",
            ),
        )
        object.__setattr__(
            self,
            "column_index_structure",
            freeze_optional_json_mapping(
                self.column_index_structure,
                field_name="table_fingerprint.column_index_structure",
            ),
        )


@dataclass(frozen=True, slots=True)
class EnvironmentProvenance:
    """Runtime environment fingerprint for reproducibility."""

    package_name: str
    package_version: str
    python_version: str
    dependency_versions: dict[str, str | None]
    schema_version: int = ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2
    platform: dict[str, str] = field(default_factory=_empty_platform_provenance)
    blas_lapack: dict[str, JsonValue] = field(default_factory=dict)
    thread_environment: dict[str, str | None] = field(default_factory=dict)
    timezone: str | None = None
    locale: dict[str, str | None] = field(default_factory=dict)
    constraints_fingerprint: dict[str, str | None] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "dependency_versions",
            freeze_json_mapping(
                self.dependency_versions,
                field_name="environment_provenance.dependency_versions",
            ),
        )
        object.__setattr__(
            self,
            "platform",
            freeze_json_mapping(
                self.platform,
                field_name="environment_provenance.platform",
            ),
        )
        object.__setattr__(
            self,
            "blas_lapack",
            freeze_json_mapping(
                self.blas_lapack,
                field_name="environment_provenance.blas_lapack",
            ),
        )
        object.__setattr__(
            self,
            "thread_environment",
            freeze_json_mapping(
                self.thread_environment,
                field_name="environment_provenance.thread_environment",
            ),
        )
        object.__setattr__(
            self,
            "locale",
            freeze_json_mapping(
                self.locale,
                field_name="environment_provenance.locale",
            ),
        )
        object.__setattr__(
            self,
            "constraints_fingerprint",
            freeze_json_mapping(
                self.constraints_fingerprint,
                field_name="environment_provenance.constraints_fingerprint",
            ),
        )


@dataclass(frozen=True, slots=True)
class PreprocessingStageProvenance:
    """Executed preprocessing-stage provenance record."""

    stage: str
    operation: str
    parameters: Mapping[str, object]
    input_shape: tuple[int, int]
    output_shape: tuple[int, int]
    input_hash: str
    output_hash: str
    dropped_row_ids: tuple[str, ...]
    dropped_row_count: int
    phospho_input_hash: str | None = None
    phospho_output_hash: str | None = None
    schema_version: int = PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3
    consumed_input_tables: tuple[TableFingerprint, ...] = ()
    produced_output_tables: tuple[TableFingerprint, ...] = ()
    backend: str | None = None
    random_seed: int | None = None
    determinism: DeterminismKind | str = DeterminismKind.DETERMINISTIC
    reproducibility_caveats: tuple[ReproducibilityCaveat, ...] = ()
    imputed_cell_count: int = 0
    imputed_row_ids: tuple[str, ...] = ()
    notes: str | None = None
    diagnostics: dict[str, JsonValue] | None = None
    batch_correction_provenance: BatchCorrectionProvenance | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "parameters",
            freeze_json_mapping(
                self.parameters,
                field_name="preprocessing_stage.parameters",
            ),
        )
        object.__setattr__(
            self,
            "input_shape",
            _required_shape(
                self.input_shape, field_name="preprocessing_stage.input_shape"
            ),
        )
        object.__setattr__(
            self,
            "output_shape",
            _required_shape(
                self.output_shape,
                field_name="preprocessing_stage.output_shape",
            ),
        )
        object.__setattr__(
            self,
            "dropped_row_ids",
            _provenance_string_tuple(
                self.dropped_row_ids,
                field_name="preprocessing_stage.dropped_row_ids",
            ),
        )
        object.__setattr__(
            self,
            "consumed_input_tables",
            _required_table_fingerprint_tuple(
                self.consumed_input_tables,
                field_name="preprocessing_stage.consumed_input_tables",
            ),
        )
        object.__setattr__(
            self,
            "produced_output_tables",
            _required_table_fingerprint_tuple(
                self.produced_output_tables,
                field_name="preprocessing_stage.produced_output_tables",
            ),
        )
        object.__setattr__(
            self,
            "reproducibility_caveats",
            _required_reproducibility_caveat_tuple(
                self.reproducibility_caveats,
                field_name="preprocessing_stage.reproducibility_caveats",
            ),
        )
        object.__setattr__(
            self,
            "imputed_row_ids",
            _provenance_string_tuple(
                self.imputed_row_ids,
                field_name="preprocessing_stage.imputed_row_ids",
            ),
        )
        object.__setattr__(
            self,
            "diagnostics",
            freeze_optional_json_mapping(
                self.diagnostics,
                field_name="preprocessing_stage.diagnostics",
            ),
        )


@dataclass(frozen=True, slots=True)
class BatchCorrectionRejectedEntity:
    """Rejected row, site, or sample recorded for correction provenance."""

    entity_type: str
    identifier: str
    reason: str
    details: Mapping[str, JsonValue] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "details",
            freeze_optional_json_mapping(
                self.details,
                field_name="batch_correction_rejected_entity.details",
            ),
        )


@dataclass(frozen=True, slots=True)
class BatchCorrectionProvenance:
    """Planned or executed batch-correction provenance record.

    This model is an audit structure only. It records requested intent, resolved
    inputs, fingerprints, diagnostics, and rejection reasons without selecting
    controls, validating scientific eligibility, or modifying matrices.
    """

    requested_method: str
    resolved_parameters: Mapping[str, JsonValue]
    preprocessing_stage_order: tuple[str, ...]
    control_site_source: Mapping[str, JsonValue]
    selected_site_key_rows: tuple[str, ...]
    batch_metadata: Mapping[str, JsonValue]
    replicate_metadata: Mapping[str, JsonValue] | None
    design_metadata: Mapping[str, JsonValue]
    missing_value_policy: Mapping[str, JsonValue]
    observation_masks: tuple[TableFingerprint, ...]
    input_matrix_fingerprint: TableFingerprint
    output_matrix_fingerprint: TableFingerprint | None
    diagnostics: Mapping[str, JsonValue] = field(default_factory=_empty_json_mapping)
    warnings: tuple[str, ...] = ()
    rejected_entities: tuple[BatchCorrectionRejectedEntity, ...] = ()
    phospy_version: str = "unknown"
    python_version: str = "unknown"
    dependency_versions: Mapping[str, str | None] = field(default_factory=dict)
    schema_version: int = BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1
    imputation_policy: Mapping[str, JsonValue] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "resolved_parameters",
            freeze_json_mapping(
                self.resolved_parameters,
                field_name="batch_correction_provenance.resolved_parameters",
            ),
        )
        object.__setattr__(
            self,
            "preprocessing_stage_order",
            _provenance_string_tuple(
                self.preprocessing_stage_order,
                field_name="batch_correction_provenance.preprocessing_stage_order",
            ),
        )
        object.__setattr__(
            self,
            "control_site_source",
            freeze_json_mapping(
                self.control_site_source,
                field_name="batch_correction_provenance.control_site_source",
            ),
        )
        object.__setattr__(
            self,
            "selected_site_key_rows",
            _provenance_string_tuple(
                self.selected_site_key_rows,
                field_name="batch_correction_provenance.selected_site_key_rows",
            ),
        )
        object.__setattr__(
            self,
            "batch_metadata",
            freeze_json_mapping(
                self.batch_metadata,
                field_name="batch_correction_provenance.batch_metadata",
            ),
        )
        object.__setattr__(
            self,
            "replicate_metadata",
            freeze_optional_json_mapping(
                self.replicate_metadata,
                field_name="batch_correction_provenance.replicate_metadata",
            ),
        )
        object.__setattr__(
            self,
            "design_metadata",
            freeze_json_mapping(
                self.design_metadata,
                field_name="batch_correction_provenance.design_metadata",
            ),
        )
        object.__setattr__(
            self,
            "missing_value_policy",
            freeze_json_mapping(
                self.missing_value_policy,
                field_name="batch_correction_provenance.missing_value_policy",
            ),
        )
        object.__setattr__(
            self,
            "observation_masks",
            _required_table_fingerprint_tuple(
                self.observation_masks,
                field_name="batch_correction_provenance.observation_masks",
            ),
        )
        if not isinstance(self.input_matrix_fingerprint, TableFingerprint):
            raise PhosPyInputError(
                "batch_correction_provenance.input_matrix_fingerprint must be "
                "a TableFingerprint"
            )
        if self.output_matrix_fingerprint is not None and not isinstance(
            self.output_matrix_fingerprint,
            TableFingerprint,
        ):
            raise PhosPyInputError(
                "batch_correction_provenance.output_matrix_fingerprint must be "
                "a TableFingerprint or None"
            )
        object.__setattr__(
            self,
            "diagnostics",
            freeze_json_mapping(
                self.diagnostics,
                field_name="batch_correction_provenance.diagnostics",
            ),
        )
        object.__setattr__(
            self,
            "warnings",
            _provenance_string_tuple(
                self.warnings,
                field_name="batch_correction_provenance.warnings",
            ),
        )
        object.__setattr__(
            self,
            "rejected_entities",
            _required_batch_correction_rejected_entity_tuple(
                self.rejected_entities,
            ),
        )
        object.__setattr__(
            self,
            "dependency_versions",
            freeze_json_mapping(
                self.dependency_versions,
                field_name="batch_correction_provenance.dependency_versions",
            ),
        )
        object.__setattr__(
            self,
            "imputation_policy",
            freeze_json_mapping(
                self.imputation_policy,
                field_name="batch_correction_provenance.imputation_policy",
            ),
        )


@dataclass(frozen=True, slots=True)
class ReferenceProvenance:
    """Resolved reference identity and table fingerprints."""

    source_type: str
    organism: str
    bundle_id: str | None
    table_fingerprints: tuple[TableFingerprint, ...]
    source_name: str | None = None
    source_version: str | None = None
    retrieved_at: str | None = None
    identifier_namespace: str | None = None
    sequence_window: Mapping[str, JsonValue] | None = None
    manifest: Mapping[str, JsonValue] | None = None
    identifier_normalisation: ReferenceIdentifierNormalisationReport | None = None
    reference_context: ReferenceContextProtocol | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "table_fingerprints",
            _required_table_fingerprint_tuple(
                self.table_fingerprints,
                field_name="reference_provenance.table_fingerprints",
            ),
        )
        object.__setattr__(
            self,
            "sequence_window",
            freeze_optional_json_mapping(
                self.sequence_window,
                field_name="reference_provenance.sequence_window",
            ),
        )
        object.__setattr__(
            self,
            "manifest",
            freeze_optional_json_mapping(
                self.manifest,
                field_name="reference_provenance.manifest",
            ),
        )
        source_version = _optional_provenance_text(self.source_version)
        object.__setattr__(self, "source_version", source_version)
        validate_reference_source_version_agreement(
            (
                ("provenance.source_version", source_version),
                (
                    "reference_context.source_version",
                    _reference_context_source_version(self.reference_context),
                ),
                (
                    "manifest.source_version",
                    _manifest_source_version(self.manifest),
                ),
            )
        )


@dataclass(frozen=True, slots=True)
class KinaseLibraryResourceProvenance:
    """Resolved provenance for local Kinase Library-style motif resources."""

    source_type: str
    source_name: str
    source_version: str
    license: str
    score_scale: str
    organisms: tuple[str, ...]
    sequence_window: Mapping[str, JsonValue]
    source_files: Mapping[str, JsonValue]
    table_fingerprints: tuple[TableFingerprint, ...]
    retrieved_at: str | None = None
    manifest: Mapping[str, JsonValue] | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "organisms", tuple(self.organisms))
        object.__setattr__(
            self,
            "sequence_window",
            freeze_json_mapping(
                self.sequence_window,
                field_name="kinase_library_resource_provenance.sequence_window",
            ),
        )
        object.__setattr__(
            self,
            "source_files",
            freeze_json_mapping(
                self.source_files,
                field_name="kinase_library_resource_provenance.source_files",
            ),
        )
        object.__setattr__(
            self,
            "table_fingerprints",
            _required_table_fingerprint_tuple(
                self.table_fingerprints,
                field_name="kinase_library_resource_provenance.table_fingerprints",
            ),
        )
        object.__setattr__(
            self,
            "manifest",
            freeze_optional_json_mapping(
                self.manifest,
                field_name="kinase_library_resource_provenance.manifest",
            ),
        )


@dataclass(frozen=True, slots=True)
class RunProvenance:
    """Machine-readable run provenance payload.

    Workflow-specific audit details that are not table-transforming stages, such
    as protein-aware preparation summaries, belong in `workflow_parameters`.
    `reference_context` records the biological reference context of the input
    dataset for downstream compatibility checks. `reference` remains the
    workflow reference resource provenance when a workflow consumes an explicit
    reference bundle.
    """

    environment: EnvironmentProvenance
    input_tables: tuple[TableFingerprint, ...]
    preprocessing_stages: tuple[PreprocessingStageProvenance, ...]
    reference: ReferenceProvenance | None
    workflow_name: str | None
    workflow_parameters: Mapping[str, object]
    random_state: int | None
    random_seed_policy: str | None
    output_tables: tuple[TableFingerprint, ...]
    scientific_policies: tuple[ScientificPolicyRecord, ...] = ()
    reference_context: ReferenceContextProtocol | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.environment, EnvironmentProvenance):
            raise PhosPyInputError(
                "run_provenance.environment must be EnvironmentProvenance"
            )
        object.__setattr__(
            self,
            "input_tables",
            _required_table_fingerprint_tuple(
                self.input_tables,
                field_name="run_provenance.input_tables",
            ),
        )
        object.__setattr__(
            self,
            "preprocessing_stages",
            _required_preprocessing_stage_tuple(self.preprocessing_stages),
        )
        if self.reference is not None and not isinstance(
            self.reference,
            ReferenceProvenance,
        ):
            raise PhosPyInputError(
                "run_provenance.reference must be ReferenceProvenance or None"
            )
        object.__setattr__(
            self,
            "workflow_parameters",
            freeze_json_mapping(
                self.workflow_parameters,
                field_name="run_provenance.workflow_parameters",
            ),
        )
        object.__setattr__(
            self,
            "output_tables",
            _required_table_fingerprint_tuple(
                self.output_tables,
                field_name="run_provenance.output_tables",
            ),
        )
        object.__setattr__(
            self,
            "scientific_policies",
            _required_scientific_policy_tuple(self.scientific_policies),
        )


def _required_non_negative_row_count(value: object, *, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise PhosPyInputError(f"{field_name} must be a non-negative integer")
    count = int(value)
    if count < 0:
        raise PhosPyInputError(f"{field_name} must be non-negative")
    return count


def _required_provenance_text_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(
        values,
        Sequence,
    ):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    return tuple(
        _required_provenance_text(value, field_name=f"{field_name}[]")
        for value in values
    )


def _provenance_string_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(
        values,
        Sequence,
    ):
        raise PhosPyInputError(f"{field_name} must be a sequence of strings")
    return tuple(str(value) for value in values)


def _required_shape(value: object, *, field_name: str) -> tuple[int, int]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise PhosPyInputError(f"{field_name} must contain exactly two integers")
    sequence = tuple(value)
    if len(sequence) != 2:
        raise PhosPyInputError(f"{field_name} must contain exactly two integers")
    return (
        _required_non_negative_row_count(sequence[0], field_name=f"{field_name}[0]"),
        _required_non_negative_row_count(sequence[1], field_name=f"{field_name}[1]"),
    )


def _required_row_attrition_record_tuple(
    records: object,
) -> tuple[RowAttritionRecord, ...]:
    if isinstance(records, (str, bytes, bytearray)) or not isinstance(
        records,
        Sequence,
    ):
        raise PhosPyInputError(
            "row_attrition_report.records must be a sequence of "
            "RowAttritionRecord values"
        )
    record_tuple = tuple(records)
    for record in record_tuple:
        if not isinstance(record, RowAttritionRecord):
            raise PhosPyInputError(
                "row_attrition_report.records must contain only "
                "RowAttritionRecord values"
            )
    return record_tuple


def _required_table_fingerprint_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    fingerprints = tuple(values)
    for fingerprint in fingerprints:
        if not isinstance(fingerprint, TableFingerprint):
            raise PhosPyInputError(
                f"{field_name} must contain only TableFingerprint values"
            )
    return fingerprints


def _required_reproducibility_caveat_tuple(
    values: object,
    *,
    field_name: str,
) -> tuple[ReproducibilityCaveat, ...]:
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    caveats = tuple(values)
    for caveat in caveats:
        if not isinstance(caveat, ReproducibilityCaveat):
            raise PhosPyInputError(
                f"{field_name} must contain only ReproducibilityCaveat values"
            )
    return caveats


def _required_batch_correction_rejected_entity_tuple(
    values: object,
) -> tuple[BatchCorrectionRejectedEntity, ...]:
    field_name = "batch_correction_provenance.rejected_entities"
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    entities = tuple(values)
    for entity in entities:
        if not isinstance(entity, BatchCorrectionRejectedEntity):
            raise PhosPyInputError(
                f"{field_name} must contain only BatchCorrectionRejectedEntity values"
            )
    return entities


def _required_preprocessing_stage_tuple(
    values: object,
) -> tuple[PreprocessingStageProvenance, ...]:
    field_name = "run_provenance.preprocessing_stages"
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    stages = tuple(values)
    for stage in stages:
        if not isinstance(stage, PreprocessingStageProvenance):
            raise PhosPyInputError(
                f"{field_name} must contain only PreprocessingStageProvenance values"
            )
    return stages


def _required_scientific_policy_tuple(
    values: object,
) -> tuple[ScientificPolicyRecord, ...]:
    field_name = "run_provenance.scientific_policies"
    if isinstance(values, (str, bytes, bytearray)) or not isinstance(values, Sequence):
        raise PhosPyInputError(f"{field_name} must be a sequence")
    policies = tuple(values)
    for policy in policies:
        if not isinstance(policy, ScientificPolicyRecord):
            raise PhosPyInputError(
                f"{field_name} must contain only ScientificPolicyRecord values"
            )
    return policies


def validate_reference_source_version_agreement(
    entries: Sequence[tuple[str, object | None]],
) -> None:
    known: list[tuple[str, str]] = []
    for label, value in entries:
        normalized = _known_reference_source_version(value)
        if normalized is not None:
            known.append((label, normalized))
    if len(known) < 2:
        return
    baseline_label, baseline_value = known[0]
    for label, value in known[1:]:
        if value == baseline_value:
            continue
        raise ReferenceValidationError(
            "Reference provenance source-version mismatch:\n"
            f"{baseline_label}={baseline_value!r},\n"
            f"{label}={value!r}"
        )


def _known_reference_source_version(value: object | None) -> str | None:
    if not isinstance(value, str):
        return None
    text = value.strip()
    return text or None


def _reference_context_source_version(reference_context: object | None) -> str | None:
    if reference_context is None:
        return None
    return _known_reference_source_version(
        getattr(reference_context, "source_version", None)
    )


def _manifest_source_version(
    manifest: Mapping[str, JsonValue] | None,
) -> str | None:
    if not isinstance(manifest, Mapping):
        return None
    return _known_reference_source_version(manifest.get("source_version"))


__all__ = [
    "BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1",
    "BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS",
    "BatchCorrectionProvenance",
    "BatchCorrectionRejectedEntity",
    "DeterminismKind",
    "EnvironmentProvenance",
    "InputIntensityScaleEvidence",
    "JsonValue",
    "KinaseLibraryResourceProvenance",
    "PREPROCESSING_EXTERNAL_NONDETERMINISM_CAVEAT_CODE",
    "PREPROCESSING_STAGE_DETERMINISM_DETERMINISTIC",
    "PREPROCESSING_STAGE_DETERMINISM_EXTERNALLY_NONDETERMINISTIC",
    "PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY",
    "PREPROCESSING_STAGE_DETERMINISM_PURE",
    "PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC",
    "PreprocessingStageProvenance",
    "ReferenceProvenance",
    "ReferenceContextProtocol",
    "ReproducibilityCaveat",
    "RowAttritionRecord",
    "RowAttritionReport",
    "RunProvenance",
    "TableFingerprint",
    "TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V1",
    "TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V2",
    "TrustedDatasetConstructionAssertions",
    "TrustedDatasetConstructionEvidence",
]
