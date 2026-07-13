"""Typed machine-readable provenance models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from numbers import Integral
from typing import TYPE_CHECKING

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord

if TYPE_CHECKING:
    from phospy.science.references.identifiers import (
        ReferenceIdentifierNormalisationReport,
    )
    from phospy.science.references.models import ReferenceContext

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V1 = 1
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2 = 2
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3 = 3
ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V1 = 1
ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2 = 2
BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1 = 1


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
JsonValue = (
    JsonPrimitive | list["JsonValue"] | tuple["JsonValue", ...] | dict[str, "JsonValue"]
)


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
        if not isinstance(self.details, Mapping):
            raise PhosPyInputError("reproducibility_caveat.details must be a mapping")
        object.__setattr__(
            self,
            "details",
            {str(key): value for key, value in self.details.items()},
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible caveat payload."""

        return {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "details": dict(self.details),
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


def _optional_provenance_text(value: object | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


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


@dataclass(frozen=True, slots=True)
class BatchCorrectionRejectedEntity:
    """Rejected row, site, or sample recorded for correction provenance."""

    entity_type: str
    identifier: str
    reason: str
    details: Mapping[str, JsonValue] | None = None


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
    reference_context: ReferenceContext | None = None

    def __post_init__(self) -> None:
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
    reference_context: ReferenceContext | None = None


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
    "ReproducibilityCaveat",
    "RowAttritionRecord",
    "RowAttritionReport",
    "RunProvenance",
    "TableFingerprint",
]
