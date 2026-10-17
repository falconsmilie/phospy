"""Immutable structural audit records for batch correction."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import (
    freeze_json_mapping,
    freeze_optional_json_mapping,
    thaw_json_mapping,
)
from phospy.provenance.models._shared import (
    JsonValue,
    _empty_json_mapping,
    _provenance_string_tuple,
)
from phospy.provenance.models.tables import (
    TableFingerprint,
    _required_table_fingerprint_tuple,
)

BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1 = 1

BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS = frozenset(
    {"not_provided", "unknown", "none", "null"}
)

_BATCH_CORRECTION_PROVENANCE_JSON_MAPPING_FIELDS = frozenset(
    {
        "resolved_parameters",
        "control_site_source",
        "batch_metadata",
        "design_metadata",
        "missing_value_policy",
        "diagnostics",
        "dependency_versions",
        "imputation_policy",
    }
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
    """Planned or executed batch-correction audit record.

    This immutable wire contract records requested intent, resolved values,
    fingerprints, diagnostics and rejection reasons. It does not select
    controls, determine scientific eligibility or modify matrices.
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

    def __getattribute__(self, name: str) -> object:
        if name in _BATCH_CORRECTION_PROVENANCE_JSON_MAPPING_FIELDS:
            value = object.__getattribute__(self, name)
            return thaw_json_mapping(
                value,
                field_name=f"batch_correction_provenance.{name}",
            )
        if name == "replicate_metadata":
            value = object.__getattribute__(self, name)
            if value is None:
                return None
            return thaw_json_mapping(
                value,
                field_name="batch_correction_provenance.replicate_metadata",
            )
        return object.__getattribute__(self, name)

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


__all__ = [
    "BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1",
    "BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS",
    "BatchCorrectionProvenance",
    "BatchCorrectionRejectedEntity",
]
