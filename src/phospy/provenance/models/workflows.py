"""Workflow and preprocessing-stage provenance models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import (
    freeze_json_mapping,
    freeze_optional_json_mapping,
    thaw_json_mapping,
)
from phospy.provenance.models._shared import (
    JsonValue,
    _empty_json_mapping,
    _optional_provenance_text,
    _provenance_string_tuple,
    _required_provenance_text,
    _required_shape,
)
from phospy.provenance.models.environment import EnvironmentProvenance
from phospy.provenance.models.references import (
    ReferenceContextProtocol,
    ReferenceProvenance,
    _require_run_provenance_reference_context_organism_coherence,
)
from phospy.provenance.models.tables import (
    TableFingerprint,
    _required_table_fingerprint_tuple,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V1 = 1

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2 = 2

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3 = 3

BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1 = 1

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
    _frozen_parameters: Mapping[str, JsonValue] = field(
        init=False,
        repr=False,
        compare=False,
    )
    _frozen_diagnostics: Mapping[str, JsonValue] | None = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __getattribute__(self, name: str) -> object:
        if name == "parameters":
            try:
                frozen = object.__getattribute__(self, "_frozen_parameters")
            except AttributeError:
                return object.__getattribute__(self, name)
            return thaw_json_mapping(
                frozen, field_name="preprocessing_stage.parameters"
            )
        if name == "diagnostics":
            try:
                frozen_optional = object.__getattribute__(
                    self,
                    "_frozen_diagnostics",
                )
            except AttributeError:
                return object.__getattribute__(self, name)
            if frozen_optional is None:
                return None
            return thaw_json_mapping(
                frozen_optional,
                field_name="preprocessing_stage.diagnostics",
            )
        return object.__getattribute__(self, name)

    def __post_init__(self) -> None:
        raw_parameters = object.__getattribute__(self, "parameters")
        frozen_parameters = freeze_json_mapping(
            raw_parameters,
            field_name="preprocessing_stage.parameters",
        )
        object.__setattr__(self, "_frozen_parameters", frozen_parameters)
        object.__setattr__(
            self,
            "parameters",
            frozen_parameters,
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
        raw_diagnostics = object.__getattribute__(self, "diagnostics")
        frozen_diagnostics = freeze_optional_json_mapping(
            raw_diagnostics,
            field_name="preprocessing_stage.diagnostics",
        )
        object.__setattr__(self, "_frozen_diagnostics", frozen_diagnostics)
        object.__setattr__(self, "diagnostics", frozen_diagnostics)


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
    _frozen_workflow_parameters: Mapping[str, JsonValue] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __getattribute__(self, name: str) -> object:
        if name == "workflow_parameters":
            try:
                frozen = object.__getattribute__(self, "_frozen_workflow_parameters")
            except AttributeError:
                return object.__getattribute__(self, name)
            return thaw_json_mapping(
                frozen,
                field_name="run_provenance.workflow_parameters",
            )
        return object.__getattribute__(self, name)

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
        raw_workflow_parameters = object.__getattribute__(self, "workflow_parameters")
        frozen_workflow_parameters = freeze_json_mapping(
            raw_workflow_parameters,
            field_name="run_provenance.workflow_parameters",
        )
        object.__setattr__(
            self,
            "_frozen_workflow_parameters",
            frozen_workflow_parameters,
        )
        object.__setattr__(
            self,
            "workflow_parameters",
            frozen_workflow_parameters,
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
        _require_run_provenance_reference_context_organism_coherence(
            reference=self.reference,
            reference_context=self.reference_context,
        )


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
