"""Typed machine-readable provenance models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from phospy.provenance.scientific_policy_models import ScientificPolicyRecord

if TYPE_CHECKING:
    from phospy.science.references.identifiers import (
        ReferenceIdentifierNormalisationReport,
    )

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V1 = 1
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2 = 2
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3 = 3
ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V1 = 1
ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2 = 2
BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1 = 1
PREPROCESSING_STAGE_DETERMINISM_PURE = "pure"
PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC = "seeded_stochastic"
PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY = "external_dependency"
BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS = frozenset(
    {"not_provided", "unknown", "none", "null"}
)

JsonPrimitive = str | int | float | bool | None
JsonValue = (
    JsonPrimitive | list["JsonValue"] | tuple["JsonValue", ...] | dict[str, "JsonValue"]
)


def _empty_platform_provenance() -> dict[str, str]:
    return {}


def _empty_json_mapping() -> dict[str, JsonValue]:
    return {}


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
    determinism: str = PREPROCESSING_STAGE_DETERMINISM_PURE
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


__all__ = [
    "BATCH_CORRECTION_PROVENANCE_SCHEMA_VERSION_V1",
    "BATCH_CORRECTION_SELECTED_SITE_KEY_ROW_SENTINELS",
    "BatchCorrectionProvenance",
    "BatchCorrectionRejectedEntity",
    "EnvironmentProvenance",
    "JsonValue",
    "KinaseLibraryResourceProvenance",
    "PreprocessingStageProvenance",
    "ReferenceProvenance",
    "RunProvenance",
    "TableFingerprint",
]
