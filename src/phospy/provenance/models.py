"""Typed machine-readable provenance models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from phospy.provenance.scientific_policy_models import ScientificPolicyRecord

if TYPE_CHECKING:
    from phospy.references.identifiers import ReferenceIdentifierNormalisationReport

PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V1 = 1
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V2 = 2
PREPROCESSING_STAGE_PROVENANCE_SCHEMA_VERSION_V3 = 3
ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V1 = 1
ENVIRONMENT_PROVENANCE_SCHEMA_VERSION_V2 = 2
PREPROCESSING_STAGE_DETERMINISM_PURE = "pure"
PREPROCESSING_STAGE_DETERMINISM_SEEDED_STOCHASTIC = "seeded_stochastic"
PREPROCESSING_STAGE_DETERMINISM_EXTERNAL_DEPENDENCY = "external_dependency"

JsonPrimitive = str | int | float | bool | None
JsonValue = (
    JsonPrimitive | list["JsonValue"] | tuple["JsonValue", ...] | dict[str, "JsonValue"]
)


def _empty_platform_provenance() -> dict[str, str]:
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
    hash_algorithm: str
    hash_value: str
    exact_hash_algorithm: str | None = None
    exact_hash_value: str | None = None
    tolerance_hash_algorithm: str | None = None
    tolerance_hash_value: str | None = None
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
    is_deterministic: bool = True
    imputed_cell_count: int = 0
    imputed_row_ids: tuple[str, ...] = ()
    notes: str | None = None
    diagnostics: dict[str, JsonValue] | None = None


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
class RunProvenance:
    """Machine-readable run provenance payload."""

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
    "EnvironmentProvenance",
    "JsonValue",
    "PreprocessingStageProvenance",
    "ReferenceProvenance",
    "RunProvenance",
    "TableFingerprint",
]
