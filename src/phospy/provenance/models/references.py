"""Reference-resource provenance and compatibility models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Protocol

from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.immutability import (
    freeze_json_mapping,
    freeze_optional_json_mapping,
)
from phospy.provenance.models._shared import JsonValue, _optional_provenance_text
from phospy.provenance.models.tables import (
    TableFingerprint,
    _required_table_fingerprint_tuple,
)
from phospy.provenance.organisms import Organism, normalize_organism
from phospy.provenance.reference_identifiers import (
    ReferenceIdentifierNormalisationReport,
)


class ReferenceContextProtocol(Protocol):
    """Structural protocol for reference-context provenance values."""

    @property
    def organism(self) -> Organism: ...

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


@dataclass(frozen=True, slots=True, init=False)
class ReferenceProvenance:
    """Resolved reference identity and table fingerprints."""

    source_type: str
    organism: Organism
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

    def __init__(
        self,
        source_type: str,
        organism: object,
        bundle_id: str | None,
        table_fingerprints: tuple[TableFingerprint, ...],
        source_name: str | None = None,
        source_version: str | None = None,
        retrieved_at: str | None = None,
        identifier_namespace: str | None = None,
        sequence_window: Mapping[str, JsonValue] | None = None,
        manifest: Mapping[str, JsonValue] | None = None,
        identifier_normalisation: ReferenceIdentifierNormalisationReport | None = None,
        reference_context: ReferenceContextProtocol | None = None,
    ) -> None:
        object.__setattr__(self, "source_type", source_type)
        object.__setattr__(self, "organism", organism)
        object.__setattr__(self, "bundle_id", bundle_id)
        object.__setattr__(self, "table_fingerprints", table_fingerprints)
        object.__setattr__(self, "source_name", source_name)
        object.__setattr__(self, "source_version", source_version)
        object.__setattr__(self, "retrieved_at", retrieved_at)
        object.__setattr__(self, "identifier_namespace", identifier_namespace)
        object.__setattr__(self, "sequence_window", sequence_window)
        object.__setattr__(self, "manifest", manifest)
        object.__setattr__(self, "identifier_normalisation", identifier_normalisation)
        object.__setattr__(self, "reference_context", reference_context)
        self.__post_init__()

    def __post_init__(self) -> None:
        organism = normalize_organism(
            self.organism,
            field_name="reference_provenance.organism",
            error_type=ReferenceValidationError,
        )
        object.__setattr__(self, "organism", organism)
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
        _require_reference_provenance_organism_coherence(
            organism=organism,
            reference_context=self.reference_context,
            manifest=self.manifest,
        )
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


def _require_reference_provenance_organism_coherence(
    *,
    organism: Organism,
    reference_context: ReferenceContextProtocol | None,
    manifest: Mapping[str, JsonValue] | None,
) -> None:
    entries: list[tuple[str, object]] = [("reference_provenance.organism", organism)]
    if reference_context is not None:
        entries.append(
            (
                "reference_provenance.reference_context.organism",
                reference_context.organism,
            )
        )
    manifest_organism = _manifest_organism(manifest)
    if manifest_organism is not None:
        entries.append(("reference_provenance.manifest.organism", manifest_organism))
    _require_organism_identity_agreement(
        entries=entries,
        conflict_prefix="Reference provenance organism mismatch",
    )


def _require_run_provenance_reference_context_organism_coherence(
    *,
    reference: ReferenceProvenance | None,
    reference_context: ReferenceContextProtocol | None,
) -> None:
    if reference is None or reference_context is None:
        return
    _require_organism_identity_agreement(
        entries=[
            ("run_provenance.reference.organism", reference.organism),
            ("run_provenance.reference_context.organism", reference_context.organism),
        ],
        conflict_prefix="Run provenance reference-context organism mismatch",
    )


def _require_organism_identity_agreement(
    *,
    entries: list[tuple[str, object]],
    conflict_prefix: str,
) -> None:
    if not entries:
        return
    normalized = [
        (
            field_name,
            normalize_organism(
                value,
                field_name=field_name,
                error_type=ReferenceValidationError,
            ),
            value,
        )
        for field_name, value in entries
    ]
    expected_field, expected_organism, _ = normalized[0]
    conflicts = [
        (field_name, organism, raw_value)
        for field_name, organism, raw_value in normalized[1:]
        if organism is not expected_organism
    ]
    if not conflicts:
        return
    conflict_text = "; ".join(
        f"{field_name}={_format_organism_value(raw_value)!r}"
        f" resolved_to={organism.value!r}"
        for field_name, organism, raw_value in conflicts
    )
    raise ReferenceValidationError(
        f"{conflict_prefix}: {expected_field}={expected_organism.value!r}; "
        f"{conflict_text}"
    )


def _manifest_organism(manifest: Mapping[str, JsonValue] | None) -> object | None:
    if not isinstance(manifest, Mapping):
        return None
    organism_common_name = manifest.get("organism_common_name")
    if isinstance(organism_common_name, str) and organism_common_name.strip():
        return organism_common_name
    organism = manifest.get("organism")
    if isinstance(organism, str) and organism.strip():
        return organism
    return None


def _format_organism_value(value: object) -> str:
    if isinstance(value, Organism):
        return value.value
    return str(value)


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
