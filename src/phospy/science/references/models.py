"""Reference domain models."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import InitVar, dataclass, field
from datetime import date
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import cast

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.provenance.hashing import fingerprint_table, hash_json_payload
from phospy.provenance.models import JsonValue, ReferenceProvenance
from phospy.science.references.manifest import (
    RedistributionStatus as RedistributionStatus,
)
from phospy.science.references.manifest import (
    ReferenceFileManifest as ReferenceFileManifest,
)
from phospy.science.references.manifest import (
    ReferenceManifest,
    SequenceWindowDefinition,
)


class Organism(str, Enum):
    """Public organism identifiers used in dataset/reference contracts.

    Enum membership defines contract syntax. Bundled runtime scientific support
    may be narrower in a given release.
    """

    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"


class ReferencePreset(str, Enum):
    """Built-in organism presets for bundled-reference resolution.

    Enum values define public organism lanes accepted by request contracts.
    Bundled runtime references may cover only a subset in a given release.
    """

    AUTO = "auto"
    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"


ReferenceBuildPath = str | Path | PathLike[str]
_REFERENCE_CONTEXT_ID_PREFIX = "reference-context-v1:"


def _empty_json_dict() -> dict[str, JsonValue]:
    return {}


@dataclass(frozen=True, slots=True)
class ReferenceContext:
    """Comparable biological reference identity context."""

    organism: str
    protein_namespace: str
    source_name: str
    source_version: str
    proteome_version: str | None
    reference_table_sha256: str | None
    reference_context_id: str = field(init=False, compare=False)

    def __post_init__(self) -> None:
        organism = _required_reference_context_text(
            self.organism,
            field_name="reference_context.organism",
        ).lower()
        protein_namespace = _required_reference_context_text(
            self.protein_namespace,
            field_name="reference_context.protein_namespace",
        )
        source_name = _required_reference_context_text(
            self.source_name,
            field_name="reference_context.source_name",
        )
        source_version = _required_reference_context_text(
            self.source_version,
            field_name="reference_context.source_version",
        )
        proteome_version = _optional_reference_context_text(self.proteome_version)
        reference_table_sha256 = _optional_reference_context_text(
            self.reference_table_sha256
        )
        if reference_table_sha256 is not None:
            reference_table_sha256 = reference_table_sha256.lower()
        object.__setattr__(self, "organism", organism)
        object.__setattr__(self, "protein_namespace", protein_namespace)
        object.__setattr__(self, "source_name", source_name)
        object.__setattr__(self, "source_version", source_version)
        object.__setattr__(self, "proteome_version", proteome_version)
        object.__setattr__(
            self,
            "reference_table_sha256",
            reference_table_sha256,
        )
        object.__setattr__(
            self,
            "reference_context_id",
            _REFERENCE_CONTEXT_ID_PREFIX + hash_json_payload(self._identity_payload()),
        )

    @classmethod
    def from_manifest(cls, manifest: ReferenceManifest) -> ReferenceContext:
        """Build a reference context from manifest identity metadata."""

        source_version = manifest.source_version
        if source_version is None:
            raise ReferenceValidationError(
                "reference_context.source_version must be non-empty"
            )
        return cls(
            organism=manifest.organism_common_name or manifest.organism,
            protein_namespace=manifest.protein_namespace,
            source_name=manifest.source_name,
            source_version=source_version,
            proteome_version=None,
            reference_table_sha256=manifest.table_sha256,
        )

    @classmethod
    def from_manifest_payload(
        cls,
        payload: Mapping[str, JsonValue],
    ) -> ReferenceContext | None:
        """Build a reference context from serialized manifest identity metadata."""

        return _reference_context_from_manifest_payload(payload)

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> ReferenceContext:
        """Deserialize a reference context payload."""

        return cls(
            organism=_payload_required_text(
                payload,
                "organism",
                field_name="reference_context.organism",
            ),
            protein_namespace=_payload_required_text(
                payload,
                "protein_namespace",
                field_name="reference_context.protein_namespace",
            ),
            source_name=_payload_required_text(
                payload,
                "source_name",
                field_name="reference_context.source_name",
            ),
            source_version=_payload_required_text(
                payload,
                "source_version",
                field_name="reference_context.source_version",
            ),
            proteome_version=_payload_optional_text(payload.get("proteome_version")),
            reference_table_sha256=_payload_optional_text(
                payload.get("reference_table_sha256")
            ),
        )

    def to_payload(self) -> dict[str, JsonValue]:
        """Return a JSON-compatible reference-context payload."""

        payload = self._identity_payload()
        payload["reference_context_id"] = self.reference_context_id
        return payload

    def _identity_payload(self) -> dict[str, JsonValue]:
        return {
            "organism": self.organism,
            "protein_namespace": self.protein_namespace,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "proteome_version": self.proteome_version,
            "reference_table_sha256": self.reference_table_sha256,
        }


def reference_context_from_provenance(
    provenance: ReferenceProvenance,
    *,
    manifest: ReferenceManifest | None = None,
) -> ReferenceContext | None:
    """Resolve a reference context from provenance or richer manifest metadata."""

    if provenance.reference_context is not None:
        return provenance.reference_context
    if manifest is not None:
        return reference_context_from_manifest_if_complete(manifest)
    if provenance.manifest is not None:
        return _reference_context_from_manifest_payload(provenance.manifest)
    if (
        provenance.source_name is None
        or provenance.source_version is None
        or provenance.identifier_namespace is None
    ):
        return None
    return ReferenceContext(
        organism=provenance.organism,
        protein_namespace=provenance.identifier_namespace,
        source_name=provenance.source_name,
        source_version=provenance.source_version,
        proteome_version=None,
        reference_table_sha256=None,
    )


def reference_context_from_manifest_if_complete(
    manifest: ReferenceManifest,
) -> ReferenceContext | None:
    """Build reference context only when manifest identity metadata is complete."""

    source_version = manifest.source_version
    if not isinstance(source_version, str) or not source_version.strip():
        return None
    return ReferenceContext.from_manifest(manifest)


def _reference_context_from_manifest_payload(
    payload: Mapping[str, JsonValue],
) -> ReferenceContext | None:
    organism = _payload_optional_text(
        payload.get("organism_common_name")
    ) or _payload_optional_text(payload.get("organism"))
    protein_namespace = _payload_optional_text(
        payload.get("protein_namespace")
    ) or _payload_optional_text(payload.get("identifier_namespace"))
    source_name = _payload_optional_text(payload.get("source_name"))
    source_version = _payload_optional_text(payload.get("source_version"))
    if (
        organism is None
        or protein_namespace is None
        or source_name is None
        or source_version is None
    ):
        return None
    return ReferenceContext(
        organism=organism,
        protein_namespace=protein_namespace,
        source_name=source_name,
        source_version=source_version,
        proteome_version=_payload_optional_text(payload.get("proteome_version")),
        reference_table_sha256=_payload_optional_text(payload.get("table_sha256")),
    )


def _required_reference_context_text(value: object, *, field_name: str) -> str:
    if value is None:
        raise ReferenceValidationError(f"{field_name} must be non-empty")
    text = str(value).strip()
    if not text:
        raise ReferenceValidationError(f"{field_name} must be non-empty")
    return text


def _optional_reference_context_text(value: object | None) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _payload_required_text(
    payload: Mapping[str, object],
    key: str,
    *,
    field_name: str,
) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ReferenceValidationError(f"{field_name} must be non-empty")
    return value.strip()


def _payload_optional_text(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        return None
    text = value.strip()
    return text or None


@dataclass(frozen=True, slots=True)
class ReferenceBundleBuildRequest:
    """Request for building a local-source reference bundle.

    Construction stores caller intent only. File existence, source metadata,
    column mapping, organism compatibility, and reference validity are enforced
    by ``ReferenceBundleBuilder.run(...)``.
    """

    organism: Organism
    kinase_substrate_path: ReferenceBuildPath
    site_sequence_path: ReferenceBuildPath
    source_name: str
    source_version: str
    retrieved_at: date | str
    license: str
    redistribution_status: str
    identifier_namespace: str
    sequence_window: SequenceWindowDefinition | None = None
    bundle_id: str | None = None
    organism_common_name: str | None = None
    supports: tuple[str, ...] = (
        "kinase_workflow",
        "site_sequence_derivation",
    )
    limitations: tuple[str, ...] = (
        "caller-supplied local source files; redistribution governed by request metadata",
    )


@dataclass(frozen=True, slots=True)
class ReferenceBundleMissingValueCount:
    """Missing-value count for one important reference field."""

    table_name: str
    field_name: str
    missing_count: int

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "table_name": self.table_name,
            "field_name": self.field_name,
            "missing_count": int(self.missing_count),
        }


@dataclass(frozen=True, slots=True)
class ReferenceBundleTableValidationReport:
    """Validation status for one required reference table."""

    table_name: str
    present: bool
    row_count: int | None
    required_columns: tuple[str, ...]
    present_columns: tuple[str, ...]
    missing_required_columns: tuple[str, ...] = ()
    missing_values: tuple[ReferenceBundleMissingValueCount, ...] = ()

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "table_name": self.table_name,
            "present": bool(self.present),
            "row_count": self.row_count,
            "required_columns": self.required_columns,
            "present_columns": self.present_columns,
            "missing_required_columns": self.missing_required_columns,
            "missing_values": tuple(item.to_payload() for item in self.missing_values),
        }


@dataclass(frozen=True, slots=True)
class ReferenceBundleSourceFileValidationReport:
    """Manifest-declared source-file status for one expected reference role."""

    role: str
    present: bool
    path: str | None = None

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "role": self.role,
            "present": bool(self.present),
            "path": self.path,
        }


@dataclass(frozen=True, slots=True)
class ReferenceBundleValidationReport:
    """Structured reference-bundle validation and provenance summary."""

    bundle_name: str | None
    bundle_version: str | None
    organism: str | None
    organism_common_name: str | None
    identifier_namespace: str | None
    required_tables: tuple[ReferenceBundleTableValidationReport, ...]
    required_source_files: tuple[ReferenceBundleSourceFileValidationReport, ...]
    kinase_substrate_record_count: int
    duplicate_record_count: int
    duplicate_records: tuple[tuple[str, str], ...] = ()
    provenance_fields: dict[str, JsonValue] = field(default_factory=_empty_json_dict)
    compatibility_warnings: tuple[str, ...] = ()

    @property
    def warnings(self) -> tuple[str, ...]:
        """Return compatibility warnings in the report."""

        return self.compatibility_warnings

    def to_payload(self) -> dict[str, JsonValue]:
        required_tables: list[JsonValue] = [
            item.to_payload() for item in self.required_tables
        ]
        required_source_files: list[JsonValue] = [
            item.to_payload() for item in self.required_source_files
        ]
        duplicate_records: list[JsonValue] = []
        for kinase, substrate_site in self.duplicate_records:
            duplicate_records.append(
                {
                    "kinase": kinase,
                    "substrate_site": substrate_site,
                }
            )
        payload: dict[str, JsonValue] = {
            "bundle_name": self.bundle_name,
            "bundle_version": self.bundle_version,
            "organism": self.organism,
            "organism_common_name": self.organism_common_name,
            "identifier_namespace": self.identifier_namespace,
            "required_tables": required_tables,
            "required_source_files": required_source_files,
            "kinase_substrate_record_count": int(self.kinase_substrate_record_count),
            "duplicate_record_count": int(self.duplicate_record_count),
            "duplicate_records": duplicate_records,
            "provenance_fields": dict(self.provenance_fields),
            "compatibility_warnings": self.compatibility_warnings,
        }
        return payload


@dataclass(frozen=True, slots=True)
class BundledReferenceLane:
    """Inventory metadata for one packaged bundled reference lane."""

    organism: Organism
    bundle_id: str
    source_name: str
    source_version: str
    retrieved_at: date
    redistribution_status: str
    supports: tuple[str, ...]
    limitations: tuple[str, ...]

    def to_payload(self) -> dict[str, JsonValue]:
        return {
            "organism": self.organism.value,
            "bundle_id": self.bundle_id,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "retrieved_at": self.retrieved_at.isoformat(),
            "redistribution_status": self.redistribution_status,
            "supports": self.supports,
            "limitations": self.limitations,
        }


@dataclass(frozen=True, slots=True)
class ReferenceBundle:
    """Resolved workflow reference resources.

    Large kinase-substrate maps are supported. Runtime in downstream workflows
    is primarily controlled by dataset/reference overlap after interpreter and
    scoring-lane filtering, not only by raw map row count.
    """

    organism: Organism
    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
    provenance: ReferenceProvenance | None = None
    manifest: ReferenceManifest | None = None
    _assume_owned: InitVar[bool] = False
    _validation_report: ReferenceBundleValidationReport = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self, _assume_owned: bool) -> None:
        if not isinstance(cast(object, self.organism), Organism):
            raise ReferenceValidationError(
                "references.organism must be an Organism enum value"
            )
        kinase_substrate_map = own_dataframe(
            self.kinase_substrate_map,
            field_name="references.kinase_substrate_map",
            error_type=ReferenceValidationError,
            assume_owned=_assume_owned,
        )
        site_sequences = own_dataframe(
            self.site_sequences,
            field_name="references.site_sequences",
            error_type=ReferenceValidationError,
            assume_owned=_assume_owned,
        )
        provenance = self.provenance
        manifest = self.manifest
        reference_context = (
            reference_context_from_manifest_if_complete(manifest)
            if manifest is not None
            else None
        )
        from phospy.validation.references.bundle import ReferenceBundleValidator

        validation = ReferenceBundleValidator().run(
            organism=self.organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
        )
        identifier_normalisation = validation.identifier_normalisation
        if provenance is None:
            provenance = ReferenceProvenance(
                source_type="explicit",
                organism=self.organism.value,
                bundle_id=None,
                table_fingerprints=(
                    fingerprint_table(
                        validation.kinase_substrate_map,
                        name="references.kinase_substrate_map",
                    ),
                    fingerprint_table(
                        validation.site_sequences,
                        name="references.site_sequences",
                    ),
                ),
                identifier_normalisation=identifier_normalisation,
                reference_context=reference_context,
            )
        else:
            resolved_reference_context = reference_context_from_provenance(
                provenance,
                manifest=manifest,
            )
            should_attach_identifier_normalisation = (
                provenance.source_type == "explicit"
                and provenance.identifier_normalisation is None
            )
            should_attach_reference_context = (
                provenance.reference_context is None
                and resolved_reference_context is not None
            )
            if (
                should_attach_identifier_normalisation
                or should_attach_reference_context
            ):
                provenance = ReferenceProvenance(
                    source_type=provenance.source_type,
                    organism=provenance.organism,
                    bundle_id=provenance.bundle_id,
                    source_name=provenance.source_name,
                    source_version=provenance.source_version,
                    retrieved_at=provenance.retrieved_at,
                    identifier_namespace=provenance.identifier_namespace,
                    sequence_window=provenance.sequence_window,
                    manifest=provenance.manifest,
                    table_fingerprints=provenance.table_fingerprints,
                    identifier_normalisation=(
                        identifier_normalisation
                        if should_attach_identifier_normalisation
                        else provenance.identifier_normalisation
                    ),
                    reference_context=(
                        resolved_reference_context
                        if should_attach_reference_context
                        else provenance.reference_context
                    ),
                )
        object.__setattr__(
            self,
            "kinase_substrate_map",
            validation.kinase_substrate_map,
        )
        object.__setattr__(
            self,
            "site_sequences",
            validation.site_sequences,
        )
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "manifest", manifest)
        object.__setattr__(self, "_validation_report", validation.report)

    @classmethod
    def _from_owned(
        cls,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
        provenance: ReferenceProvenance | None = None,
        manifest: ReferenceManifest | None = None,
    ) -> ReferenceBundle:
        return cls(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
            _assume_owned=True,
        )

    def kinase_substrate_map_dataframe(self) -> pd.DataFrame:
        """Return a kinase-substrate map snapshot isolated from this bundle."""

        return export_dataframe(self.kinase_substrate_map)

    def site_sequences_dataframe(self) -> pd.DataFrame:
        """Return a site-sequence snapshot isolated from this bundle."""

        return export_dataframe(self.site_sequences)

    @property
    def validation_report(self) -> ReferenceBundleValidationReport:
        """Return the structured validation report for this bundle."""

        return self._validation_report
