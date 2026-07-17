"""Leaf reference-context value model used by provenance and science domains."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.models import JsonValue
from phospy.provenance.organisms import Organism, normalize_organism

_REFERENCE_CONTEXT_ID_PREFIX = "reference-context-v1:"


@dataclass(frozen=True, slots=True, init=False)
class ReferenceContext:
    """Comparable biological reference identity context."""

    organism: Organism
    protein_namespace: str
    source_name: str
    source_version: str
    proteome_version: str | None
    reference_table_sha256: str | None
    reference_context_id: str = field(init=False, compare=False)

    def __init__(
        self,
        organism: object,
        protein_namespace: object,
        source_name: object,
        source_version: object,
        proteome_version: object | None,
        reference_table_sha256: object | None,
    ) -> None:
        organism_input = (
            organism
            if isinstance(organism, Organism)
            else _required_reference_context_text(
                organism,
                field_name="reference_context.organism",
            )
        )
        organism = normalize_organism(
            organism_input,
            field_name="reference_context.organism",
            error_type=ReferenceValidationError,
        )
        protein_namespace = _required_reference_context_text(
            protein_namespace,
            field_name="reference_context.protein_namespace",
        )
        source_name = _required_reference_context_text(
            source_name,
            field_name="reference_context.source_name",
        )
        source_version = _required_reference_context_text(
            source_version,
            field_name="reference_context.source_version",
        )
        proteome_version = _optional_reference_context_text(proteome_version)
        reference_table_sha256 = _optional_reference_context_text(
            reference_table_sha256
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
    def from_manifest(cls, manifest: object) -> ReferenceContext:
        """Build a reference context from manifest identity metadata."""

        source_version = getattr(manifest, "source_version", None)
        if source_version is None:
            raise ReferenceValidationError(
                "reference_context.source_version must be non-empty"
            )
        organism = _required_reference_context_text(
            getattr(manifest, "organism_common_name", None)
            or getattr(manifest, "organism", None),
            field_name="reference_context.organism",
        )
        protein_namespace = _required_reference_context_text(
            getattr(manifest, "protein_namespace", None),
            field_name="reference_context.protein_namespace",
        )
        source_name = _required_reference_context_text(
            getattr(manifest, "source_name", None),
            field_name="reference_context.source_name",
        )
        return cls(
            organism=organism,
            protein_namespace=protein_namespace,
            source_name=source_name,
            source_version=_required_reference_context_text(
                source_version,
                field_name="reference_context.source_version",
            ),
            proteome_version=None,
            reference_table_sha256=getattr(manifest, "table_sha256", None),
        )

    @classmethod
    def from_manifest_payload(
        cls,
        payload: Mapping[str, JsonValue],
    ) -> ReferenceContext | None:
        """Build a reference context from serialized manifest identity metadata."""

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
        return cls(
            organism=organism,
            protein_namespace=protein_namespace,
            source_name=source_name,
            source_version=source_version,
            proteome_version=_payload_optional_text(payload.get("proteome_version")),
            reference_table_sha256=_payload_optional_text(payload.get("table_sha256")),
        )

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> ReferenceContext:
        """Deserialize a reference context payload."""

        payload = _require_payload_mapping(payload, field_name="reference_context")
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
            "organism": self.organism.value,
            "protein_namespace": self.protein_namespace,
            "source_name": self.source_name,
            "source_version": self.source_version,
            "proteome_version": self.proteome_version,
            "reference_table_sha256": self.reference_table_sha256,
        }


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


def _require_payload_mapping(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> Mapping[str, object]:
    result: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str):
            raise ReferenceValidationError(
                f"{field_name} JSON object keys must be strings; "
                f"got {type(key).__name__}"
            )
        if key in result:
            raise ReferenceValidationError(
                f"{field_name} contains duplicate JSON object key {key!r}"
            )
        result[key] = item
    return result


__all__ = ["ReferenceContext"]
