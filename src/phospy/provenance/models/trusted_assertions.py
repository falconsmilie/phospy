"""Trusted direct dataset-construction assertion provenance models."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import cast

from phospy.errors.input import PhosPyInputError
from phospy.provenance.immutability import freeze_json_mapping, thaw_json_mapping
from phospy.provenance.models._shared import (
    JsonValue,
    _optional_provenance_float,
    _optional_provenance_text,
    _required_non_negative_row_count,
    _required_provenance_bool,
    _required_provenance_text,
)

TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V1 = 1

TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V2 = 2

TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V3 = 3

TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V4 = 4

_TRUSTED_DATASET_CONSTRUCTION_REQUIRED_DIMENSIONS = (
    "identity",
    "intensity_scale",
    "quantitative_meaning",
    "aligned_structure",
    "localisation",
    "sequence",
    "reference_context",
)

_TRUSTED_DATASET_CONSTRUCTION_OPTIONAL_DIMENSIONS = ("numeric_semantic_domain",)

_TRUSTED_DATASET_CONSTRUCTION_ASSERTION_DIMENSIONS = (
    _TRUSTED_DATASET_CONSTRUCTION_REQUIRED_DIMENSIONS
    + _TRUSTED_DATASET_CONSTRUCTION_OPTIONAL_DIMENSIONS
)

_TRUSTED_DATASET_CONSTRUCTION_MISSING_ASSERTION_NAMES = {
    "identity": "identity_user_asserted",
    "intensity_scale": "intensity_scale_user_asserted",
    "quantitative_meaning": "quantitative_meaning_user_asserted",
    "aligned_structure": "aligned_structure_user_asserted",
    "localisation": "localisation_user_asserted",
    "sequence": "sequence_user_asserted",
    "reference_context": "reference_context_user_asserted",
}

_TRUSTED_DATASET_CONSTRUCTION_EVIDENCE_KINDS = frozenset({"evidence", "waiver"})

_TRUSTED_DATASET_CONSTRUCTION_EVIDENCE_PAYLOAD_KEYS = frozenset(
    {
        "kind",
        "source",
        "policy",
        "threshold",
        "waiver_reason",
        "details",
    }
)

_TRUSTED_DATASET_CONSTRUCTION_ASSERTION_PAYLOAD_KEYS = frozenset(
    {
        "schema_version",
        "assertion_metadata_provided",
        "identity",
        "intensity_scale",
        "quantitative_meaning",
        "aligned_structure",
        "localisation",
        "sequence",
        "reference_context",
        "numeric_semantic_domain",
        "identity_user_asserted",
        "intensity_scale_user_asserted",
        "quantitative_meaning_user_asserted",
        "aligned_structure_user_asserted",
        "localisation_user_asserted",
        "sequence_user_asserted",
        "reference_context_user_asserted",
        "numeric_semantic_domain_user_asserted",
        "waived_assertions",
        "missing_assertions",
        "asserted_by",
        "assertion_source",
        "notes",
        "assertion_fingerprint",
    }
)


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

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str = "trusted_dataset_construction_evidence",
    ) -> TrustedDatasetConstructionEvidence:
        """Deserialize a canonical trusted construction evidence payload."""

        payload = _required_trusted_payload_mapping(payload, field_name=field_name)
        _require_exact_payload_keys(
            payload,
            expected_keys=_TRUSTED_DATASET_CONSTRUCTION_EVIDENCE_PAYLOAD_KEYS,
            field_name=field_name,
        )
        details = cast(
            Mapping[str, JsonValue],
            _required_trusted_payload_mapping(
                payload.get("details"),
                field_name=f"{field_name}.details",
            ),
        )
        evidence = cls(
            kind=_required_trusted_payload_text(
                payload.get("kind"),
                field_name=f"{field_name}.kind",
            ),
            source=_optional_provenance_text(payload.get("source")),
            policy=_optional_provenance_text(payload.get("policy")),
            threshold=_optional_provenance_float(
                payload.get("threshold"),
                field_name=f"{field_name}.threshold",
            ),
            waiver_reason=_optional_provenance_text(payload.get("waiver_reason")),
            details=details,
        )
        _require_canonical_payload(
            observed=payload,
            expected=evidence.to_payload(),
            field_name=field_name,
        )
        return evidence

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
    explicit waiver for identity, intensity scale, quantitative meaning,
    aligned table structure, localisation, sequence, and reference context.
    Numeric-semantic domain evidence is optional unless a trusted construction
    needs to make an explicit visible waiver for a scale/meaning/value-domain
    conflict.
    A missing assertion bundle is reserved for legacy direct-construction audit
    markers.
    """

    identity: TrustedDatasetConstructionEvidence | None = None
    intensity_scale: TrustedDatasetConstructionEvidence | None = None
    quantitative_meaning: TrustedDatasetConstructionEvidence | None = None
    aligned_structure: TrustedDatasetConstructionEvidence | None = None
    localisation: TrustedDatasetConstructionEvidence | None = None
    sequence: TrustedDatasetConstructionEvidence | None = None
    reference_context: TrustedDatasetConstructionEvidence | None = None
    numeric_semantic_domain: TrustedDatasetConstructionEvidence | None = None
    assertion_metadata_provided: bool = True
    asserted_by: str | None = None
    assertion_source: str | None = None
    notes: str | None = None
    schema_version: int = TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V4

    def __post_init__(self) -> None:
        identity = _optional_trusted_construction_evidence(
            self.identity,
            field_name="trusted_dataset_construction_assertions.identity",
        )
        intensity_scale = _optional_trusted_construction_evidence(
            self.intensity_scale,
            field_name="trusted_dataset_construction_assertions.intensity_scale",
        )
        quantitative_meaning = _optional_trusted_construction_evidence(
            self.quantitative_meaning,
            field_name="trusted_dataset_construction_assertions.quantitative_meaning",
        )
        aligned_structure = _optional_trusted_construction_evidence(
            self.aligned_structure,
            field_name="trusted_dataset_construction_assertions.aligned_structure",
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
        numeric_semantic_domain = _optional_trusted_construction_evidence(
            self.numeric_semantic_domain,
            field_name=(
                "trusted_dataset_construction_assertions.numeric_semantic_domain"
            ),
        )
        assertion_metadata_provided = _required_provenance_bool(
            self.assertion_metadata_provided,
            field_name=(
                "trusted_dataset_construction_assertions.assertion_metadata_provided"
            ),
        )
        supplied_assertions = {
            "identity": identity,
            "intensity_scale": intensity_scale,
            "quantitative_meaning": quantitative_meaning,
            "aligned_structure": aligned_structure,
            "localisation": localisation,
            "sequence": sequence,
            "reference_context": reference_context,
        }
        optional_assertions = {
            "numeric_semantic_domain": numeric_semantic_domain,
        }
        if not assertion_metadata_provided and any(
            value is not None
            for value in (
                *supplied_assertions.values(),
                *optional_assertions.values(),
            )
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
            _require_aligned_structure_evidence(aligned_structure)
            _require_localisation_evidence_or_waiver(localisation)
        schema_version = _required_non_negative_row_count(
            self.schema_version,
            field_name="trusted_dataset_construction_assertions.schema_version",
        )
        if schema_version != TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V4:
            raise PhosPyInputError(
                "trusted_dataset_construction_assertions.schema_version must be "
                f"{TRUSTED_DATASET_CONSTRUCTION_ASSERTIONS_SCHEMA_VERSION_V4}"
            )
        object.__setattr__(self, "identity", identity)
        object.__setattr__(self, "intensity_scale", intensity_scale)
        object.__setattr__(self, "quantitative_meaning", quantitative_meaning)
        object.__setattr__(self, "aligned_structure", aligned_structure)
        object.__setattr__(self, "localisation", localisation)
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "reference_context", reference_context)
        object.__setattr__(
            self,
            "numeric_semantic_domain",
            numeric_semantic_domain,
        )
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
    def intensity_scale_user_asserted(self) -> bool:
        """Return whether intensity scale has typed evidence or a waiver."""

        return self.intensity_scale is not None

    @property
    def quantitative_meaning_user_asserted(self) -> bool:
        """Return whether quantitative meaning has typed evidence or a waiver."""

        return self.quantitative_meaning is not None

    @property
    def aligned_structure_user_asserted(self) -> bool:
        """Return whether aligned table structure has typed evidence or a waiver."""

        return self.aligned_structure is not None

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
    def numeric_semantic_domain_user_asserted(self) -> bool:
        """Return whether numeric-semantic domain has evidence or a waiver."""

        return self.numeric_semantic_domain is not None

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
        for dimension in _TRUSTED_DATASET_CONSTRUCTION_ASSERTION_DIMENSIONS:
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

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
        *,
        field_name: str = "trusted_dataset_construction_assertions",
    ) -> TrustedDatasetConstructionAssertions:
        """Deserialize a canonical trusted construction assertion payload."""

        payload = _required_trusted_payload_mapping(payload, field_name=field_name)
        _require_exact_payload_keys(
            payload,
            expected_keys=_TRUSTED_DATASET_CONSTRUCTION_ASSERTION_PAYLOAD_KEYS,
            field_name=field_name,
        )
        assertions = cls(
            identity=_trusted_evidence_from_payload(
                payload.get("identity"),
                field_name=f"{field_name}.identity",
            ),
            intensity_scale=_trusted_evidence_from_payload(
                payload.get("intensity_scale"),
                field_name=f"{field_name}.intensity_scale",
            ),
            quantitative_meaning=_trusted_evidence_from_payload(
                payload.get("quantitative_meaning"),
                field_name=f"{field_name}.quantitative_meaning",
            ),
            aligned_structure=_trusted_evidence_from_payload(
                payload.get("aligned_structure"),
                field_name=f"{field_name}.aligned_structure",
            ),
            localisation=_trusted_evidence_from_payload(
                payload.get("localisation"),
                field_name=f"{field_name}.localisation",
            ),
            sequence=_trusted_evidence_from_payload(
                payload.get("sequence"),
                field_name=f"{field_name}.sequence",
            ),
            reference_context=_trusted_evidence_from_payload(
                payload.get("reference_context"),
                field_name=f"{field_name}.reference_context",
            ),
            numeric_semantic_domain=_trusted_evidence_from_payload(
                payload.get("numeric_semantic_domain"),
                field_name=f"{field_name}.numeric_semantic_domain",
            ),
            assertion_metadata_provided=_required_provenance_bool(
                payload.get("assertion_metadata_provided"),
                field_name=f"{field_name}.assertion_metadata_provided",
            ),
            asserted_by=_optional_provenance_text(payload.get("asserted_by")),
            assertion_source=_optional_provenance_text(payload.get("assertion_source")),
            notes=_optional_provenance_text(payload.get("notes")),
            schema_version=_required_non_negative_row_count(
                payload.get("schema_version"),
                field_name=f"{field_name}.schema_version",
            ),
        )
        _require_canonical_payload(
            observed=payload,
            expected=assertions.to_payload(),
            field_name=field_name,
        )
        return assertions

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
            "intensity_scale": _trusted_evidence_payload(self.intensity_scale),
            "quantitative_meaning": _trusted_evidence_payload(
                self.quantitative_meaning
            ),
            "aligned_structure": _trusted_evidence_payload(self.aligned_structure),
            "localisation": _trusted_evidence_payload(self.localisation),
            "sequence": _trusted_evidence_payload(self.sequence),
            "reference_context": _trusted_evidence_payload(self.reference_context),
            "numeric_semantic_domain": _trusted_evidence_payload(
                self.numeric_semantic_domain
            ),
            "identity_user_asserted": bool(self.identity_user_asserted),
            "intensity_scale_user_asserted": bool(self.intensity_scale_user_asserted),
            "quantitative_meaning_user_asserted": bool(
                self.quantitative_meaning_user_asserted
            ),
            "aligned_structure_user_asserted": bool(
                self.aligned_structure_user_asserted
            ),
            "localisation_user_asserted": bool(self.localisation_user_asserted),
            "sequence_user_asserted": bool(self.sequence_user_asserted),
            "reference_context_user_asserted": bool(
                self.reference_context_user_asserted
            ),
            "numeric_semantic_domain_user_asserted": bool(
                self.numeric_semantic_domain_user_asserted
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


def _trusted_evidence_from_payload(
    value: object | None,
    *,
    field_name: str,
) -> TrustedDatasetConstructionEvidence | None:
    if value is None:
        return None
    payload = _required_trusted_payload_mapping(value, field_name=field_name)
    return TrustedDatasetConstructionEvidence.from_payload(
        payload,
        field_name=field_name,
    )


def _required_trusted_payload_mapping(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise PhosPyInputError(f"{field_name} must be a mapping")
    return thaw_json_mapping(value, field_name=field_name)


def _required_trusted_payload_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must be a string")
    return _required_provenance_text(value, field_name=field_name)


def _require_exact_payload_keys(
    payload: Mapping[str, object],
    *,
    expected_keys: frozenset[str],
    field_name: str,
) -> None:
    actual_keys = set(payload)
    missing = sorted(expected_keys - actual_keys)
    unexpected = sorted(actual_keys - expected_keys)
    if not missing and not unexpected:
        return
    details: list[str] = []
    if missing:
        details.append("missing keys: " + ", ".join(missing))
    if unexpected:
        details.append("unexpected keys: " + ", ".join(unexpected))
    raise PhosPyInputError(f"{field_name} is not canonical; " + "; ".join(details))


def _require_canonical_payload(
    *,
    observed: Mapping[str, object],
    expected: Mapping[str, object],
    field_name: str,
) -> None:
    if _canonical_json_payload(observed) == _canonical_json_payload(expected):
        return
    raise PhosPyInputError(
        f"{field_name} is not canonical; assertion_fingerprint or derived "
        "assertion fields do not match the recomputed canonical payload"
    )


def _canonical_json_payload(payload: Mapping[str, object]) -> str:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


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


def _require_aligned_structure_evidence(
    value: TrustedDatasetConstructionEvidence | None,
) -> None:
    if value is None:
        return
    if value.kind == "waiver":
        raise PhosPyInputError(
            "trusted_dataset_construction_assertions.aligned_structure cannot be "
            "waived; table shape and alignment are mechanically verified"
        )
