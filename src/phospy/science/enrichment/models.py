"""Domain contracts for enrichment workflows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, cast

from phospy.errors.validation import WorkflowValidationError

ENRICHMENT_METHOD_OVER_REPRESENTATION = "over_representation"
EnrichmentMethod = Literal["over_representation"]
SUPPORTED_ENRICHMENT_METHODS: tuple[EnrichmentMethod, ...] = (
    ENRICHMENT_METHOD_OVER_REPRESENTATION,
)

MULTIPLE_TESTING_CORRECTION_NONE = "none"
MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG = "benjamini_hochberg"
MultipleTestingCorrection = Literal["none", "benjamini_hochberg"]
SUPPORTED_MULTIPLE_TESTING_CORRECTIONS: tuple[MultipleTestingCorrection, ...] = (
    MULTIPLE_TESTING_CORRECTION_NONE,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
)

ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL = "gene_symbol"
ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID = "protein_id"
ENRICHMENT_IDENTIFIER_KIND_SITE_KEY = "site_key"
ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID = "display_id"
ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE = "phosphosite"
EnrichmentIdentifierKind = Literal[
    "gene_symbol",
    "protein_id",
    "site_key",
    "display_id",
    "phosphosite",
]
SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS: tuple[EnrichmentIdentifierKind, ...] = (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID,
    ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID,
    ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE,
)
GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS: frozenset[EnrichmentIdentifierKind] = frozenset(
    {
        ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID,
    }
)
PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS: frozenset[EnrichmentIdentifierKind] = frozenset(
    {
        ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID,
        ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE,
    }
)

ENRICHMENT_COLLECTION_KIND_GENE_SET = "gene_set"
ENRICHMENT_COLLECTION_KIND_PTM_SET = "ptm_set"
EnrichmentCollectionKind = Literal["gene_set", "ptm_set"]
SUPPORTED_ENRICHMENT_COLLECTION_KINDS: tuple[EnrichmentCollectionKind, ...] = (
    ENRICHMENT_COLLECTION_KIND_GENE_SET,
    ENRICHMENT_COLLECTION_KIND_PTM_SET,
)


@dataclass(frozen=True, slots=True)
class EnrichmentSet:
    """One named enrichment set with explicit identifier semantics."""

    set_id: str
    name: str
    identifiers: Sequence[str]
    identifier_kind: EnrichmentIdentifierKind
    source_name: str | None = None
    source_version: str | None = None
    description: str | None = None

    def __post_init__(self) -> None:
        set_id = _require_non_empty_string(
            self.set_id,
            field_name="enrichment_set.set_id",
        )
        name = _require_non_empty_string(
            self.name,
            field_name=f"enrichment_set[{set_id!r}].name",
        )
        identifiers = _normalise_member_identifier_sequence(
            self.identifiers,
            field_name=f"enrichment_set[{set_id!r}].identifiers",
        )
        identifier_kind = _require_identifier_kind(
            self.identifier_kind,
            field_name=f"enrichment_set[{set_id!r}].identifier_kind",
        )
        object.__setattr__(self, "set_id", set_id)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "identifiers", identifiers)
        object.__setattr__(self, "identifier_kind", identifier_kind)
        object.__setattr__(
            self,
            "source_name",
            _normalise_optional_string(
                self.source_name,
                field_name=f"enrichment_set[{set_id!r}].source_name",
            ),
        )
        object.__setattr__(
            self,
            "source_version",
            _normalise_optional_string(
                self.source_version,
                field_name=f"enrichment_set[{set_id!r}].source_version",
            ),
        )
        object.__setattr__(
            self,
            "description",
            _normalise_optional_string(
                self.description,
                field_name=f"enrichment_set[{set_id!r}].description",
            ),
        )


@dataclass(frozen=True, slots=True, init=False)
class EnrichmentSetCollection:
    """Offline enrichment set collection with homogeneous identifier semantics."""

    enrichment_sets: tuple[EnrichmentSet, ...]
    identifier_kind: EnrichmentIdentifierKind
    collection_kind: EnrichmentCollectionKind
    source_name: str | None
    source_version: str | None

    def __init__(
        self,
        *,
        sets: Sequence[EnrichmentSet],
        identifier_kind: EnrichmentIdentifierKind | None = None,
        collection_kind: EnrichmentCollectionKind | None = None,
        source_name: str | None = None,
        source_version: str | None = None,
    ) -> None:
        enrichment_sets = _normalise_enrichment_set_sequence(
            sets,
            field_name="enrichment_set_collection.sets",
        )
        resolved_identifier_kind = _resolve_collection_identifier_kind(
            enrichment_sets=enrichment_sets,
            identifier_kind=identifier_kind,
            field_name="enrichment_set_collection.identifier_kind",
        )
        resolved_collection_kind = _resolve_collection_kind(
            collection_kind=collection_kind,
            identifier_kind=resolved_identifier_kind,
            field_name="enrichment_set_collection.collection_kind",
        )
        object.__setattr__(self, "enrichment_sets", enrichment_sets)
        object.__setattr__(self, "identifier_kind", resolved_identifier_kind)
        object.__setattr__(self, "collection_kind", resolved_collection_kind)
        object.__setattr__(
            self,
            "source_name",
            _normalise_optional_string(
                source_name,
                field_name="enrichment_set_collection.source_name",
            ),
        )
        object.__setattr__(
            self,
            "source_version",
            _normalise_optional_string(
                source_version,
                field_name="enrichment_set_collection.source_version",
            ),
        )

    @property
    def sets(self) -> tuple[EnrichmentSet, ...] | dict[str, tuple[str, ...]]:
        return self.enrichment_sets

    @property
    def set_ids(self) -> tuple[str, ...]:
        return tuple(enrichment_set.set_id for enrichment_set in self.enrichment_sets)

    @property
    def set_by_id(self) -> dict[str, EnrichmentSet]:
        return {
            enrichment_set.set_id: enrichment_set
            for enrichment_set in self.enrichment_sets
        }

    @property
    def term_names(self) -> dict[str, str]:
        return {
            enrichment_set.set_id: enrichment_set.name
            for enrichment_set in self.enrichment_sets
        }

    @property
    def members_by_set_id(self) -> dict[str, tuple[str, ...]]:
        return {
            enrichment_set.set_id: tuple(enrichment_set.identifiers)
            for enrichment_set in self.enrichment_sets
        }


@dataclass(frozen=True, slots=True, init=False)
class GeneSetCollection(EnrichmentSetCollection):
    """Explicit gene-level enrichment collection.

    The identifiers inside ``sets`` are gene/protein-level identifiers. This
    object is intentionally separate from ``PtmSetCollection`` so a request
    cannot silently mix gene-set and phosphosite/PTM-set semantics.
    """

    def __init__(
        self,
        *,
        sets: Mapping[str, Sequence[str]] | Sequence[EnrichmentSet],
        identifier_kind: EnrichmentIdentifierKind = (
            ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
        ),
        term_names: Mapping[str, str] | None = None,
        source_name: str | None = "user",
        source_version: str | None = None,
        descriptions: Mapping[str, str | None] | None = None,
    ) -> None:
        identifier_kind = _require_identifier_kind(
            identifier_kind,
            field_name="gene_set_collection.identifier_kind",
        )
        if identifier_kind not in GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
            allowed = ", ".join(sorted(GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS))
            raise WorkflowValidationError(
                f"gene_set_collection.identifier_kind must be gene-level: {allowed}"
            )
        enrichment_sets = _coerce_legacy_collection_sets(
            sets,
            identifier_kind=identifier_kind,
            term_names=term_names,
            descriptions=descriptions,
            source_name=source_name,
            source_version=source_version,
            field_name="gene_set_collection.sets",
        )
        EnrichmentSetCollection.__init__(
            self,
            sets=enrichment_sets,
            identifier_kind=identifier_kind,
            collection_kind=ENRICHMENT_COLLECTION_KIND_GENE_SET,
            source_name=source_name,
            source_version=source_version,
        )

    @property
    def sets(self) -> dict[str, tuple[str, ...]]:
        return self.members_by_set_id


@dataclass(frozen=True, slots=True, init=False)
class PtmSetCollection(EnrichmentSetCollection):
    """Explicit phosphosite/PTM-level enrichment collection."""

    def __init__(
        self,
        *,
        sets: Mapping[str, Sequence[str]] | Sequence[EnrichmentSet],
        identifier_kind: EnrichmentIdentifierKind = ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        term_names: Mapping[str, str] | None = None,
        source_name: str | None = "user",
        source_version: str | None = None,
        descriptions: Mapping[str, str | None] | None = None,
    ) -> None:
        identifier_kind = _require_identifier_kind(
            identifier_kind,
            field_name="ptm_set_collection.identifier_kind",
        )
        if identifier_kind not in PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
            allowed = ", ".join(sorted(PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS))
            raise WorkflowValidationError(
                f"ptm_set_collection.identifier_kind must be PTM-level: {allowed}"
            )
        enrichment_sets = _coerce_legacy_collection_sets(
            sets,
            identifier_kind=identifier_kind,
            term_names=term_names,
            descriptions=descriptions,
            source_name=source_name,
            source_version=source_version,
            field_name="ptm_set_collection.sets",
        )
        EnrichmentSetCollection.__init__(
            self,
            sets=enrichment_sets,
            identifier_kind=identifier_kind,
            collection_kind=ENRICHMENT_COLLECTION_KIND_PTM_SET,
            source_name=source_name,
            source_version=source_version,
        )

    @property
    def sets(self) -> dict[str, tuple[str, ...]]:
        return self.members_by_set_id


@dataclass(frozen=True, slots=True)
class EnrichmentResultRecord:
    """One enrichment term row.

    This is a result-shape contract only. It stores counts and optional
    p-values produced by a future executor; it does not calculate them.
    """

    term_id: str
    collection_kind: EnrichmentCollectionKind
    identifier_kind: EnrichmentIdentifierKind
    input_overlap_count: int
    background_overlap_count: int
    set_size: int
    term_name: str | None = None
    overlap_identifiers: tuple[str, ...] = ()
    p_value: float | None = None
    adjusted_p_value: float | None = None

    def __post_init__(self) -> None:
        term_id = _require_non_empty_string(
            self.term_id,
            field_name="enrichment_result_record.term_id",
        )
        collection_kind = _require_supported_value(
            self.collection_kind,
            field_name="enrichment_result_record.collection_kind",
            supported=SUPPORTED_ENRICHMENT_COLLECTION_KINDS,
        )
        identifier_kind = _require_identifier_kind(
            self.identifier_kind,
            field_name="enrichment_result_record.identifier_kind",
        )
        object.__setattr__(self, "term_id", term_id)
        object.__setattr__(
            self,
            "collection_kind",
            cast(EnrichmentCollectionKind, collection_kind),
        )
        object.__setattr__(self, "identifier_kind", identifier_kind)
        object.__setattr__(
            self,
            "input_overlap_count",
            _require_int_at_least(
                self.input_overlap_count,
                field_name="enrichment_result_record.input_overlap_count",
                minimum=0,
            ),
        )
        object.__setattr__(
            self,
            "background_overlap_count",
            _require_int_at_least(
                self.background_overlap_count,
                field_name="enrichment_result_record.background_overlap_count",
                minimum=0,
            ),
        )
        object.__setattr__(
            self,
            "set_size",
            _require_int_at_least(
                self.set_size,
                field_name="enrichment_result_record.set_size",
                minimum=0,
            ),
        )
        object.__setattr__(
            self,
            "term_name",
            (
                term_id
                if self.term_name is None
                else _require_non_empty_string(
                    self.term_name,
                    field_name="enrichment_result_record.term_name",
                )
            ),
        )
        object.__setattr__(
            self,
            "overlap_identifiers",
            _normalise_identifier_sequence(
                self.overlap_identifiers,
                field_name="enrichment_result_record.overlap_identifiers",
                allow_empty=True,
            ),
        )
        object.__setattr__(
            self,
            "p_value",
            _normalise_optional_unit_interval(
                self.p_value,
                field_name="enrichment_result_record.p_value",
            ),
        )
        object.__setattr__(
            self,
            "adjusted_p_value",
            _normalise_optional_unit_interval(
                self.adjusted_p_value,
                field_name="enrichment_result_record.adjusted_p_value",
            ),
        )


def _require_identifier_kind(
    value: object,
    *,
    field_name: str,
) -> EnrichmentIdentifierKind:
    return cast(
        EnrichmentIdentifierKind,
        _require_supported_value(
            value,
            field_name=field_name,
            supported=SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS,
        ),
    )


def _normalise_identifier_sequence(
    value: object,
    *,
    field_name: str,
    allow_empty: bool = False,
) -> tuple[str, ...]:
    if isinstance(value, str | bytes | bytearray) or not isinstance(value, Sequence):
        raise WorkflowValidationError(f"{field_name} must be a sequence of strings")
    normalised = tuple(
        _require_non_empty_string(item, field_name=f"{field_name}[]") for item in value
    )
    if not allow_empty and not normalised:
        raise WorkflowValidationError(f"{field_name} must not be empty")
    return normalised


def _normalise_member_identifier_sequence(
    value: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    identifiers = _normalise_identifier_sequence(value, field_name=field_name)
    seen: set[str] = set()
    distinct: list[str] = []
    for identifier in identifiers:
        if identifier in seen:
            continue
        seen.add(identifier)
        distinct.append(identifier)
    return tuple(distinct)


def _normalise_enrichment_set_sequence(
    value: object,
    *,
    field_name: str,
) -> tuple[EnrichmentSet, ...]:
    if isinstance(value, str | bytes | bytearray) or not isinstance(value, Sequence):
        raise WorkflowValidationError(f"{field_name} must be a sequence")
    enrichment_sets = tuple(value)
    if not enrichment_sets:
        raise WorkflowValidationError(f"{field_name} must not be empty")
    normalised: list[EnrichmentSet] = []
    observed_ids: set[str] = set()
    for enrichment_set in enrichment_sets:
        if not isinstance(enrichment_set, EnrichmentSet):
            raise WorkflowValidationError(
                f"{field_name} must contain EnrichmentSet values"
            )
        if enrichment_set.set_id in observed_ids:
            raise WorkflowValidationError(f"{field_name} set_id values must be unique")
        observed_ids.add(enrichment_set.set_id)
        normalised.append(enrichment_set)
    return tuple(normalised)


def _resolve_collection_identifier_kind(
    *,
    enrichment_sets: tuple[EnrichmentSet, ...],
    identifier_kind: EnrichmentIdentifierKind | None,
    field_name: str,
) -> EnrichmentIdentifierKind:
    resolved_identifier_kind = (
        enrichment_sets[0].identifier_kind
        if identifier_kind is None
        else _require_identifier_kind(identifier_kind, field_name=field_name)
    )
    for enrichment_set in enrichment_sets:
        if enrichment_set.identifier_kind != resolved_identifier_kind:
            raise WorkflowValidationError(
                f"{field_name} cannot mix identifier_kind values; "
                f"set_id={enrichment_set.set_id!r}, "
                f"observed={enrichment_set.identifier_kind!r}, "
                f"expected={resolved_identifier_kind!r}"
            )
    return resolved_identifier_kind


def _resolve_collection_kind(
    *,
    collection_kind: EnrichmentCollectionKind | None,
    identifier_kind: EnrichmentIdentifierKind,
    field_name: str,
) -> EnrichmentCollectionKind:
    expected_collection_kind = _collection_kind_for_identifier_kind(identifier_kind)
    if collection_kind is None:
        return expected_collection_kind
    resolved_collection_kind = cast(
        EnrichmentCollectionKind,
        _require_supported_value(
            collection_kind,
            field_name=field_name,
            supported=SUPPORTED_ENRICHMENT_COLLECTION_KINDS,
        ),
    )
    if resolved_collection_kind != expected_collection_kind:
        raise WorkflowValidationError(
            f"{field_name} must match identifier_kind {identifier_kind!r}; "
            f"observed={resolved_collection_kind!r}, "
            f"expected={expected_collection_kind!r}"
        )
    return resolved_collection_kind


def _collection_kind_for_identifier_kind(
    identifier_kind: EnrichmentIdentifierKind,
) -> EnrichmentCollectionKind:
    if identifier_kind in GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
        return ENRICHMENT_COLLECTION_KIND_GENE_SET
    if identifier_kind in PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
        return ENRICHMENT_COLLECTION_KIND_PTM_SET
    raise WorkflowValidationError(
        f"enrichment identifier_kind {identifier_kind!r} is not supported"
    )


def _coerce_legacy_collection_sets(
    value: object,
    *,
    identifier_kind: EnrichmentIdentifierKind,
    term_names: Mapping[str, str] | None,
    descriptions: Mapping[str, str | None] | None,
    source_name: str | None,
    source_version: str | None,
    field_name: str,
) -> tuple[EnrichmentSet, ...]:
    if isinstance(value, Mapping):
        sets = _normalise_set_mapping(value, field_name=field_name)
        normalised_term_names = _normalise_term_names(
            term_names or {},
            allowed_term_ids=frozenset(sets),
            field_name=field_name.replace(".sets", ".term_names"),
        )
        normalised_descriptions = _normalise_optional_term_descriptions(
            descriptions or {},
            allowed_term_ids=frozenset(sets),
            field_name=field_name.replace(".sets", ".descriptions"),
        )
        return tuple(
            EnrichmentSet(
                set_id=set_id,
                name=normalised_term_names.get(set_id, set_id),
                identifiers=identifiers,
                identifier_kind=identifier_kind,
                source_name=source_name,
                source_version=source_version,
                description=normalised_descriptions.get(set_id),
            )
            for set_id, identifiers in sets.items()
        )
    if term_names is not None:
        raise WorkflowValidationError(
            f"{field_name} term_names are only valid when sets is a mapping"
        )
    if descriptions is not None:
        raise WorkflowValidationError(
            f"{field_name} descriptions are only valid when sets is a mapping"
        )
    return _normalise_enrichment_set_sequence(value, field_name=field_name)


def _require_collection_matches_identifier_kind(
    collection: object,
    *,
    identifier_kind: EnrichmentIdentifierKind,
    field_name: str,
) -> EnrichmentSetCollection:
    if not isinstance(collection, EnrichmentSetCollection):
        raise WorkflowValidationError(
            f"{field_name} must be EnrichmentSetCollection, "
            "GeneSetCollection, or PtmSetCollection"
        )
    if collection.identifier_kind != identifier_kind:
        raise WorkflowValidationError(
            f"{field_name}.identifier_kind must match request identifier_kind; "
            f"observed={collection.identifier_kind!r}, expected={identifier_kind!r}"
        )
    expected_collection_kind = _collection_kind_for_identifier_kind(identifier_kind)
    if collection.collection_kind != expected_collection_kind:
        raise WorkflowValidationError(
            f"{field_name}.collection_kind must match identifier_kind; "
            f"observed={collection.collection_kind!r}, "
            f"expected={expected_collection_kind!r}"
        )
    return collection


def _require_supported_value(
    value: object,
    *,
    field_name: str,
    supported: tuple[str, ...],
) -> str:
    if not isinstance(value, str) or value not in supported:
        allowed = ", ".join(repr(item) for item in supported)
        raise WorkflowValidationError(f"{field_name} must be one of: {allowed}")
    return value


def _normalise_set_mapping(
    value: object,
    *,
    field_name: str,
) -> dict[str, tuple[str, ...]]:
    if not isinstance(value, Mapping):
        raise WorkflowValidationError(f"{field_name} must be a mapping")
    if not value:
        raise WorkflowValidationError(f"{field_name} must not be empty")
    normalised: dict[str, tuple[str, ...]] = {}
    for raw_term_id, raw_identifiers in value.items():
        term_id = _require_non_empty_string(
            raw_term_id,
            field_name=f"{field_name}.term_id",
        )
        if term_id in normalised:
            raise WorkflowValidationError(f"{field_name} term IDs must be unique")
        identifiers = _normalise_identifier_sequence(
            raw_identifiers,
            field_name=f"{field_name}[{term_id!r}]",
        )
        normalised[term_id] = identifiers
    return normalised


def _normalise_term_names(
    value: object,
    *,
    allowed_term_ids: frozenset[str],
    field_name: str,
) -> dict[str, str]:
    if not isinstance(value, Mapping):
        raise WorkflowValidationError(f"{field_name} must be a mapping")
    normalised: dict[str, str] = {}
    for raw_term_id, raw_term_name in value.items():
        term_id = _require_non_empty_string(
            raw_term_id,
            field_name=f"{field_name}.term_id",
        )
        if term_id not in allowed_term_ids:
            raise WorkflowValidationError(
                f"{field_name} keys must refer to IDs present in sets"
            )
        normalised[term_id] = _require_non_empty_string(
            raw_term_name,
            field_name=f"{field_name}[{term_id!r}]",
        )
    return normalised


def _normalise_optional_term_descriptions(
    value: object,
    *,
    allowed_term_ids: frozenset[str],
    field_name: str,
) -> dict[str, str]:
    if not isinstance(value, Mapping):
        raise WorkflowValidationError(f"{field_name} must be a mapping")
    normalised: dict[str, str] = {}
    for raw_term_id, raw_description in value.items():
        term_id = _require_non_empty_string(
            raw_term_id,
            field_name=f"{field_name}.term_id",
        )
        if term_id not in allowed_term_ids:
            raise WorkflowValidationError(
                f"{field_name} keys must refer to IDs present in sets"
            )
        description = _normalise_optional_string(
            raw_description,
            field_name=f"{field_name}[{term_id!r}]",
        )
        if description is not None:
            normalised[term_id] = description
    return normalised


def _require_non_empty_string(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise WorkflowValidationError(f"{field_name} must be a non-empty string")
    return value.strip()


def _normalise_optional_string(value: object | None, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_non_empty_string(value, field_name=field_name)


def _require_int_at_least(value: object, *, field_name: str, minimum: int) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise WorkflowValidationError(f"{field_name} must be an integer >= {minimum}")
    return int(value)


def _normalise_optional_unit_interval(
    value: object | None,
    *,
    field_name: str,
) -> float | None:
    if value is None:
        return None
    if not isinstance(value, int | float) or isinstance(value, bool):
        raise WorkflowValidationError(f"{field_name} must be numeric or None")
    normalised = float(value)
    if normalised < 0.0 or normalised > 1.0:
        raise WorkflowValidationError(f"{field_name} must be within [0.0, 1.0]")
    return normalised


__all__ = [
    "ENRICHMENT_COLLECTION_KIND_GENE_SET",
    "ENRICHMENT_COLLECTION_KIND_PTM_SET",
    "ENRICHMENT_IDENTIFIER_KIND_DISPLAY_ID",
    "ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL",
    "ENRICHMENT_IDENTIFIER_KIND_PHOSPHOSITE",
    "ENRICHMENT_IDENTIFIER_KIND_PROTEIN_ID",
    "ENRICHMENT_IDENTIFIER_KIND_SITE_KEY",
    "ENRICHMENT_METHOD_OVER_REPRESENTATION",
    "GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS",
    "MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG",
    "MULTIPLE_TESTING_CORRECTION_NONE",
    "PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS",
    "SUPPORTED_ENRICHMENT_COLLECTION_KINDS",
    "SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS",
    "SUPPORTED_ENRICHMENT_METHODS",
    "SUPPORTED_MULTIPLE_TESTING_CORRECTIONS",
    "EnrichmentCollectionKind",
    "EnrichmentIdentifierKind",
    "EnrichmentMethod",
    "EnrichmentResultRecord",
    "EnrichmentSet",
    "EnrichmentSetCollection",
    "GeneSetCollection",
    "MultipleTestingCorrection",
    "PtmSetCollection",
]
