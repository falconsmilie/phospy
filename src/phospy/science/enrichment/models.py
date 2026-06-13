"""Domain contracts for enrichment workflows."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
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
class GeneSetCollection:
    """Explicit gene-level enrichment collection.

    The identifiers inside ``sets`` are gene/protein-level identifiers. This
    object is intentionally separate from ``PtmSetCollection`` so a request
    cannot silently mix gene-set and phosphosite/PTM-set semantics.
    """

    sets: Mapping[str, Sequence[str]]
    identifier_kind: EnrichmentIdentifierKind = ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    term_names: Mapping[str, str] = field(default_factory=dict)
    source_name: str = "user"
    source_version: str | None = None

    def __post_init__(self) -> None:
        identifier_kind = _require_identifier_kind(
            self.identifier_kind,
            field_name="gene_set_collection.identifier_kind",
        )
        if identifier_kind not in GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
            allowed = ", ".join(sorted(GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS))
            raise WorkflowValidationError(
                f"gene_set_collection.identifier_kind must be gene-level: {allowed}"
            )
        sets = _normalise_set_mapping(
            self.sets,
            field_name="gene_set_collection.sets",
        )
        term_names = _normalise_term_names(
            self.term_names,
            allowed_term_ids=frozenset(sets),
            field_name="gene_set_collection.term_names",
        )
        object.__setattr__(self, "sets", sets)
        object.__setattr__(self, "identifier_kind", identifier_kind)
        object.__setattr__(self, "term_names", term_names)
        object.__setattr__(
            self,
            "source_name",
            _require_non_empty_string(
                self.source_name,
                field_name="gene_set_collection.source_name",
            ),
        )
        object.__setattr__(
            self,
            "source_version",
            _normalise_optional_string(
                self.source_version,
                field_name="gene_set_collection.source_version",
            ),
        )

    @property
    def collection_kind(self) -> EnrichmentCollectionKind:
        return ENRICHMENT_COLLECTION_KIND_GENE_SET


@dataclass(frozen=True, slots=True)
class PtmSetCollection:
    """Explicit phosphosite/PTM-level enrichment collection."""

    sets: Mapping[str, Sequence[str]]
    identifier_kind: EnrichmentIdentifierKind = ENRICHMENT_IDENTIFIER_KIND_SITE_KEY
    term_names: Mapping[str, str] = field(default_factory=dict)
    source_name: str = "user"
    source_version: str | None = None

    def __post_init__(self) -> None:
        identifier_kind = _require_identifier_kind(
            self.identifier_kind,
            field_name="ptm_set_collection.identifier_kind",
        )
        if identifier_kind not in PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
            allowed = ", ".join(sorted(PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS))
            raise WorkflowValidationError(
                f"ptm_set_collection.identifier_kind must be PTM-level: {allowed}"
            )
        sets = _normalise_set_mapping(
            self.sets,
            field_name="ptm_set_collection.sets",
        )
        term_names = _normalise_term_names(
            self.term_names,
            allowed_term_ids=frozenset(sets),
            field_name="ptm_set_collection.term_names",
        )
        object.__setattr__(self, "sets", sets)
        object.__setattr__(self, "identifier_kind", identifier_kind)
        object.__setattr__(self, "term_names", term_names)
        object.__setattr__(
            self,
            "source_name",
            _require_non_empty_string(
                self.source_name,
                field_name="ptm_set_collection.source_name",
            ),
        )
        object.__setattr__(
            self,
            "source_version",
            _normalise_optional_string(
                self.source_version,
                field_name="ptm_set_collection.source_version",
            ),
        )

    @property
    def collection_kind(self) -> EnrichmentCollectionKind:
        return ENRICHMENT_COLLECTION_KIND_PTM_SET


EnrichmentSetCollection = GeneSetCollection | PtmSetCollection


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


def _require_collection_matches_identifier_kind(
    collection: object,
    *,
    identifier_kind: EnrichmentIdentifierKind,
    field_name: str,
) -> EnrichmentSetCollection:
    if not isinstance(collection, GeneSetCollection | PtmSetCollection):
        raise WorkflowValidationError(
            f"{field_name} must be GeneSetCollection or PtmSetCollection"
        )
    if collection.identifier_kind != identifier_kind:
        raise WorkflowValidationError(
            f"{field_name}.identifier_kind must match request identifier_kind; "
            f"observed={collection.identifier_kind!r}, expected={identifier_kind!r}"
        )
    if (
        isinstance(collection, GeneSetCollection)
        and identifier_kind not in GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS
    ):
        raise WorkflowValidationError(
            "gene-level enrichment requires identifier_kind 'gene_symbol' or "
            "'protein_id'"
        )
    if (
        isinstance(collection, PtmSetCollection)
        and identifier_kind not in PTM_LEVEL_ENRICHMENT_IDENTIFIER_KINDS
    ):
        raise WorkflowValidationError(
            "PTM-set enrichment requires identifier_kind 'site_key', 'display_id', "
            "or 'phosphosite'"
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
    "EnrichmentSetCollection",
    "GeneSetCollection",
    "MultipleTestingCorrection",
    "PtmSetCollection",
]
