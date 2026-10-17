"""Public Enrichment workflow result contracts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import cast

import pandas as pd

from phospy.enrichment.contracts.collections import (
    EnrichmentSetCollection,
    _normalise_identifier_sequence,
    _normalise_optional_multiple_testing_correction,
    _normalise_optional_non_negative_float,
    _normalise_optional_unit_interval,
    _require_identifier_kind,
    _require_int_at_least,
    _require_non_empty_string,
    _require_supported_value,
)
from phospy.enrichment.contracts.config import EnrichmentConfig
from phospy.enrichment.contracts.identifier_sets import (
    EnrichmentIdentifierSetProvenance,
)
from phospy.enrichment.contracts.vocabulary import (
    ENRICHMENT_COLLECTION_KIND_GENE_SET,
    ENRICHMENT_COLLECTION_KIND_PTM_SET,
    GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS,
    SUPPORTED_ENRICHMENT_COLLECTION_KINDS,
    SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS,
    EnrichmentCollectionKind,
    EnrichmentIdentifierKind,
    MultipleTestingCorrection,
)
from phospy.errors.validation import ContractValidationError
from phospy.frames.comparison import dataframe_equals
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.provenance.immutability import freeze_json_mapping_with_error_type
from phospy.provenance.result_caveats import ResultCaveat, validate_result_caveats
from phospy.run_lineage.models import (
    RunProvenance,
)


def _empty_json_mapping() -> Mapping[str, object]:
    return {}


@dataclass(frozen=True, slots=True)
class EnrichmentResultRecord:
    """One enrichment term row."""

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
    correction_method: MultipleTestingCorrection | None = None
    enrichment_ratio: float | None = None

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
        object.__setattr__(
            self,
            "correction_method",
            _normalise_optional_multiple_testing_correction(
                self.correction_method,
                field_name="enrichment_result_record.correction_method",
            ),
        )
        object.__setattr__(
            self,
            "enrichment_ratio",
            _normalise_optional_non_negative_float(
                self.enrichment_ratio,
                field_name="enrichment_result_record.enrichment_ratio",
            ),
        )


@dataclass(frozen=True, slots=True, eq=False)
class EnrichmentWorkflowResult:
    """Top-level native enrichment result container.

    The result contract stores an explicit enrichment result shape. Direct
    construction validates only local container consistency; no enrichment
    statistics are calculated here.

    ``table``, ``result_table``, and ``to_dataframe()`` return in-memory
    defensive snapshots only. Exporting, formatting, plotting, and report
    generation belong to IO or presentation adapters.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit content comparison.
    """

    __hash__ = object.__hash__

    identifier_kind: EnrichmentIdentifierKind
    set_collection: EnrichmentSetCollection
    config: EnrichmentConfig
    records: tuple[EnrichmentResultRecord, ...] = ()
    unmatched_identifiers: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    caveats: tuple[ResultCaveat, ...] = ()
    diagnostics: Mapping[str, object] = field(default_factory=_empty_json_mapping)
    method_metadata: Mapping[str, object] = field(default_factory=_empty_json_mapping)
    background_summary: Mapping[str, object] = field(
        default_factory=_empty_json_mapping
    )
    set_collection_summary: Mapping[str, object] = field(
        default_factory=_empty_json_mapping
    )
    selected_identifier_provenance: EnrichmentIdentifierSetProvenance | None = None
    background_identifier_provenance: EnrichmentIdentifierSetProvenance | None = None
    provenance: RunProvenance | None = None
    _result_table: pd.DataFrame = field(init=False, repr=False)

    def __post_init__(self) -> None:
        identifier_kind = _require_enrichment_identifier_kind(
            self.identifier_kind,
            field_name="enrichment_result.identifier_kind",
        )
        set_collection = _require_enrichment_collection_matches_identifier_kind(
            self.set_collection,
            identifier_kind=identifier_kind,
            field_name="enrichment_result.set_collection",
        )
        config = _require_enrichment_config(
            self.config,
            field_name="enrichment_result.config",
        )
        records = tuple(
            _require_enrichment_result_record(
                record,
                field_name="enrichment_result.records",
            )
            for record in self.records
        )
        for record in records:
            if record.identifier_kind != identifier_kind:
                raise ContractValidationError(
                    "enrichment_result.records identifier_kind values must match "
                    "enrichment_result.identifier_kind"
                )
            if record.collection_kind != set_collection.collection_kind:
                raise ContractValidationError(
                    "enrichment_result.records collection_kind values must match "
                    "enrichment_result.set_collection"
                )
        unmatched_identifiers = _normalise_enrichment_identifier_sequence(
            self.unmatched_identifiers,
            field_name="enrichment_result.unmatched_identifiers",
            allow_empty=True,
        )
        warnings = tuple(_validate_enrichment_warning(value) for value in self.warnings)
        caveats = validate_result_caveats(
            self.caveats,
            field_name="enrichment_result.caveats",
            error_type=ContractValidationError,
        )
        diagnostics = _freeze_enrichment_json_mapping(
            self.diagnostics,
            field_name="enrichment_result.diagnostics",
        )
        method_metadata = _freeze_enrichment_json_mapping(
            self.method_metadata,
            field_name="enrichment_result.method_metadata",
        )
        background_summary = _freeze_enrichment_json_mapping(
            self.background_summary,
            field_name="enrichment_result.background_summary",
        )
        set_collection_summary = _freeze_enrichment_json_mapping(
            self.set_collection_summary,
            field_name="enrichment_result.set_collection_summary",
        )
        selected_identifier_provenance = _validate_optional_identifier_set_provenance(
            self.selected_identifier_provenance,
            field_name="enrichment_result.selected_identifier_provenance",
        )
        background_identifier_provenance = _validate_optional_identifier_set_provenance(
            self.background_identifier_provenance,
            field_name="enrichment_result.background_identifier_provenance",
        )
        provenance = _validate_optional_run_provenance(
            self.provenance,
            field_name="enrichment_result.provenance",
        )
        result_table = own_dataframe(
            _enrichment_records_to_dataframe(records),
            field_name="enrichment_result.table",
            error_type=ContractValidationError,
            assume_owned=True,
        )
        object.__setattr__(self, "identifier_kind", identifier_kind)
        object.__setattr__(self, "set_collection", set_collection)
        object.__setattr__(self, "config", config)
        object.__setattr__(self, "records", records)
        object.__setattr__(self, "unmatched_identifiers", unmatched_identifiers)
        object.__setattr__(self, "warnings", warnings)
        object.__setattr__(self, "caveats", caveats)
        object.__setattr__(self, "diagnostics", diagnostics)
        object.__setattr__(self, "method_metadata", method_metadata)
        object.__setattr__(self, "background_summary", background_summary)
        object.__setattr__(
            self,
            "set_collection_summary",
            set_collection_summary,
        )
        object.__setattr__(
            self,
            "selected_identifier_provenance",
            selected_identifier_provenance,
        )
        object.__setattr__(
            self,
            "background_identifier_provenance",
            background_identifier_provenance,
        )
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "_result_table", result_table)

    @property
    def table(self) -> pd.DataFrame:
        """Return a defensive snapshot of the enrichment result table."""

        return export_dataframe(self._result_table)

    @property
    def result_table(self) -> pd.DataFrame:
        """Return a defensive snapshot of the enrichment result table."""

        return export_dataframe(self._result_table)

    def to_dataframe(self) -> pd.DataFrame:
        """Return a defensive snapshot of the enrichment result table."""

        return export_dataframe(self._result_table)

    def scientifically_equals(
        self,
        other: object,
        *,
        include_provenance: bool = True,
    ) -> bool:
        """Return ``True`` when another enrichment result has the same content."""

        if not isinstance(other, EnrichmentWorkflowResult):
            return False
        same_content = (
            self.identifier_kind == other.identifier_kind
            and self.set_collection == other.set_collection
            and self.config == other.config
            and self.records == other.records
            and self.unmatched_identifiers == other.unmatched_identifiers
            and self.warnings == other.warnings
            and self.caveats == other.caveats
            and self.diagnostics == other.diagnostics
            and self.method_metadata == other.method_metadata
            and self.background_summary == other.background_summary
            and self.set_collection_summary == other.set_collection_summary
            and self.selected_identifier_provenance
            == other.selected_identifier_provenance
            and self.background_identifier_provenance
            == other.background_identifier_provenance
            and dataframe_equals(self._result_table, other._result_table)
        )
        if not same_content:
            return False
        if include_provenance and self.provenance != other.provenance:
            return False
        return True


def _enrichment_records_to_dataframe(
    records: tuple[EnrichmentResultRecord, ...],
) -> pd.DataFrame:
    columns = [
        "term_id",
        "term_name",
        "collection_kind",
        "identifier_kind",
        "input_overlap_count",
        "background_overlap_count",
        "set_size",
        "overlap_identifiers",
        "p_value",
        "adjusted_p_value",
        "correction_method",
        "enrichment_ratio",
    ]
    rows = [
        {
            "term_id": record.term_id,
            "term_name": record.term_name,
            "collection_kind": record.collection_kind,
            "identifier_kind": record.identifier_kind,
            "input_overlap_count": record.input_overlap_count,
            "background_overlap_count": record.background_overlap_count,
            "set_size": record.set_size,
            "overlap_identifiers": record.overlap_identifiers,
            "p_value": record.p_value,
            "adjusted_p_value": record.adjusted_p_value,
            "correction_method": record.correction_method,
            "enrichment_ratio": record.enrichment_ratio,
        }
        for record in records
    ]
    return pd.DataFrame(rows, columns=columns)


def _validate_enrichment_warning(value: object) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise ContractValidationError(
            "enrichment_result.warnings must contain non-empty strings"
        )
    return value.strip()


def _require_enrichment_identifier_kind(
    value: object,
    *,
    field_name: str,
) -> EnrichmentIdentifierKind:
    if value in SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS:
        return value
    allowed = ", ".join(sorted(SUPPORTED_ENRICHMENT_IDENTIFIER_KINDS))
    raise ContractValidationError(f"{field_name} must be one of: {allowed}")


def _normalise_enrichment_identifier_sequence(
    value: object,
    *,
    field_name: str,
    allow_empty: bool = False,
) -> tuple[str, ...]:
    if isinstance(value, str | bytes | bytearray) or not isinstance(value, Sequence):
        raise ContractValidationError(f"{field_name} must be a sequence of strings")
    sequence = cast(Sequence[object], value)
    normalised = tuple(
        _require_non_empty_text(item, field_name=f"{field_name}[]") for item in sequence
    )
    if not allow_empty and not normalised:
        raise ContractValidationError(f"{field_name} must not be empty")
    return normalised


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str) or value.strip() == "":
        raise ContractValidationError(f"{field_name} must be a non-empty string")
    return value


def _require_enrichment_collection_matches_identifier_kind(
    collection: object,
    *,
    identifier_kind: EnrichmentIdentifierKind,
    field_name: str,
) -> EnrichmentSetCollection:
    if not isinstance(collection, EnrichmentSetCollection):
        raise ContractValidationError(
            f"{field_name} must be EnrichmentSetCollection, "
            "GeneSetCollection, or PtmSetCollection"
        )
    if collection.identifier_kind != identifier_kind:
        raise ContractValidationError(
            f"{field_name}.identifier_kind must match result identifier_kind; "
            f"observed={collection.identifier_kind!r}, expected={identifier_kind!r}"
        )
    expected_collection_kind = _collection_kind_for_identifier_kind(identifier_kind)
    if collection.collection_kind != expected_collection_kind:
        raise ContractValidationError(
            f"{field_name}.collection_kind must match identifier_kind; "
            f"observed={collection.collection_kind!r}, "
            f"expected={expected_collection_kind!r}"
        )
    return collection


def _collection_kind_for_identifier_kind(
    identifier_kind: EnrichmentIdentifierKind,
) -> str:
    if identifier_kind in GENE_LEVEL_ENRICHMENT_IDENTIFIER_KINDS:
        return ENRICHMENT_COLLECTION_KIND_GENE_SET
    return ENRICHMENT_COLLECTION_KIND_PTM_SET


def _require_enrichment_config(
    value: object,
    *,
    field_name: str,
) -> EnrichmentConfig:
    if not isinstance(value, EnrichmentConfig):
        raise ContractValidationError(f"{field_name} must be EnrichmentConfig")
    return value


def _require_enrichment_result_record(
    value: object,
    *,
    field_name: str,
) -> EnrichmentResultRecord:
    if not isinstance(value, EnrichmentResultRecord):
        raise ContractValidationError(
            f"{field_name} must contain EnrichmentResultRecord values"
        )
    return value


def _freeze_enrichment_json_mapping(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object]:
    return freeze_json_mapping_with_error_type(
        value,
        field_name=field_name,
        error_type=ContractValidationError,
    )


def _validate_optional_identifier_set_provenance(
    value: object | None,
    *,
    field_name: str,
) -> EnrichmentIdentifierSetProvenance | None:
    if value is None or isinstance(value, EnrichmentIdentifierSetProvenance):
        return value
    raise ContractValidationError(
        f"{field_name} must be EnrichmentIdentifierSetProvenance or None"
    )


def _validate_optional_run_provenance(
    value: object | None,
    *,
    field_name: str,
) -> RunProvenance | None:
    if value is None or isinstance(value, RunProvenance):
        return value
    raise ContractValidationError(f"{field_name} must be RunProvenance or None")


__all__ = [
    "EnrichmentResultRecord",
    "EnrichmentWorkflowResult",
]
