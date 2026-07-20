"""Public enrichment workflow result contracts."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import pandas as pd

from phospy.contracts.configs import EnrichmentConfig
from phospy.contracts.enrichment_identifier_sets import (
    EnrichmentIdentifierSetProvenance,
)
from phospy.contracts.result_caveats import ResultCaveat, validate_result_caveats
from phospy.errors.validation import ContractValidationError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.provenance.immutability import freeze_json_mapping_with_error_type
from phospy.provenance.models import RunProvenance
from phospy.science.enrichment.models import (
    EnrichmentIdentifierKind,
    EnrichmentResultRecord,
    EnrichmentSetCollection,
    _normalise_identifier_sequence,
    _require_collection_matches_identifier_kind,
    _require_identifier_kind,
)


@dataclass(frozen=True, slots=True)
class EnrichmentWorkflowResult:
    """Top-level native enrichment result container.

    The result contract stores an explicit enrichment result shape. Direct
    construction validates only local container consistency; no enrichment
    statistics are calculated here.

    ``table``, ``result_table``, and ``to_dataframe()`` return in-memory
    defensive snapshots only. Exporting, formatting, plotting, and report
    generation belong to IO or presentation adapters.
    """

    identifier_kind: EnrichmentIdentifierKind
    set_collection: EnrichmentSetCollection
    config: EnrichmentConfig
    records: tuple[EnrichmentResultRecord, ...] = ()
    unmatched_identifiers: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    caveats: tuple[ResultCaveat, ...] = ()
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    method_metadata: Mapping[str, object] = field(default_factory=dict)
    background_summary: Mapping[str, object] = field(default_factory=dict)
    set_collection_summary: Mapping[str, object] = field(default_factory=dict)
    selected_identifier_provenance: EnrichmentIdentifierSetProvenance | None = None
    background_identifier_provenance: EnrichmentIdentifierSetProvenance | None = None
    provenance: RunProvenance | None = None
    _result_table: pd.DataFrame = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        identifier_kind = _require_identifier_kind(
            self.identifier_kind,
            field_name="enrichment_result.identifier_kind",
        )
        set_collection = _require_collection_matches_identifier_kind(
            self.set_collection,
            identifier_kind=identifier_kind,
            field_name="enrichment_result.set_collection",
        )
        if not isinstance(self.config, EnrichmentConfig):
            raise ContractValidationError(
                "enrichment_result.config must be EnrichmentConfig"
            )
        records = tuple(self.records)
        for record in records:
            if not isinstance(record, EnrichmentResultRecord):
                raise ContractValidationError(
                    "enrichment_result.records must contain "
                    "EnrichmentResultRecord values"
                )
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
        unmatched_identifiers = _normalise_identifier_sequence(
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
        if self.provenance is not None and not isinstance(
            self.provenance, RunProvenance
        ):
            raise ContractValidationError(
                "enrichment_result.provenance must be RunProvenance or None"
            )
        result_table = own_dataframe(
            _enrichment_records_to_dataframe(records),
            field_name="enrichment_result.table",
            error_type=ContractValidationError,
            assume_owned=True,
        )
        object.__setattr__(self, "identifier_kind", identifier_kind)
        object.__setattr__(self, "set_collection", set_collection)
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


__all__ = [
    "EnrichmentResultRecord",
    "EnrichmentWorkflowResult",
]
