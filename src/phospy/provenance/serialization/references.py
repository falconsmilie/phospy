"""Reference provenance payload serialization."""

from __future__ import annotations

from collections.abc import Mapping

from phospy.provenance.models import ReferenceProvenance
from phospy.provenance.reference_context import ReferenceContext
from phospy.provenance.reference_identifiers import (
    ReferenceIdentifierNormalisationRecord,
    ReferenceIdentifierNormalisationReport,
)
from phospy.provenance.serialization._payload import (
    optional_mapping,
    optional_raw_str,
    optional_str,
    require_int,
    require_mapping,
    require_raw_str,
    require_sequence,
    require_str,
    to_json_safe,
    to_json_value,
)
from phospy.provenance.serialization.tables import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)


def reference_to_payload(reference: ReferenceProvenance) -> dict[str, object]:
    return {
        "source_type": reference.source_type,
        "organism": reference.organism.value,
        "bundle_id": reference.bundle_id,
        "source_name": reference.source_name,
        "source_version": reference.source_version,
        "retrieved_at": reference.retrieved_at,
        "identifier_namespace": reference.identifier_namespace,
        "sequence_window": (
            None
            if reference.sequence_window is None
            else to_json_safe(reference.sequence_window)
        ),
        "manifest": (
            None if reference.manifest is None else to_json_safe(reference.manifest)
        ),
        "reference_context": (
            None
            if reference.reference_context is None
            else reference.reference_context.to_payload()
        ),
        "table_fingerprints": [
            table_fingerprint_to_payload(item) for item in reference.table_fingerprints
        ],
        "identifier_normalisation": (
            None
            if reference.identifier_normalisation is None
            else _reference_identifier_normalisation_to_payload(
                reference.identifier_normalisation
            )
        ),
    }


def reference_from_payload(payload: Mapping[str, object]) -> ReferenceProvenance:
    table_payload = require_sequence(
        payload.get("table_fingerprints"),
        field_name="reference_provenance.table_fingerprints",
    )
    sequence_window_payload = optional_mapping(
        payload.get("sequence_window"),
        field_name="reference_provenance.sequence_window",
    )
    manifest_payload = optional_mapping(
        payload.get("manifest"),
        field_name="reference_provenance.manifest",
    )
    reference_context_payload = optional_mapping(
        payload.get("reference_context"),
        field_name="reference_provenance.reference_context",
    )
    return ReferenceProvenance(
        source_type=require_str(
            payload.get("source_type"),
            field_name="reference_provenance.source_type",
        ),
        organism=require_str(
            payload.get("organism"),
            field_name="reference_provenance.organism",
        ),
        bundle_id=optional_str(
            payload.get("bundle_id"),
            field_name="reference_provenance.bundle_id",
        ),
        source_name=optional_str(
            payload.get("source_name"),
            field_name="reference_provenance.source_name",
        ),
        source_version=optional_str(
            payload.get("source_version"),
            field_name="reference_provenance.source_version",
        ),
        retrieved_at=optional_str(
            payload.get("retrieved_at"),
            field_name="reference_provenance.retrieved_at",
        ),
        identifier_namespace=optional_str(
            payload.get("identifier_namespace"),
            field_name="reference_provenance.identifier_namespace",
        ),
        sequence_window=None
        if sequence_window_payload is None
        else {
            key: to_json_value(value) for key, value in sequence_window_payload.items()
        },
        manifest=None
        if manifest_payload is None
        else {key: to_json_value(value) for key, value in manifest_payload.items()},
        reference_context=None
        if reference_context_payload is None
        else ReferenceContext.from_payload(reference_context_payload),
        table_fingerprints=tuple(
            table_fingerprint_from_payload(
                require_mapping(
                    item,
                    field_name=f"reference_provenance.table_fingerprints[{position}]",
                )
            )
            for position, item in enumerate(table_payload)
        ),
        identifier_normalisation=_optional_reference_identifier_normalisation_from_payload(
            payload.get("identifier_normalisation"),
            field_name="reference_provenance.identifier_normalisation",
        ),
    )


def _reference_identifier_normalisation_to_payload(
    report: ReferenceIdentifierNormalisationReport,
) -> dict[str, object]:
    return {
        "schema_version": int(report.schema_version),
        "original_row_count": int(report.original_row_count),
        "normalised_row_count": int(report.normalised_row_count),
        "invalid_identifier_count": int(report.invalid_identifier_count),
        "changed_identifier_count": int(report.changed_identifier_count),
        "duplicate_identifier_count": int(report.duplicate_identifier_count),
        "conflict_count": int(report.conflict_count),
        "records": [
            _reference_identifier_normalisation_record_to_payload(record)
            for record in report.records
        ],
    }


def _reference_identifier_normalisation_record_to_payload(
    record: ReferenceIdentifierNormalisationRecord,
) -> dict[str, object]:
    return {
        "table_name": record.table_name,
        "column_name": record.column_name,
        "row_position": int(record.row_position),
        "identifier_kind": record.identifier_kind,
        "original_value": record.original_value,
        "normalised_value": record.normalised_value,
        "status": record.status,
        "reason": record.reason,
    }


def _optional_reference_identifier_normalisation_from_payload(
    value: object,
    *,
    field_name: str,
) -> ReferenceIdentifierNormalisationReport | None:
    if value is None:
        return None
    return _reference_identifier_normalisation_from_payload(
        require_mapping(value, field_name=field_name)
    )


def _reference_identifier_normalisation_from_payload(
    payload: Mapping[str, object],
) -> ReferenceIdentifierNormalisationReport:
    records_payload = require_sequence(
        payload.get("records"),
        field_name="reference_identifier_normalisation.records",
    )
    return ReferenceIdentifierNormalisationReport(
        schema_version=require_int(
            payload.get("schema_version"),
            field_name="reference_identifier_normalisation.schema_version",
        ),
        original_row_count=require_int(
            payload.get("original_row_count"),
            field_name="reference_identifier_normalisation.original_row_count",
        ),
        normalised_row_count=require_int(
            payload.get("normalised_row_count"),
            field_name="reference_identifier_normalisation.normalised_row_count",
        ),
        invalid_identifier_count=require_int(
            payload.get("invalid_identifier_count"),
            field_name="reference_identifier_normalisation.invalid_identifier_count",
        ),
        changed_identifier_count=require_int(
            payload.get("changed_identifier_count"),
            field_name="reference_identifier_normalisation.changed_identifier_count",
        ),
        duplicate_identifier_count=require_int(
            payload.get("duplicate_identifier_count"),
            field_name=(
                "reference_identifier_normalisation.duplicate_identifier_count"
            ),
        ),
        conflict_count=require_int(
            payload.get("conflict_count"),
            field_name="reference_identifier_normalisation.conflict_count",
        ),
        records=tuple(
            _reference_identifier_normalisation_record_from_payload(
                require_mapping(
                    item,
                    field_name=(
                        f"reference_identifier_normalisation.records[{position}]"
                    ),
                )
            )
            for position, item in enumerate(records_payload)
        ),
    )


def _reference_identifier_normalisation_record_from_payload(
    payload: Mapping[str, object],
) -> ReferenceIdentifierNormalisationRecord:
    return ReferenceIdentifierNormalisationRecord(
        table_name=require_str(
            payload.get("table_name"),
            field_name="reference_identifier_normalisation_record.table_name",
        ),
        column_name=require_str(
            payload.get("column_name"),
            field_name="reference_identifier_normalisation_record.column_name",
        ),
        row_position=require_int(
            payload.get("row_position"),
            field_name="reference_identifier_normalisation_record.row_position",
        ),
        identifier_kind=require_str(
            payload.get("identifier_kind"),
            field_name="reference_identifier_normalisation_record.identifier_kind",
        ),
        original_value=require_raw_str(
            payload.get("original_value"),
            field_name="reference_identifier_normalisation_record.original_value",
        ),
        normalised_value=optional_raw_str(
            payload.get("normalised_value"),
            field_name="reference_identifier_normalisation_record.normalised_value",
        ),
        status=require_str(
            payload.get("status"),
            field_name="reference_identifier_normalisation_record.status",
        ),
        reason=optional_str(
            payload.get("reason"),
            field_name="reference_identifier_normalisation_record.reason",
        ),
    )
