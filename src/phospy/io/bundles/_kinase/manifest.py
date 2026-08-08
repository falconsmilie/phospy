"""Current-contract manifest serialization and parsing for kinase bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import NoReturn

from phospy.contracts.results import KinaseWorkflowResult
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._kinase.constants import (
    KINASE_BUNDLE_KIND,
    KINASE_BUNDLE_MANIFEST_VERSION,
)
from phospy.io.bundles._shared.integrity import (
    require_file_entry,
    require_optional_table_entry,
    require_table_entry,
)
from phospy.io.bundles._shared.intensity_scale_state import (
    intensity_scale_state_to_payload,
)
from phospy.io.bundles._shared.primitives import (
    require_bool,
    require_int,
    require_mapping,
    require_str,
)
from phospy.io.bundles._shared.processing_state import (
    processing_state_to_payload,
)
from phospy.provenance.serialization import to_payload as provenance_to_payload


@dataclass(frozen=True, slots=True)
class KinaseManifestSections:
    """Decoded current manifest sections needed to load kinase bundles."""

    manifest_version: int
    dataset_metadata: Mapping[str, object]
    dataset_tables: Mapping[str, object]
    references_metadata: Mapping[str, object]
    reference_tables: Mapping[str, object]
    scoring_tables: Mapping[str, object]
    prediction_tables: Mapping[str, object]
    activity_enabled: bool
    activity_method_metadata: Mapping[str, object] | None
    activity_method_summary: Mapping[str, object] | None
    activity_input_semantics: Mapping[str, object] | None
    activity_profile_metadata: Mapping[str, object] | None
    activity_membership_selection: Mapping[str, object] | None
    activity_tables: Mapping[str, object]
    provenance_payload: Mapping[str, object]
    config_snapshot_entry: Mapping[str, object]
    caveats_payload: object = ()


@dataclass(frozen=True, slots=True)
class _KinaseManifestSchema:
    """Schema-parsed manifest payload before file-record validation."""

    manifest_version: int
    dataset_metadata: Mapping[str, object]
    references_metadata: Mapping[str, object]
    dataset_tables: Mapping[str, object]
    reference_tables: Mapping[str, object]
    scoring_tables: Mapping[str, object]
    prediction_tables: Mapping[str, object]
    activity_enabled: bool
    activity_method_metadata: Mapping[str, object] | None
    activity_method_summary: Mapping[str, object] | None
    activity_input_semantics: Mapping[str, object] | None
    activity_profile_metadata: Mapping[str, object] | None
    activity_membership_selection: Mapping[str, object] | None
    activity_tables: Mapping[str, object]
    provenance_payload: Mapping[str, object]
    config_snapshot_raw: object
    caveats_payload: object


@dataclass(frozen=True, slots=True)
class _KinaseManifestFileRecords:
    """Manifest table/file records after path and record validation."""

    dataset_tables: Mapping[str, object]
    reference_tables: Mapping[str, object]
    scoring_tables: Mapping[str, object]
    prediction_tables: Mapping[str, object]
    activity_tables: Mapping[str, object]
    config_snapshot_entry: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class _KinaseManifestTopLevel:
    """Top-level manifest fields after schema and version validation."""

    manifest_version: int
    dataset_payload: Mapping[str, object]
    references_payload: Mapping[str, object]
    outputs_payload: Mapping[str, object]
    provenance_payload: Mapping[str, object]
    config_snapshot_raw: object
    caveats_payload: object


@dataclass(frozen=True, slots=True)
class _KinaseMetadataSections:
    """Scientific metadata sections referenced by the manifest."""

    dataset_metadata: Mapping[str, object]
    references_metadata: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class _KinaseOutputPayloads:
    """Output subsection payloads after output schema validation."""

    scoring_payload: Mapping[str, object]
    prediction_payload: Mapping[str, object]
    activity_payload: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class _KinaseActivityMetadata:
    """Activity semantic metadata after enabled/disabled consistency validation."""

    enabled: bool
    method_metadata: Mapping[str, object] | None
    method_summary: Mapping[str, object] | None
    input_semantics: Mapping[str, object] | None
    profile_metadata: Mapping[str, object] | None
    membership_selection: Mapping[str, object] | None


@dataclass(frozen=True, slots=True)
class _KinaseTableSections:
    """Manifest table sections after table-key schema validation."""

    dataset_tables: Mapping[str, object]
    reference_tables: Mapping[str, object]
    scoring_tables: Mapping[str, object]
    prediction_tables: Mapping[str, object]
    activity_tables: Mapping[str, object]


_LEGACY_KINASE_BUNDLE_SCHEMA_ERROR = (
    "Legacy kinase bundle schemas are no longer supported. Regenerate the bundle "
    "with the current PhosPy version."
)
_MANIFEST_ALLOWED_FIELDS = frozenset(
    {
        "bundle_type",
        "manifest_version",
        "table_format",
        "dataset",
        "resolved_references",
        "outputs",
        "provenance",
        "config_snapshot",
        "caveats",
    }
)
_MANIFEST_REQUIRED_FIELDS = frozenset(
    {
        "bundle_type",
        "manifest_version",
        "table_format",
        "dataset",
        "resolved_references",
        "outputs",
        "provenance",
        "config_snapshot",
    }
)
_DATASET_ALLOWED_FIELDS = frozenset({"metadata", "tables"})
_DATASET_TABLE_KEYS = frozenset(
    {"phospho", "site_metadata", "sample_metadata", "total"}
)
_REFERENCES_ALLOWED_FIELDS = frozenset({"metadata", "tables"})
_REFERENCE_TABLE_KEYS = frozenset({"kinase_substrate_map", "site_sequences"})
_OUTPUTS_ALLOWED_FIELDS = frozenset({"scoring", "prediction", "activity"})
_SCORING_ALLOWED_FIELDS = frozenset({"tables"})
_SCORING_TABLE_KEYS = frozenset(
    {
        "profile_scores",
        "motif_scores",
        "rank_weighted_fusion_scores",
        "kinase_library_motif_scores",
        "combined_profile_motif_scores",
        "score_fusion_weights",
        "kinase_library_site_diagnostics",
        "kinase_library_kinase_diagnostics",
        "substrate_contributions",
    }
)
_SCORING_TABLE_REQUIRED_KEYS = frozenset(
    {
        "profile_scores",
        "motif_scores",
        "rank_weighted_fusion_scores",
        "score_fusion_weights",
    }
)
_PREDICTION_ALLOWED_FIELDS = frozenset({"tables"})
_PREDICTION_TABLE_KEYS = frozenset({"pred_mat", "substrate_list"})
_ACTIVITY_ALLOWED_FIELDS = frozenset(
    {
        "enabled",
        "method",
        "summary",
        "input_semantics",
        "profile_metadata",
        "membership_selection",
        "tables",
    }
)
_ACTIVITY_REQUIRED_FIELDS = _ACTIVITY_ALLOWED_FIELDS - frozenset(
    {"membership_selection"}
)
_ACTIVITY_TABLE_KEYS = frozenset(
    {
        "weighted_activity",
        "thresholded_substrate_mean_activity",
        "thresholded_substrate_counts",
        "activity_substrate_counts",
        "target_counts",
        "target_table",
        "statistics_table",
    }
)


def build_manifest(
    *,
    result: KinaseWorkflowResult,
    table_format: str,
    dataset_tables: Mapping[str, object],
    reference_tables: Mapping[str, object],
    scoring_tables: Mapping[str, object],
    prediction_tables: Mapping[str, object],
    activity_tables: Mapping[str, object],
    config_snapshot_entry: Mapping[str, object],
) -> dict[str, object]:
    """Build the current manifest payload from bundle contract data."""

    if result.provenance is None:
        raise PhosPyInputError(
            "kinase bundle saving requires result.provenance; "
            "bundle manifests must include explicit provenance metadata"
        )

    return {
        "bundle_type": KINASE_BUNDLE_KIND,
        "manifest_version": KINASE_BUNDLE_MANIFEST_VERSION,
        "table_format": table_format,
        "dataset": {
            "metadata": {
                "organism": (
                    None
                    if result.dataset.organism is None
                    else result.dataset.organism.value
                ),
                "intensity_scale_state": intensity_scale_state_to_payload(
                    result.dataset.intensity_scale_state
                ),
                "processing_state": processing_state_to_payload(
                    result.dataset.processing_state
                ),
            },
            "tables": dict(dataset_tables),
        },
        "resolved_references": {
            "metadata": {
                "organism": result.references.organism.value,
            },
            "tables": dict(reference_tables),
        },
        "outputs": {
            "scoring": {
                "tables": dict(scoring_tables),
            },
            "prediction": {
                "tables": dict(prediction_tables),
            },
            "activity": {
                "enabled": result.activity_result is not None,
                "method": (
                    None
                    if result.activity_result is None
                    else result.activity_result.activity_method.to_payload()
                ),
                "summary": (
                    None
                    if result.activity_result is None
                    or result.activity_result.method_summary is None
                    else result.activity_result.method_summary.to_payload()
                ),
                "input_semantics": (
                    None
                    if result.activity_result is None
                    else result.activity_result.input_semantics.to_payload()
                ),
                "profile_metadata": (
                    None
                    if result.activity_result is None
                    else result.activity_result.profile_metadata.to_payload()
                ),
                "membership_selection": (
                    None
                    if result.activity_result is None
                    or result.activity_result.membership_selection is None
                    else result.activity_result.membership_selection.to_payload()
                ),
                "tables": dict(activity_tables),
            },
        },
        "provenance": provenance_to_payload(result.provenance),
        "caveats": [caveat.to_payload() for caveat in result.caveats],
        "config_snapshot": dict(config_snapshot_entry),
    }


def parse_manifest(payload: Mapping[str, object]) -> KinaseManifestSections:
    """Parse and validate current-contract kinase manifest payload."""

    schema = _parse_manifest_schema(payload)
    _validate_manifest_paths(schema)
    file_records = _validate_manifest_file_records(schema)
    return _assemble_manifest_sections(schema=schema, file_records=file_records)


def _parse_manifest_schema(payload: Mapping[str, object]) -> _KinaseManifestSchema:
    """Parse manifest object schemas through typed section stages."""

    top_level = _parse_top_level_schema(payload)
    metadata = _parse_scientific_metadata_sections(top_level)
    outputs = _parse_output_payloads(top_level.outputs_payload)
    activity = _parse_activity_metadata(outputs.activity_payload)
    tables = _parse_table_sections(
        top_level=top_level,
        outputs=outputs,
    )
    return _KinaseManifestSchema(
        manifest_version=top_level.manifest_version,
        dataset_metadata=metadata.dataset_metadata,
        references_metadata=metadata.references_metadata,
        dataset_tables=tables.dataset_tables,
        reference_tables=tables.reference_tables,
        scoring_tables=tables.scoring_tables,
        prediction_tables=tables.prediction_tables,
        activity_enabled=activity.enabled,
        activity_method_metadata=activity.method_metadata,
        activity_method_summary=activity.method_summary,
        activity_input_semantics=activity.input_semantics,
        activity_profile_metadata=activity.profile_metadata,
        activity_membership_selection=activity.membership_selection,
        activity_tables=tables.activity_tables,
        provenance_payload=top_level.provenance_payload,
        config_snapshot_raw=top_level.config_snapshot_raw,
        caveats_payload=top_level.caveats_payload,
    )


def _parse_top_level_schema(
    payload: Mapping[str, object],
) -> _KinaseManifestTopLevel:
    """Parse top-level manifest schema, kind, version, and JSON references."""

    _reject_unsupported_fields(
        payload,
        field_name="bundle manifest",
        allowed_fields=_MANIFEST_ALLOWED_FIELDS,
    )
    _require_fields(
        payload,
        field_name="bundle manifest",
        required_fields=_MANIFEST_REQUIRED_FIELDS,
        unsupported_shape=True,
    )
    bundle_type = require_str(
        payload.get("bundle_type"),
        field_name="bundle manifest.bundle_type",
    )
    if bundle_type != KINASE_BUNDLE_KIND:
        _raise_unsupported_manifest_shape(
            "unsupported bundle manifest bundle_type "
            f"'{bundle_type}'; expected '{KINASE_BUNDLE_KIND}'"
        )
    manifest_version = require_int(
        payload.get("manifest_version"),
        field_name="bundle manifest.manifest_version",
    )
    if manifest_version != KINASE_BUNDLE_MANIFEST_VERSION:
        _reject_unsupported_manifest_version(manifest_version)
    require_str(
        payload.get("table_format"),
        field_name="bundle manifest.table_format",
    )
    provenance_raw = payload.get("provenance")
    if provenance_raw is None:
        _raise_unsupported_manifest_shape("bundle manifest.provenance is required")
    return _KinaseManifestTopLevel(
        manifest_version=manifest_version,
        dataset_payload=_parse_dataset_payload(payload.get("dataset")),
        references_payload=_parse_references_payload(
            payload.get("resolved_references")
        ),
        outputs_payload=_parse_outputs_payload(payload.get("outputs")),
        provenance_payload=require_mapping(
            provenance_raw,
            field_name="bundle manifest.provenance",
        ),
        config_snapshot_raw=payload.get("config_snapshot"),
        caveats_payload=payload.get("caveats", []),
    )


def _reject_unsupported_manifest_version(manifest_version: int) -> NoReturn:
    if manifest_version == 2:
        _raise_unsupported_manifest_shape(
            "bundle manifest.manifest_version=2 is a legacy kinase bundle "
            "schema; activity input semantics, profile identity, and "
            "condition-summary aggregation metadata were not part of schema "
            "version 2, so the bundle must be regenerated with the current "
            "PhosPy version"
        )
    _raise_unsupported_manifest_shape(
        "unsupported bundle manifest version "
        f"'{manifest_version}'; expected {KINASE_BUNDLE_MANIFEST_VERSION}"
    )


def _parse_dataset_payload(value: object) -> Mapping[str, object]:
    dataset_payload = require_mapping(
        value,
        field_name="bundle manifest.dataset",
    )
    _reject_unsupported_fields(
        dataset_payload,
        field_name="bundle manifest.dataset",
        allowed_fields=_DATASET_ALLOWED_FIELDS,
    )
    _require_fields(
        dataset_payload,
        field_name="bundle manifest.dataset",
        required_fields=_DATASET_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    return dataset_payload


def _parse_references_payload(value: object) -> Mapping[str, object]:
    references_payload = require_mapping(
        value,
        field_name="bundle manifest.resolved_references",
    )
    _reject_unsupported_fields(
        references_payload,
        field_name="bundle manifest.resolved_references",
        allowed_fields=_REFERENCES_ALLOWED_FIELDS,
    )
    _require_fields(
        references_payload,
        field_name="bundle manifest.resolved_references",
        required_fields=_REFERENCES_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    return references_payload


def _parse_outputs_payload(value: object) -> Mapping[str, object]:
    outputs_payload = require_mapping(
        value,
        field_name="bundle manifest.outputs",
    )
    _reject_unsupported_fields(
        outputs_payload,
        field_name="bundle manifest.outputs",
        allowed_fields=_OUTPUTS_ALLOWED_FIELDS,
    )
    _require_fields(
        outputs_payload,
        field_name="bundle manifest.outputs",
        required_fields=_OUTPUTS_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    return outputs_payload


def _parse_scientific_metadata_sections(
    top_level: _KinaseManifestTopLevel,
) -> _KinaseMetadataSections:
    return _KinaseMetadataSections(
        dataset_metadata=require_mapping(
            top_level.dataset_payload.get("metadata"),
            field_name="bundle manifest.dataset.metadata",
        ),
        references_metadata=require_mapping(
            top_level.references_payload.get("metadata"),
            field_name="bundle manifest.resolved_references.metadata",
        ),
    )


def _parse_output_payloads(
    outputs_payload: Mapping[str, object],
) -> _KinaseOutputPayloads:
    scoring_payload = require_mapping(
        outputs_payload.get("scoring"),
        field_name="bundle manifest.outputs.scoring",
    )
    _reject_unsupported_fields(
        scoring_payload,
        field_name="bundle manifest.outputs.scoring",
        allowed_fields=_SCORING_ALLOWED_FIELDS,
    )
    _require_fields(
        scoring_payload,
        field_name="bundle manifest.outputs.scoring",
        required_fields=_SCORING_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    prediction_payload = require_mapping(
        outputs_payload.get("prediction"),
        field_name="bundle manifest.outputs.prediction",
    )
    _reject_unsupported_fields(
        prediction_payload,
        field_name="bundle manifest.outputs.prediction",
        allowed_fields=_PREDICTION_ALLOWED_FIELDS,
    )
    _require_fields(
        prediction_payload,
        field_name="bundle manifest.outputs.prediction",
        required_fields=_PREDICTION_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    activity_payload = require_mapping(
        outputs_payload.get("activity"),
        field_name="bundle manifest.outputs.activity",
    )
    _reject_unsupported_fields(
        activity_payload,
        field_name="bundle manifest.outputs.activity",
        allowed_fields=_ACTIVITY_ALLOWED_FIELDS,
    )
    _require_fields(
        activity_payload,
        field_name="bundle manifest.outputs.activity",
        required_fields=_ACTIVITY_REQUIRED_FIELDS,
        unsupported_shape=True,
    )
    return _KinaseOutputPayloads(
        scoring_payload=scoring_payload,
        prediction_payload=prediction_payload,
        activity_payload=activity_payload,
    )


def _parse_activity_metadata(
    activity_payload: Mapping[str, object],
) -> _KinaseActivityMetadata:
    activity_enabled = require_bool(
        activity_payload.get("enabled"),
        field_name="bundle manifest.outputs.activity.enabled",
    )
    activity_method_metadata = _optional_mapping_payload(
        activity_payload.get("method"),
        field_name="bundle manifest.outputs.activity.method",
    )
    activity_method_summary = _optional_mapping_payload(
        activity_payload.get("summary"),
        field_name="bundle manifest.outputs.activity.summary",
    )
    activity_input_semantics = _optional_mapping_payload(
        activity_payload.get("input_semantics"),
        field_name="bundle manifest.outputs.activity.input_semantics",
    )
    activity_profile_metadata = _optional_mapping_payload(
        activity_payload.get("profile_metadata"),
        field_name="bundle manifest.outputs.activity.profile_metadata",
    )
    activity_membership_selection = _optional_mapping_payload(
        activity_payload.get("membership_selection"),
        field_name="bundle manifest.outputs.activity.membership_selection",
    )
    _validate_activity_metadata_switches(
        activity_enabled=activity_enabled,
        activity_method_metadata=activity_method_metadata,
        activity_method_summary=activity_method_summary,
        activity_input_semantics=activity_input_semantics,
        activity_profile_metadata=activity_profile_metadata,
        activity_membership_selection=activity_membership_selection,
    )
    return _KinaseActivityMetadata(
        enabled=activity_enabled,
        method_metadata=activity_method_metadata,
        method_summary=activity_method_summary,
        input_semantics=activity_input_semantics,
        profile_metadata=activity_profile_metadata,
        membership_selection=activity_membership_selection,
    )


def _optional_mapping_payload(
    value: object,
    *,
    field_name: str,
) -> Mapping[str, object] | None:
    if value is None:
        return None
    return require_mapping(value, field_name=field_name)


def _validate_activity_metadata_switches(
    *,
    activity_enabled: bool,
    activity_method_metadata: Mapping[str, object] | None,
    activity_method_summary: Mapping[str, object] | None,
    activity_input_semantics: Mapping[str, object] | None,
    activity_profile_metadata: Mapping[str, object] | None,
    activity_membership_selection: Mapping[str, object] | None,
) -> None:
    if activity_enabled:
        _require_enabled_activity_metadata(
            activity_method_metadata=activity_method_metadata,
            activity_input_semantics=activity_input_semantics,
            activity_profile_metadata=activity_profile_metadata,
        )
        return
    _reject_disabled_activity_metadata(
        activity_method_metadata=activity_method_metadata,
        activity_method_summary=activity_method_summary,
        activity_input_semantics=activity_input_semantics,
        activity_profile_metadata=activity_profile_metadata,
        activity_membership_selection=activity_membership_selection,
    )


def _require_enabled_activity_metadata(
    *,
    activity_method_metadata: Mapping[str, object] | None,
    activity_input_semantics: Mapping[str, object] | None,
    activity_profile_metadata: Mapping[str, object] | None,
) -> None:
    if activity_method_metadata is None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.method is required when activity is enabled"
        )
    if activity_input_semantics is None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.input_semantics is required when "
            "activity is enabled; regenerate the bundle from the original "
            "KinaseActivityResult"
        )
    if activity_profile_metadata is None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.profile_metadata is required when "
            "activity is enabled; regenerate the bundle from the original "
            "KinaseActivityResult"
        )


def _reject_disabled_activity_metadata(
    *,
    activity_method_metadata: Mapping[str, object] | None,
    activity_method_summary: Mapping[str, object] | None,
    activity_input_semantics: Mapping[str, object] | None,
    activity_profile_metadata: Mapping[str, object] | None,
    activity_membership_selection: Mapping[str, object] | None,
) -> None:
    if activity_method_metadata is not None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.method must be null when activity is disabled"
        )
    if activity_method_summary is not None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.summary must be null when activity is disabled"
        )
    if activity_input_semantics is not None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.input_semantics must be null when "
            "activity is disabled; remove the semantic payload or regenerate the "
            "bundle"
        )
    if activity_profile_metadata is not None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.profile_metadata must be null when "
            "activity is disabled; remove the semantic payload or regenerate the "
            "bundle"
        )
    if activity_membership_selection is not None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.membership_selection must be null "
            "when activity is disabled; remove the membership payload or "
            "regenerate the bundle"
        )


def _parse_table_sections(
    *,
    top_level: _KinaseManifestTopLevel,
    outputs: _KinaseOutputPayloads,
) -> _KinaseTableSections:
    return _KinaseTableSections(
        dataset_tables=_parse_dataset_tables(top_level.dataset_payload),
        reference_tables=_parse_reference_tables(top_level.references_payload),
        scoring_tables=_parse_scoring_tables(outputs.scoring_payload),
        prediction_tables=_parse_prediction_tables(outputs.prediction_payload),
        activity_tables=_parse_activity_tables(outputs.activity_payload),
    )


def _parse_dataset_tables(
    dataset_payload: Mapping[str, object],
) -> Mapping[str, object]:
    dataset_tables = require_mapping(
        dataset_payload.get("tables"),
        field_name="bundle manifest.dataset.tables",
    )
    _reject_unsupported_fields(
        dataset_tables,
        field_name="bundle manifest.dataset.tables",
        allowed_fields=_DATASET_TABLE_KEYS,
    )
    _require_fields(
        dataset_tables,
        field_name="bundle manifest.dataset.tables",
        required_fields=_DATASET_TABLE_KEYS,
        unsupported_shape=True,
    )
    return dataset_tables


def _parse_reference_tables(
    references_payload: Mapping[str, object],
) -> Mapping[str, object]:
    reference_tables = require_mapping(
        references_payload.get("tables"),
        field_name="bundle manifest.resolved_references.tables",
    )
    _reject_unsupported_fields(
        reference_tables,
        field_name="bundle manifest.resolved_references.tables",
        allowed_fields=_REFERENCE_TABLE_KEYS,
    )
    _require_fields(
        reference_tables,
        field_name="bundle manifest.resolved_references.tables",
        required_fields=_REFERENCE_TABLE_KEYS,
        unsupported_shape=True,
    )
    return reference_tables


def _parse_scoring_tables(
    scoring_payload: Mapping[str, object],
) -> Mapping[str, object]:
    scoring_tables = require_mapping(
        scoring_payload.get("tables"),
        field_name="bundle manifest.outputs.scoring.tables",
    )
    _reject_unsupported_fields(
        scoring_tables,
        field_name="bundle manifest.outputs.scoring.tables",
        allowed_fields=_SCORING_TABLE_KEYS,
    )
    _require_fields(
        scoring_tables,
        field_name="bundle manifest.outputs.scoring.tables",
        required_fields=_SCORING_TABLE_REQUIRED_KEYS,
        unsupported_shape=True,
    )
    return scoring_tables


def _parse_prediction_tables(
    prediction_payload: Mapping[str, object],
) -> Mapping[str, object]:
    prediction_tables = require_mapping(
        prediction_payload.get("tables"),
        field_name="bundle manifest.outputs.prediction.tables",
    )
    _reject_unsupported_fields(
        prediction_tables,
        field_name="bundle manifest.outputs.prediction.tables",
        allowed_fields=_PREDICTION_TABLE_KEYS,
    )
    _require_fields(
        prediction_tables,
        field_name="bundle manifest.outputs.prediction.tables",
        required_fields=_PREDICTION_TABLE_KEYS,
        unsupported_shape=True,
    )
    return prediction_tables


def _parse_activity_tables(
    activity_payload: Mapping[str, object],
) -> Mapping[str, object]:
    activity_tables = require_mapping(
        activity_payload.get("tables"),
        field_name="bundle manifest.outputs.activity.tables",
    )
    _reject_unsupported_fields(
        activity_tables,
        field_name="bundle manifest.outputs.activity.tables",
        allowed_fields=_ACTIVITY_TABLE_KEYS,
    )
    _require_fields(
        activity_tables,
        field_name="bundle manifest.outputs.activity.tables",
        required_fields=_ACTIVITY_TABLE_KEYS,
        unsupported_shape=True,
    )
    return activity_tables


def _validate_manifest_paths(schema: _KinaseManifestSchema) -> None:
    """Validate manifest-declared path strings before file-record validation."""

    _validate_table_paths(
        schema.dataset_tables,
        field_name="bundle manifest.dataset.tables",
    )
    _validate_table_paths(
        schema.reference_tables,
        field_name="bundle manifest.resolved_references.tables",
    )
    _validate_table_paths(
        schema.scoring_tables,
        field_name="bundle manifest.outputs.scoring.tables",
    )
    _validate_table_paths(
        schema.prediction_tables,
        field_name="bundle manifest.outputs.prediction.tables",
    )
    _validate_table_paths(
        schema.activity_tables,
        field_name="bundle manifest.outputs.activity.tables",
    )
    _validate_file_path(
        schema.config_snapshot_raw,
        field_name="bundle manifest.config_snapshot",
    )


def _validate_manifest_file_records(
    schema: _KinaseManifestSchema,
) -> _KinaseManifestFileRecords:
    """Validate manifest table/file records after path validation."""

    _require_table_entries(
        schema.dataset_tables,
        field_name="bundle manifest.dataset.tables",
        optional_keys=frozenset({"sample_metadata", "total"}),
    )
    _require_table_entries(
        schema.reference_tables,
        field_name="bundle manifest.resolved_references.tables",
        optional_keys=frozenset(),
    )
    _require_table_entries(
        schema.scoring_tables,
        field_name="bundle manifest.outputs.scoring.tables",
        optional_keys=frozenset(
            {
                "motif_scores",
                "rank_weighted_fusion_scores",
                "kinase_library_motif_scores",
                "combined_profile_motif_scores",
                "score_fusion_weights",
                "kinase_library_site_diagnostics",
                "kinase_library_kinase_diagnostics",
                "substrate_contributions",
            }
        ),
    )
    _require_table_entries(
        schema.prediction_tables,
        field_name="bundle manifest.outputs.prediction.tables",
        optional_keys=frozenset({"substrate_list"}),
    )
    _require_table_entries(
        schema.activity_tables,
        field_name="bundle manifest.outputs.activity.tables",
        optional_keys=_ACTIVITY_TABLE_KEYS,
    )
    try:
        config_snapshot_entry = require_file_entry(
            schema.config_snapshot_raw,
            field_name="bundle manifest.config_snapshot",
            expected_logical_type="config_snapshot",
        )
    except PhosPyInputError as exc:
        _raise_unsupported_manifest_shape(str(exc))
    return _KinaseManifestFileRecords(
        dataset_tables=schema.dataset_tables,
        reference_tables=schema.reference_tables,
        scoring_tables=schema.scoring_tables,
        prediction_tables=schema.prediction_tables,
        activity_tables=schema.activity_tables,
        config_snapshot_entry=config_snapshot_entry,
    )


def _assemble_manifest_sections(
    *,
    schema: _KinaseManifestSchema,
    file_records: _KinaseManifestFileRecords,
) -> KinaseManifestSections:
    """Assemble the public decoded manifest section model."""

    return KinaseManifestSections(
        manifest_version=schema.manifest_version,
        dataset_metadata=schema.dataset_metadata,
        dataset_tables=file_records.dataset_tables,
        references_metadata=schema.references_metadata,
        reference_tables=file_records.reference_tables,
        scoring_tables=file_records.scoring_tables,
        prediction_tables=file_records.prediction_tables,
        activity_enabled=schema.activity_enabled,
        activity_method_metadata=schema.activity_method_metadata,
        activity_method_summary=schema.activity_method_summary,
        activity_input_semantics=schema.activity_input_semantics,
        activity_profile_metadata=schema.activity_profile_metadata,
        activity_membership_selection=schema.activity_membership_selection,
        activity_tables=file_records.activity_tables,
        provenance_payload=schema.provenance_payload,
        config_snapshot_entry=file_records.config_snapshot_entry,
        caveats_payload=schema.caveats_payload,
    )


def _validate_table_paths(
    tables: Mapping[str, object],
    *,
    field_name: str,
) -> None:
    for table_key, value in tables.items():
        if value is None:
            continue
        entry_field_name = f"{field_name}.{str(table_key)}"
        _validate_file_path(value, field_name=entry_field_name)


def _validate_file_path(value: object, *, field_name: str) -> None:
    try:
        entry = require_mapping(value, field_name=field_name)
        relative_path = require_str(entry.get("path"), field_name=f"{field_name}.path")
        _validate_relative_manifest_path(
            relative_path,
            field_name=f"{field_name}.path",
        )
    except PhosPyInputError as exc:
        _raise_unsupported_manifest_shape(str(exc))


def _validate_relative_manifest_path(path: str, *, field_name: str) -> None:
    candidate = PurePosixPath(path)
    if candidate.is_absolute():
        _raise_unsupported_manifest_shape(f"{field_name} must be a relative path")
    if ".." in candidate.parts:
        _raise_unsupported_manifest_shape(
            f"{field_name} must not contain parent-directory traversal"
        )


def _raise_unsupported_manifest_shape(message: str) -> NoReturn:
    raise PhosPyInputError(f"{_LEGACY_KINASE_BUNDLE_SCHEMA_ERROR} {message}.")


def _require_table_entries(
    tables: Mapping[str, object],
    *,
    field_name: str,
    optional_keys: frozenset[str],
) -> None:
    for table_key, value in tables.items():
        entry_field_name = f"{field_name}.{str(table_key)}"
        try:
            if str(table_key) in optional_keys:
                require_optional_table_entry(value, field_name=entry_field_name)
            else:
                require_table_entry(value, field_name=entry_field_name)
        except PhosPyInputError as exc:
            _raise_unsupported_manifest_shape(str(exc))


def _require_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    required_fields: frozenset[str],
    unsupported_shape: bool,
) -> None:
    missing_fields = sorted(
        str(key) for key in required_fields if str(key) not in payload
    )
    if not missing_fields:
        return
    missing = ", ".join(missing_fields)
    message = f"{field_name} is missing required field(s): {missing}"
    if unsupported_shape:
        _raise_unsupported_manifest_shape(message)
    raise PhosPyInputError(message)


def _reject_unsupported_fields(
    payload: Mapping[str, object],
    *,
    field_name: str,
    allowed_fields: frozenset[str],
) -> None:
    unknown_fields = sorted(
        str(key) for key in payload.keys() if str(key) not in allowed_fields
    )
    if unknown_fields:
        unknown = ", ".join(unknown_fields)
        _raise_unsupported_manifest_shape(
            f"{field_name} contains unsupported field(s): {unknown}"
        )
