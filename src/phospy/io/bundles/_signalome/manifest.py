"""Current-contract manifest serialization and parsing for signalome bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import NoReturn

from phospy.contracts.results import SignalomeWorkflowResult
from phospy.errors.input import PhosPyInputError
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
from phospy.io.bundles._signalome.constants import (
    SIGNALOME_BUNDLE_KIND,
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
)
from phospy.io.bundles._signalome.diagnostics import (
    signalome_alignment_diagnostics_to_payload,
    signalome_clustering_preparation_diagnostics_to_payload,
    signalome_module_selection_diagnostics_to_payload,
    signalome_network_correlation_diagnostics_to_payload,
    signalome_score_preconditioning_diagnostics_to_payload,
)
from phospy.provenance.serialization import to_payload as provenance_to_payload


@dataclass(frozen=True, slots=True)
class SignalomeManifestSections:
    """Decoded current manifest sections needed to load signalome bundles."""

    manifest_version: int
    dataset_metadata: Mapping[str, object]
    dataset_tables: Mapping[str, object]
    references_metadata: Mapping[str, object]
    reference_tables: Mapping[str, object]
    scoring_tables: Mapping[str, object]
    prediction_tables: Mapping[str, object]
    upstream_activity_enabled: bool
    activity_tables: Mapping[str, object]
    signalome_metadata: Mapping[str, object]
    signalome_tables: Mapping[str, object]
    provenance_payload: Mapping[str, object]
    config_snapshot_entry: Mapping[str, object]
    signalome_caveats_payload: object = ()
    upstream_kinase_caveats_payload: object = ()


@dataclass(frozen=True, slots=True)
class _SignalomeManifestSchema:
    """Schema-parsed signalome manifest payload before file-record validation."""

    manifest_version: int
    dataset_payload: Mapping[str, object]
    references_payload: Mapping[str, object]
    upstream_payload: Mapping[str, object]
    scoring_payload: Mapping[str, object]
    prediction_payload: Mapping[str, object]
    activity_payload: Mapping[str, object]
    signalome_outputs_payload: Mapping[str, object]
    dataset_tables: Mapping[str, object]
    reference_tables: Mapping[str, object]
    scoring_tables: Mapping[str, object]
    prediction_tables: Mapping[str, object]
    upstream_activity_enabled: bool
    activity_tables: Mapping[str, object]
    signalome_metadata: Mapping[str, object]
    signalome_tables: Mapping[str, object]
    provenance_payload: Mapping[str, object]
    config_snapshot_raw: object
    signalome_caveats_payload: object
    upstream_kinase_caveats_payload: object


@dataclass(frozen=True, slots=True)
class _SignalomeManifestFileRecords:
    """Signalome manifest table/file records after validation."""

    dataset_tables: Mapping[str, object]
    reference_tables: Mapping[str, object]
    scoring_tables: Mapping[str, object]
    prediction_tables: Mapping[str, object]
    activity_tables: Mapping[str, object]
    signalome_tables: Mapping[str, object]
    config_snapshot_entry: Mapping[str, object]


_LEGACY_SIGNALOME_BUNDLE_SCHEMA_ERROR = (
    "Legacy signalome bundle schemas are no longer supported. Regenerate the bundle "
    "with the current PhosPy version."
)
_MANIFEST_ALLOWED_FIELDS = frozenset(
    {
        "bundle_type",
        "manifest_version",
        "table_format",
        "dataset",
        "resolved_references",
        "upstream_kinase_outputs",
        "signalome_outputs",
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
        "upstream_kinase_outputs",
        "signalome_outputs",
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
_UPSTREAM_OUTPUTS_ALLOWED_FIELDS = frozenset(
    {"scoring", "prediction", "activity", "caveats"}
)
_UPSTREAM_OUTPUTS_REQUIRED_FIELDS = frozenset({"scoring", "prediction", "activity"})
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
_ACTIVITY_ALLOWED_FIELDS = frozenset({"enabled", "tables"})
_ACTIVITY_TABLE_KEYS = frozenset(
    {
        "weighted_activity",
        "thresholded_substrate_mean_activity",
        "thresholded_substrate_counts",
        "activity_substrate_counts",
        "target_counts",
        "target_table",
    }
)
_SIGNALOME_OUTPUTS_ALLOWED_FIELDS = frozenset({"metadata", "tables"})
_SIGNALOME_METADATA_ALLOWED_KEYS = frozenset(
    {
        "kinase_network_nodes_present",
        "expanded_signalome_present",
        "module_selection_diagnostics",
        "clustering_preparation_diagnostics",
        "score_preconditioning_diagnostics",
        "network_correlation_diagnostics",
        "alignment_diagnostics",
    }
)
_SIGNALOME_METADATA_REQUIRED_KEYS = frozenset(
    {
        "kinase_network_nodes_present",
        "expanded_signalome_present",
        "module_selection_diagnostics",
        "score_preconditioning_diagnostics",
        "network_correlation_diagnostics",
    }
)
_SIGNALOME_TABLE_KEYS = frozenset(
    {
        "module_assignments",
        "signalome_modules",
        "kinase_network_edges",
        "kinase_network_nodes",
        "kinase_network_candidate_correlations",
        "expanded_signalome",
        "site_membership",
        "protein_site_context",
    }
)


def build_manifest(
    *,
    result: SignalomeWorkflowResult,
    table_format: str,
    dataset_tables: Mapping[str, object],
    reference_tables: Mapping[str, object],
    scoring_tables: Mapping[str, object],
    prediction_tables: Mapping[str, object],
    activity_tables: Mapping[str, object],
    signalome_tables: Mapping[str, object],
    config_snapshot_entry: Mapping[str, object],
) -> dict[str, object]:
    """Build the current signalome manifest payload."""

    if result.provenance is None:
        raise PhosPyInputError(
            "signalome bundle saving requires result.provenance; "
            "bundle manifests must include explicit provenance metadata"
        )

    return {
        "bundle_type": SIGNALOME_BUNDLE_KIND,
        "manifest_version": SIGNALOME_BUNDLE_MANIFEST_VERSION,
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
                "organism": result.kinase_result.references.organism.value,
            },
            "tables": dict(reference_tables),
        },
        "upstream_kinase_outputs": {
            "scoring": {
                "tables": dict(scoring_tables),
            },
            "prediction": {
                "tables": dict(prediction_tables),
            },
            "activity": {
                "enabled": result.kinase_result.activity_result is not None,
                "tables": dict(activity_tables),
            },
            "caveats": [caveat.to_payload() for caveat in result.kinase_result.caveats],
        },
        "signalome_outputs": {
            "metadata": {
                "kinase_network_nodes_present": result.kinase_network.nodes is not None,
                "expanded_signalome_present": result.expanded_signalome is not None,
                "module_selection_diagnostics": signalome_module_selection_diagnostics_to_payload(
                    result.module_selection_diagnostics
                ),
                "clustering_preparation_diagnostics": signalome_clustering_preparation_diagnostics_to_payload(
                    result.clustering_preparation_diagnostics
                ),
                "score_preconditioning_diagnostics": signalome_score_preconditioning_diagnostics_to_payload(
                    result.score_preconditioning_diagnostics
                ),
                "alignment_diagnostics": signalome_alignment_diagnostics_to_payload(
                    result.alignment_diagnostics
                ),
                "network_correlation_diagnostics": signalome_network_correlation_diagnostics_to_payload(
                    result.kinase_network.correlation_diagnostics
                ),
            },
            "tables": dict(signalome_tables),
        },
        "caveats": [caveat.to_payload() for caveat in result.caveats],
        "provenance": provenance_to_payload(result.provenance),
        "config_snapshot": dict(config_snapshot_entry),
    }


def parse_manifest(payload: Mapping[str, object]) -> SignalomeManifestSections:
    """Parse and validate current-contract signalome manifest payload."""

    schema = _parse_manifest_schema(payload)
    _validate_manifest_paths(schema)
    file_records = _validate_manifest_file_records(schema)
    return _assemble_manifest_sections(schema=schema, file_records=file_records)


def _parse_manifest_schema(payload: Mapping[str, object]) -> _SignalomeManifestSchema:
    """Parse manifest object schemas and semantic switches."""

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
    if bundle_type != SIGNALOME_BUNDLE_KIND:
        _raise_unsupported_manifest_shape(
            "unsupported bundle manifest bundle_type "
            f"'{bundle_type}'; expected '{SIGNALOME_BUNDLE_KIND}'"
        )
    manifest_version = require_int(
        payload.get("manifest_version"),
        field_name="bundle manifest.manifest_version",
    )
    if manifest_version != SIGNALOME_BUNDLE_MANIFEST_VERSION:
        _raise_unsupported_manifest_shape(
            "unsupported bundle manifest version "
            f"'{manifest_version}'; expected {SIGNALOME_BUNDLE_MANIFEST_VERSION}"
        )
    require_str(
        payload.get("table_format"),
        field_name="bundle manifest.table_format",
    )

    dataset_payload = require_mapping(
        payload.get("dataset"),
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
    references_payload = require_mapping(
        payload.get("resolved_references"),
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
    upstream_payload = require_mapping(
        payload.get("upstream_kinase_outputs"),
        field_name="bundle manifest.upstream_kinase_outputs",
    )
    _reject_unsupported_fields(
        upstream_payload,
        field_name="bundle manifest.upstream_kinase_outputs",
        allowed_fields=_UPSTREAM_OUTPUTS_ALLOWED_FIELDS,
    )
    _require_fields(
        upstream_payload,
        field_name="bundle manifest.upstream_kinase_outputs",
        required_fields=_UPSTREAM_OUTPUTS_REQUIRED_FIELDS,
        unsupported_shape=True,
    )
    scoring_payload = require_mapping(
        upstream_payload.get("scoring"),
        field_name="bundle manifest.upstream_kinase_outputs.scoring",
    )
    _reject_unsupported_fields(
        scoring_payload,
        field_name="bundle manifest.upstream_kinase_outputs.scoring",
        allowed_fields=_SCORING_ALLOWED_FIELDS,
    )
    _require_fields(
        scoring_payload,
        field_name="bundle manifest.upstream_kinase_outputs.scoring",
        required_fields=_SCORING_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    prediction_payload = require_mapping(
        upstream_payload.get("prediction"),
        field_name="bundle manifest.upstream_kinase_outputs.prediction",
    )
    _reject_unsupported_fields(
        prediction_payload,
        field_name="bundle manifest.upstream_kinase_outputs.prediction",
        allowed_fields=_PREDICTION_ALLOWED_FIELDS,
    )
    _require_fields(
        prediction_payload,
        field_name="bundle manifest.upstream_kinase_outputs.prediction",
        required_fields=_PREDICTION_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    activity_payload = require_mapping(
        upstream_payload.get("activity"),
        field_name="bundle manifest.upstream_kinase_outputs.activity",
    )
    _reject_unsupported_fields(
        activity_payload,
        field_name="bundle manifest.upstream_kinase_outputs.activity",
        allowed_fields=_ACTIVITY_ALLOWED_FIELDS,
    )
    _require_fields(
        activity_payload,
        field_name="bundle manifest.upstream_kinase_outputs.activity",
        required_fields=_ACTIVITY_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    upstream_activity_enabled = require_bool(
        activity_payload.get("enabled"),
        field_name="bundle manifest.upstream_kinase_outputs.activity.enabled",
    )
    signalome_outputs_payload = require_mapping(
        payload.get("signalome_outputs"),
        field_name="bundle manifest.signalome_outputs",
    )
    _reject_unsupported_fields(
        signalome_outputs_payload,
        field_name="bundle manifest.signalome_outputs",
        allowed_fields=_SIGNALOME_OUTPUTS_ALLOWED_FIELDS,
    )
    _require_fields(
        signalome_outputs_payload,
        field_name="bundle manifest.signalome_outputs",
        required_fields=_SIGNALOME_OUTPUTS_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    if payload.get("provenance") is None:
        _raise_unsupported_manifest_shape("bundle manifest.provenance is required")
    provenance_payload = require_mapping(
        payload.get("provenance"),
        field_name="bundle manifest.provenance",
    )
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
    scoring_tables = require_mapping(
        scoring_payload.get("tables"),
        field_name="bundle manifest.upstream_kinase_outputs.scoring.tables",
    )
    _reject_unsupported_fields(
        scoring_tables,
        field_name="bundle manifest.upstream_kinase_outputs.scoring.tables",
        allowed_fields=_SCORING_TABLE_KEYS,
    )
    _require_fields(
        scoring_tables,
        field_name="bundle manifest.upstream_kinase_outputs.scoring.tables",
        required_fields=_SCORING_TABLE_REQUIRED_KEYS,
        unsupported_shape=True,
    )
    prediction_tables = require_mapping(
        prediction_payload.get("tables"),
        field_name="bundle manifest.upstream_kinase_outputs.prediction.tables",
    )
    _reject_unsupported_fields(
        prediction_tables,
        field_name="bundle manifest.upstream_kinase_outputs.prediction.tables",
        allowed_fields=_PREDICTION_TABLE_KEYS,
    )
    _require_fields(
        prediction_tables,
        field_name="bundle manifest.upstream_kinase_outputs.prediction.tables",
        required_fields=_PREDICTION_TABLE_KEYS,
        unsupported_shape=True,
    )
    activity_tables = require_mapping(
        activity_payload.get("tables"),
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables",
    )
    _reject_unsupported_fields(
        activity_tables,
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables",
        allowed_fields=_ACTIVITY_TABLE_KEYS,
    )
    _require_fields(
        activity_tables,
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables",
        required_fields=_ACTIVITY_TABLE_KEYS,
        unsupported_shape=True,
    )
    signalome_metadata = require_mapping(
        signalome_outputs_payload.get("metadata"),
        field_name="bundle manifest.signalome_outputs.metadata",
    )
    _reject_unsupported_fields(
        signalome_metadata,
        field_name="bundle manifest.signalome_outputs.metadata",
        allowed_fields=_SIGNALOME_METADATA_ALLOWED_KEYS,
    )
    _require_fields(
        signalome_metadata,
        field_name="bundle manifest.signalome_outputs.metadata",
        required_fields=_SIGNALOME_METADATA_REQUIRED_KEYS,
        unsupported_shape=True,
    )
    require_bool(
        signalome_metadata.get("kinase_network_nodes_present"),
        field_name="bundle manifest.signalome_outputs.metadata.kinase_network_nodes_present",
    )
    require_bool(
        signalome_metadata.get("expanded_signalome_present"),
        field_name="bundle manifest.signalome_outputs.metadata.expanded_signalome_present",
    )
    signalome_tables = require_mapping(
        signalome_outputs_payload.get("tables"),
        field_name="bundle manifest.signalome_outputs.tables",
    )
    _reject_unsupported_fields(
        signalome_tables,
        field_name="bundle manifest.signalome_outputs.tables",
        allowed_fields=_SIGNALOME_TABLE_KEYS,
    )
    _require_fields(
        signalome_tables,
        field_name="bundle manifest.signalome_outputs.tables",
        required_fields=_SIGNALOME_TABLE_KEYS,
        unsupported_shape=True,
    )
    return _SignalomeManifestSchema(
        manifest_version=manifest_version,
        dataset_payload=dataset_payload,
        references_payload=references_payload,
        upstream_payload=upstream_payload,
        scoring_payload=scoring_payload,
        prediction_payload=prediction_payload,
        activity_payload=activity_payload,
        signalome_outputs_payload=signalome_outputs_payload,
        dataset_tables=dataset_tables,
        reference_tables=reference_tables,
        scoring_tables=scoring_tables,
        prediction_tables=prediction_tables,
        upstream_activity_enabled=upstream_activity_enabled,
        activity_tables=activity_tables,
        signalome_metadata=signalome_metadata,
        signalome_tables=signalome_tables,
        provenance_payload=provenance_payload,
        config_snapshot_raw=payload.get("config_snapshot"),
        signalome_caveats_payload=payload.get("caveats", []),
        upstream_kinase_caveats_payload=upstream_payload.get("caveats", []),
    )


def _validate_manifest_paths(schema: _SignalomeManifestSchema) -> None:
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
        field_name="bundle manifest.upstream_kinase_outputs.scoring.tables",
    )
    _validate_table_paths(
        schema.prediction_tables,
        field_name="bundle manifest.upstream_kinase_outputs.prediction.tables",
    )
    _validate_table_paths(
        schema.activity_tables,
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables",
    )
    _validate_table_paths(
        schema.signalome_tables,
        field_name="bundle manifest.signalome_outputs.tables",
    )
    _validate_file_path(
        schema.config_snapshot_raw,
        field_name="bundle manifest.config_snapshot",
    )


def _validate_manifest_file_records(
    schema: _SignalomeManifestSchema,
) -> _SignalomeManifestFileRecords:
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
        field_name="bundle manifest.upstream_kinase_outputs.scoring.tables",
        optional_keys=frozenset(
            {
                "motif_scores",
                "rank_weighted_fusion_scores",
                "kinase_library_motif_scores",
                "combined_profile_motif_scores",
                "score_fusion_weights",
                "kinase_library_site_diagnostics",
                "kinase_library_kinase_diagnostics",
            }
        ),
    )
    _require_table_entries(
        schema.prediction_tables,
        field_name="bundle manifest.upstream_kinase_outputs.prediction.tables",
        optional_keys=frozenset({"substrate_list"}),
    )
    _require_table_entries(
        schema.activity_tables,
        field_name="bundle manifest.upstream_kinase_outputs.activity.tables",
        optional_keys=_ACTIVITY_TABLE_KEYS,
    )
    _require_table_entries(
        schema.signalome_tables,
        field_name="bundle manifest.signalome_outputs.tables",
        optional_keys=frozenset(
            {
                "kinase_network_nodes",
                "kinase_network_candidate_correlations",
                "expanded_signalome",
                "site_membership",
                "protein_site_context",
            }
        ),
    )
    try:
        config_snapshot_entry = require_file_entry(
            schema.config_snapshot_raw,
            field_name="bundle manifest.config_snapshot",
            expected_logical_type="config_snapshot",
        )
    except PhosPyInputError as exc:
        _raise_unsupported_manifest_shape(str(exc))
    return _SignalomeManifestFileRecords(
        dataset_tables=schema.dataset_tables,
        reference_tables=schema.reference_tables,
        scoring_tables=schema.scoring_tables,
        prediction_tables=schema.prediction_tables,
        activity_tables=schema.activity_tables,
        signalome_tables=schema.signalome_tables,
        config_snapshot_entry=config_snapshot_entry,
    )


def _assemble_manifest_sections(
    *,
    schema: _SignalomeManifestSchema,
    file_records: _SignalomeManifestFileRecords,
) -> SignalomeManifestSections:
    """Assemble the public decoded manifest section model."""

    return SignalomeManifestSections(
        manifest_version=schema.manifest_version,
        dataset_metadata=require_mapping(
            schema.dataset_payload.get("metadata"),
            field_name="bundle manifest.dataset.metadata",
        ),
        dataset_tables=file_records.dataset_tables,
        references_metadata=require_mapping(
            schema.references_payload.get("metadata"),
            field_name="bundle manifest.resolved_references.metadata",
        ),
        reference_tables=file_records.reference_tables,
        scoring_tables=file_records.scoring_tables,
        prediction_tables=file_records.prediction_tables,
        upstream_activity_enabled=schema.upstream_activity_enabled,
        activity_tables=file_records.activity_tables,
        signalome_metadata=schema.signalome_metadata,
        signalome_tables=file_records.signalome_tables,
        provenance_payload=schema.provenance_payload,
        config_snapshot_entry=file_records.config_snapshot_entry,
        signalome_caveats_payload=schema.signalome_caveats_payload,
        upstream_kinase_caveats_payload=schema.upstream_kinase_caveats_payload,
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
    raise PhosPyInputError(f"{_LEGACY_SIGNALOME_BUNDLE_SCHEMA_ERROR} {message}.")


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
