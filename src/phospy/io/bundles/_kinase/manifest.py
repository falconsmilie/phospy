"""Current-contract manifest serialization and parsing for kinase bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
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
    activity_tables: Mapping[str, object]
    provenance_payload: Mapping[str, object]
    config_snapshot_entry: Mapping[str, object]
    caveats_payload: object = ()


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
_ACTIVITY_ALLOWED_FIELDS = frozenset({"enabled", "method", "summary", "tables"})
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
                "tables": dict(activity_tables),
            },
        },
        "provenance": provenance_to_payload(result.provenance),
        "caveats": [caveat.to_payload() for caveat in result.caveats],
        "config_snapshot": dict(config_snapshot_entry),
    }


def parse_manifest(payload: Mapping[str, object]) -> KinaseManifestSections:
    """Parse and validate current-contract kinase manifest payload."""

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
        _raise_unsupported_manifest_shape(
            "unsupported bundle manifest version "
            f"'{manifest_version}'; expected {KINASE_BUNDLE_MANIFEST_VERSION}"
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
    outputs_payload = require_mapping(
        payload.get("outputs"),
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
        required_fields=_ACTIVITY_ALLOWED_FIELDS,
        unsupported_shape=True,
    )
    activity_enabled = require_bool(
        activity_payload.get("enabled"),
        field_name="bundle manifest.outputs.activity.enabled",
    )
    activity_method_raw = activity_payload.get("method")
    activity_method_metadata = (
        None
        if activity_method_raw is None
        else require_mapping(
            activity_method_raw,
            field_name="bundle manifest.outputs.activity.method",
        )
    )
    activity_summary_raw = activity_payload.get("summary")
    activity_method_summary = (
        None
        if activity_summary_raw is None
        else require_mapping(
            activity_summary_raw,
            field_name="bundle manifest.outputs.activity.summary",
        )
    )
    if activity_enabled and activity_method_metadata is None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.method is required when activity is enabled"
        )
    if not activity_enabled and activity_method_metadata is not None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.method must be null when activity is disabled"
        )
    if not activity_enabled and activity_method_summary is not None:
        _raise_unsupported_manifest_shape(
            "bundle manifest.outputs.activity.summary must be null when activity is disabled"
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
    _require_table_entries(
        dataset_tables,
        field_name="bundle manifest.dataset.tables",
        optional_keys=frozenset({"sample_metadata", "total"}),
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
    _require_table_entries(
        reference_tables,
        field_name="bundle manifest.resolved_references.tables",
        optional_keys=frozenset(),
    )
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
    _require_table_entries(
        scoring_tables,
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
    _require_table_entries(
        prediction_tables,
        field_name="bundle manifest.outputs.prediction.tables",
        optional_keys=frozenset({"substrate_list"}),
    )
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
    _require_table_entries(
        activity_tables,
        field_name="bundle manifest.outputs.activity.tables",
        optional_keys=_ACTIVITY_TABLE_KEYS,
    )
    try:
        config_snapshot_entry = require_file_entry(
            payload.get("config_snapshot"),
            field_name="bundle manifest.config_snapshot",
            expected_logical_type="config_snapshot",
        )
    except PhosPyInputError as exc:
        _raise_unsupported_manifest_shape(str(exc))

    return KinaseManifestSections(
        manifest_version=manifest_version,
        dataset_metadata=require_mapping(
            dataset_payload.get("metadata"),
            field_name="bundle manifest.dataset.metadata",
        ),
        dataset_tables=dataset_tables,
        references_metadata=require_mapping(
            references_payload.get("metadata"),
            field_name="bundle manifest.resolved_references.metadata",
        ),
        reference_tables=reference_tables,
        scoring_tables=scoring_tables,
        prediction_tables=prediction_tables,
        activity_enabled=activity_enabled,
        activity_method_metadata=activity_method_metadata,
        activity_method_summary=activity_method_summary,
        activity_tables=activity_tables,
        provenance_payload=provenance_payload,
        config_snapshot_entry=config_snapshot_entry,
        caveats_payload=payload.get("caveats", []),
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
