"""Current-contract manifest serialization and parsing for signalome bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.contracts.results import SignalomeWorkflowResult
from phospy.errors.input import PhosPyInputError
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
    CONFIG_SNAPSHOT_RELATIVE_PATH,
    SIGNALOME_BUNDLE_KIND,
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
)
from phospy.io.bundles._signalome.diagnostics import (
    signalome_alignment_diagnostics_to_payload,
    signalome_module_selection_diagnostics_to_payload,
    signalome_network_correlation_diagnostics_to_payload,
    signalome_score_preconditioning_diagnostics_to_payload,
)
from phospy.provenance.serialization import to_payload as provenance_to_payload


@dataclass(frozen=True, slots=True)
class SignalomeManifestSections:
    """Decoded v1 manifest sections needed to load signalome bundles."""

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
    config_snapshot_path: str
    signalome_caveats_payload: object = ()
    upstream_kinase_caveats_payload: object = ()


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
) -> dict[str, object]:
    """Build the v1 signalome manifest payload."""

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
        "config_snapshot": CONFIG_SNAPSHOT_RELATIVE_PATH,
    }


def parse_manifest(payload: Mapping[str, object]) -> SignalomeManifestSections:
    """Parse and validate current-contract signalome manifest payload."""

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

    return SignalomeManifestSections(
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
        upstream_activity_enabled=upstream_activity_enabled,
        activity_tables=activity_tables,
        signalome_metadata=signalome_metadata,
        signalome_tables=signalome_tables,
        provenance_payload=provenance_payload,
        config_snapshot_path=require_str(
            payload.get("config_snapshot"),
            field_name="bundle manifest.config_snapshot",
        ),
        signalome_caveats_payload=payload.get("caveats", []),
        upstream_kinase_caveats_payload=upstream_payload.get("caveats", []),
    )


def _raise_unsupported_manifest_shape(message: str) -> None:
    raise PhosPyInputError(f"{_LEGACY_SIGNALOME_BUNDLE_SCHEMA_ERROR} {message}.")


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
