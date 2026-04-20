"""Current-contract manifest serialization and parsing for signalome bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.api.results import SignalomeWorkflowResult
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._shared.primitives import (
    require_int,
    require_mapping,
    require_str,
)
from phospy.io.bundles._shared.transformation_state import (
    transformation_state_to_payload,
)
from phospy.io.bundles._signalome.compatibility import (
    signalome_module_selection_diagnostics_to_payload,
    signalome_score_preconditioning_diagnostics_to_payload,
)
from phospy.io.bundles._signalome.constants import (
    CONFIG_SNAPSHOT_RELATIVE_PATH,
    SIGNALOME_BUNDLE_KIND,
    SIGNALOME_BUNDLE_MANIFEST_VERSION,
)


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
    activity_tables: Mapping[str, object]
    signalome_metadata: Mapping[str, object]
    signalome_tables: Mapping[str, object]
    config_snapshot_path: str


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
                "transformation_state": transformation_state_to_payload(
                    result.dataset.transformation_state
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
            },
            "tables": dict(signalome_tables),
        },
        "config_snapshot": CONFIG_SNAPSHOT_RELATIVE_PATH,
    }


def parse_manifest(payload: Mapping[str, object]) -> SignalomeManifestSections:
    """Parse and validate current-contract signalome manifest payload."""

    bundle_type = require_str(
        payload.get("bundle_type"),
        field_name="bundle manifest.bundle_type",
    )
    if bundle_type != SIGNALOME_BUNDLE_KIND:
        raise PhosPyInputError(
            "unsupported bundle manifest bundle_type "
            f"'{bundle_type}'; expected '{SIGNALOME_BUNDLE_KIND}'"
        )
    manifest_version = require_int(
        payload.get("manifest_version"),
        field_name="bundle manifest.manifest_version",
    )
    if manifest_version != SIGNALOME_BUNDLE_MANIFEST_VERSION:
        raise PhosPyInputError(
            "unsupported bundle manifest version "
            f"'{manifest_version}'; expected {SIGNALOME_BUNDLE_MANIFEST_VERSION}"
        )

    dataset_payload = require_mapping(
        payload.get("dataset"),
        field_name="bundle manifest.dataset",
    )
    references_payload = require_mapping(
        payload.get("resolved_references"),
        field_name="bundle manifest.resolved_references",
    )
    upstream_payload = require_mapping(
        payload.get("upstream_kinase_outputs"),
        field_name="bundle manifest.upstream_kinase_outputs",
    )
    scoring_payload = require_mapping(
        upstream_payload.get("scoring"),
        field_name="bundle manifest.upstream_kinase_outputs.scoring",
    )
    prediction_payload = require_mapping(
        upstream_payload.get("prediction"),
        field_name="bundle manifest.upstream_kinase_outputs.prediction",
    )
    activity_payload = require_mapping(
        upstream_payload.get("activity"),
        field_name="bundle manifest.upstream_kinase_outputs.activity",
    )
    signalome_outputs_payload = require_mapping(
        payload.get("signalome_outputs"),
        field_name="bundle manifest.signalome_outputs",
    )

    return SignalomeManifestSections(
        manifest_version=manifest_version,
        dataset_metadata=require_mapping(
            dataset_payload.get("metadata"),
            field_name="bundle manifest.dataset.metadata",
        ),
        dataset_tables=require_mapping(
            dataset_payload.get("tables"),
            field_name="bundle manifest.dataset.tables",
        ),
        references_metadata=require_mapping(
            references_payload.get("metadata"),
            field_name="bundle manifest.resolved_references.metadata",
        ),
        reference_tables=require_mapping(
            references_payload.get("tables"),
            field_name="bundle manifest.resolved_references.tables",
        ),
        scoring_tables=require_mapping(
            scoring_payload.get("tables"),
            field_name="bundle manifest.upstream_kinase_outputs.scoring.tables",
        ),
        prediction_tables=require_mapping(
            prediction_payload.get("tables"),
            field_name="bundle manifest.upstream_kinase_outputs.prediction.tables",
        ),
        activity_tables=require_mapping(
            activity_payload.get("tables"),
            field_name="bundle manifest.upstream_kinase_outputs.activity.tables",
        ),
        signalome_metadata=require_mapping(
            signalome_outputs_payload.get("metadata"),
            field_name="bundle manifest.signalome_outputs.metadata",
        ),
        signalome_tables=require_mapping(
            signalome_outputs_payload.get("tables"),
            field_name="bundle manifest.signalome_outputs.tables",
        ),
        config_snapshot_path=require_str(
            payload.get("config_snapshot"),
            field_name="bundle manifest.config_snapshot",
        ),
    )
