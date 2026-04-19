"""Current-contract manifest serialization and parsing for kinase bundles."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from phospy.api.results import KinaseWorkflowResult
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._kinase.constants import (
    CONFIG_SNAPSHOT_RELATIVE_PATH,
    KINASE_BUNDLE_KIND,
    KINASE_BUNDLE_MANIFEST_VERSION,
)
from phospy.io.bundles._shared.primitives import (
    require_int,
    require_mapping,
    require_str,
)
from phospy.io.bundles._shared.transformation_state import (
    transformation_state_to_payload,
)


@dataclass(frozen=True, slots=True)
class KinaseManifestSections:
    """Decoded v1 manifest sections needed to load kinase bundles."""

    manifest_version: int
    dataset_metadata: Mapping[str, object]
    dataset_tables: Mapping[str, object]
    references_metadata: Mapping[str, object]
    reference_tables: Mapping[str, object]
    scoring_tables: Mapping[str, object]
    prediction_tables: Mapping[str, object]
    activity_tables: Mapping[str, object]
    config_snapshot_path: str


def build_manifest(
    *,
    result: KinaseWorkflowResult,
    table_format: str,
    dataset_tables: Mapping[str, object],
    reference_tables: Mapping[str, object],
    scoring_tables: Mapping[str, object],
    prediction_tables: Mapping[str, object],
    activity_tables: Mapping[str, object],
) -> dict[str, object]:
    """Build the v1 manifest payload from current bundle contract data."""

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
                "transformation_state": transformation_state_to_payload(
                    result.dataset.transformation_state
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
                "tables": dict(activity_tables),
            },
        },
        "config_snapshot": CONFIG_SNAPSHOT_RELATIVE_PATH,
    }


def parse_manifest(payload: Mapping[str, object]) -> KinaseManifestSections:
    """Parse and validate current-contract kinase manifest payload."""

    bundle_type = require_str(
        payload.get("bundle_type"),
        field_name="bundle manifest.bundle_type",
    )
    if bundle_type != KINASE_BUNDLE_KIND:
        raise PhosPyInputError(
            "unsupported bundle manifest bundle_type "
            f"'{bundle_type}'; expected '{KINASE_BUNDLE_KIND}'"
        )
    manifest_version = require_int(
        payload.get("manifest_version"),
        field_name="bundle manifest.manifest_version",
    )
    if manifest_version != KINASE_BUNDLE_MANIFEST_VERSION:
        raise PhosPyInputError(
            "unsupported bundle manifest version "
            f"'{manifest_version}'; expected {KINASE_BUNDLE_MANIFEST_VERSION}"
        )

    dataset_payload = require_mapping(
        payload.get("dataset"),
        field_name="bundle manifest.dataset",
    )
    references_payload = require_mapping(
        payload.get("resolved_references"),
        field_name="bundle manifest.resolved_references",
    )
    outputs_payload = require_mapping(
        payload.get("outputs"),
        field_name="bundle manifest.outputs",
    )
    scoring_payload = require_mapping(
        outputs_payload.get("scoring"),
        field_name="bundle manifest.outputs.scoring",
    )
    prediction_payload = require_mapping(
        outputs_payload.get("prediction"),
        field_name="bundle manifest.outputs.prediction",
    )
    activity_payload = require_mapping(
        outputs_payload.get("activity"),
        field_name="bundle manifest.outputs.activity",
    )

    return KinaseManifestSections(
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
            field_name="bundle manifest.outputs.scoring.tables",
        ),
        prediction_tables=require_mapping(
            prediction_payload.get("tables"),
            field_name="bundle manifest.outputs.prediction.tables",
        ),
        activity_tables=require_mapping(
            activity_payload.get("tables"),
            field_name="bundle manifest.outputs.activity.tables",
        ),
        config_snapshot_path=require_str(
            payload.get("config_snapshot"),
            field_name="bundle manifest.config_snapshot",
        ),
    )
