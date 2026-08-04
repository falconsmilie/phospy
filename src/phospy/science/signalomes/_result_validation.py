"""Private signalome result identity validation."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import pandas as pd

from phospy.errors.workflows import WorkflowStageError
from phospy.science.datasets.internal_view import DatasetInternalView
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.signalomes.constants import (
    DISPLAY_ID_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_COLUMN,
    EXPANDED_SIGNALOME_ROW_KIND_SITE,
    SITE_KEY_COLUMN,
)
from phospy.science.sites.validation import require_site_key_series


@dataclass(frozen=True, slots=True)
class _SignalomeDatasetIdentityLookup:
    """Typed site-key to display-id lookup for Signalome result validation."""

    display_id_by_site_key: Mapping[str, str] = field(repr=False)

    def missing_site_keys(self, site_keys: list[str]) -> list[str]:
        return [
            site_key
            for site_key in dict.fromkeys(site_keys)
            if site_key not in self.display_id_by_site_key
        ]

    def expected_display_id(self, site_key: str) -> str:
        return self.display_id_by_site_key[site_key]


def validate_signalome_result_site_level_identity(
    *,
    module_assignments: pd.DataFrame,
    expanded_signalome: pd.DataFrame | None,
    site_membership: pd.DataFrame | None,
    dataset: AnalysisReadyPhosphoDataset | None = None,
    site_metadata: pd.DataFrame | None = None,
    dataset_identity: _SignalomeDatasetIdentityLookup | None = None,
) -> None:
    """Validate signalome result site identity against dataset metadata."""

    dataset_identity = _resolve_signalome_dataset_identity(
        dataset=dataset,
        site_metadata=site_metadata,
        dataset_identity=dataset_identity,
    )
    _validate_site_level_signalome_rows(
        table=module_assignments,
        field_name="signalome_result.module_assignments.table",
        dataset_identity=dataset_identity,
    )
    _validate_expanded_signalome_identity(
        expanded_signalome,
        dataset_identity=dataset_identity,
    )
    _validate_site_level_signalome_rows(
        table=site_membership,
        field_name="signalome_result.site_membership",
        dataset_identity=dataset_identity,
    )


def _validate_expanded_signalome_identity(
    expanded_signalome: pd.DataFrame | None,
    *,
    dataset_identity: _SignalomeDatasetIdentityLookup,
) -> None:
    missing = [
        column
        for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN)
        if expanded_signalome is not None and column not in expanded_signalome.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowStageError(
            f"signalome_result.expanded_signalome is missing required columns: {joined}"
        )
    if expanded_signalome is None or expanded_signalome.empty:
        return
    site_rows = expanded_signalome
    if EXPANDED_SIGNALOME_ROW_KIND_COLUMN in expanded_signalome.columns:
        site_rows = expanded_signalome.loc[
            expanded_signalome.loc[:, EXPANDED_SIGNALOME_ROW_KIND_COLUMN].astype(str)
            == EXPANDED_SIGNALOME_ROW_KIND_SITE,
            :,
        ]
    _validate_site_level_signalome_rows(
        table=site_rows,
        field_name="signalome_result.expanded_signalome",
        dataset_identity=dataset_identity,
    )


def _validate_site_level_signalome_rows(
    *,
    table: pd.DataFrame | None,
    field_name: str,
    dataset_identity: _SignalomeDatasetIdentityLookup,
) -> None:
    if table is None:
        return
    missing = [
        column for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN) if column not in table
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowStageError(f"{field_name} is missing required columns: {joined}")
    if table.empty:
        return
    for column_name in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN):
        invalid_count = int(
            sum(
                1
                for value in table.loc[:, column_name].tolist()
                if not isinstance(value, str) or value.strip() == ""
            )
        )
        if invalid_count:
            raise WorkflowStageError(
                f"{field_name} site rows require non-empty "
                f"{column_name} values; invalid_count={invalid_count}"
            )
    site_keys = table.loc[:, SITE_KEY_COLUMN].astype(str)
    require_site_key_series(
        site_keys,
        field_name=f"{field_name}.{SITE_KEY_COLUMN}",
        error_type=WorkflowStageError,
    )
    display_ids = table.loc[:, DISPLAY_ID_COLUMN].astype(str)
    missing_site_keys = dataset_identity.missing_site_keys(site_keys.tolist())
    if missing_site_keys:
        preview = ", ".join(repr(value) for value in missing_site_keys[:5])
        suffix = "" if len(missing_site_keys) <= 5 else " ..."
        raise WorkflowStageError(
            f"{field_name}.{SITE_KEY_COLUMN} values must align to "
            "signalome_result.dataset; missing_site_keys="
            f"{preview}{suffix}"
        )
    mismatches = [
        f"{site_key!r}: observed={display_id!r}, "
        f"expected={dataset_identity.expected_display_id(site_key)!r}"
        for site_key, display_id in zip(
            site_keys.tolist(),
            display_ids.tolist(),
            strict=True,
        )
        if dataset_identity.expected_display_id(site_key) != display_id
    ]
    if mismatches:
        preview = "; ".join(mismatches[:5])
        suffix = "" if len(mismatches) <= 5 else " ..."
        raise WorkflowStageError(
            f"{field_name}.{DISPLAY_ID_COLUMN} values must match "
            "signalome_result.dataset.site_metadata.display_id for each site_key; "
            f"mismatches={preview}{suffix}"
        )


def _resolve_signalome_dataset_identity(
    *,
    dataset: AnalysisReadyPhosphoDataset | None,
    site_metadata: pd.DataFrame | None,
    dataset_identity: _SignalomeDatasetIdentityLookup | None,
) -> _SignalomeDatasetIdentityLookup:
    provided_sources = (
        int(dataset is not None)
        + int(site_metadata is not None)
        + int(dataset_identity is not None)
    )
    if provided_sources != 1:
        raise WorkflowStageError(
            "signalome result identity validation requires exactly one dataset "
            "identity source: dataset, site_metadata, or dataset_identity"
        )
    if dataset_identity is not None:
        return dataset_identity
    if site_metadata is None:
        if dataset is None:
            raise WorkflowStageError(
                "signalome result identity validation requires dataset metadata"
            )
        site_metadata = DatasetInternalView(dataset).site_metadata
    return _signalome_dataset_identity_lookup_from_site_metadata(site_metadata)


def _signalome_dataset_identity_lookup_from_site_metadata(
    site_metadata: pd.DataFrame,
) -> _SignalomeDatasetIdentityLookup:
    missing = [
        column
        for column in (SITE_KEY_COLUMN, DISPLAY_ID_COLUMN)
        if column not in site_metadata.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise WorkflowStageError(
            "signalome_result.dataset.site_metadata is missing required columns: "
            f"{joined}"
        )
    site_keys = site_metadata.loc[:, SITE_KEY_COLUMN].astype(str)
    require_site_key_series(
        site_keys,
        field_name=f"signalome_result.dataset.site_metadata.{SITE_KEY_COLUMN}",
        error_type=WorkflowStageError,
    )
    display_ids = site_metadata.loc[:, DISPLAY_ID_COLUMN]
    invalid_display_id_count = int(
        sum(
            1
            for value in display_ids.tolist()
            if not isinstance(value, str) or value.strip() == ""
        )
    )
    if invalid_display_id_count:
        raise WorkflowStageError(
            "signalome_result.dataset.site_metadata.display_id values must be "
            "non-empty strings; "
            f"invalid_count={invalid_display_id_count}"
        )
    return _SignalomeDatasetIdentityLookup(
        display_id_by_site_key=MappingProxyType(
            {
                site_key: display_id
                for site_key, display_id in zip(
                    site_keys.tolist(),
                    display_ids.astype(str).tolist(),
                    strict=True,
                )
            }
        )
    )
