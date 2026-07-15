"""Shared helpers for workflow result caveat assembly."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import TYPE_CHECKING

import pandas as pd

from phospy.contracts.configs.localisation import LocalisationRequirement
from phospy.contracts.result_caveats import ResultCaveat, ResultCaveatSeverity
from phospy.provenance.models import TrustedDatasetConstructionAssertions
from phospy.science.datasets.direct_construction import (
    DIRECT_CONSTRUCTION_SOURCE,
    DIRECT_CONSTRUCTION_WORKFLOW_NAME,
)
from phospy.validation.identity_contracts import ReferenceContextCompatibilityWarning

if TYPE_CHECKING:
    from phospy.science.datasets.models import AnalysisReadyPhosphoDataset

_TRUSTED_CONSTRUCTION_ASSERTION_FIELDS = (
    "sequence_user_asserted",
    "identity_user_asserted",
    "quantitative_meaning_user_asserted",
    "reference_context_user_asserted",
)


def build_localisation_policy_details(
    *,
    site_metadata: pd.DataFrame,
    requirement: LocalisationRequirement,
    workflow_scope: str,
) -> dict[str, object]:
    """Return compact localisation policy diagnostics for result caveats."""

    resolved_column = _resolve_localisation_column_name(site_metadata)
    site_count = int(site_metadata.shape[0])
    known_count = 0
    missing_count = site_count
    if resolved_column is not None:
        values = site_metadata.loc[:, resolved_column]
        missing_mask = values.isna() | values.map(
            lambda value: isinstance(value, str) and value.strip() == ""
        )
        missing_count = int(missing_mask.sum())
        known_count = int(site_count - missing_count)
    return {
        "workflow_scope": workflow_scope,
        "policy": str(requirement.policy),
        "require_present": bool(requirement.require_present),
        "minimum_probability": (
            None
            if requirement.minimum_probability is None
            else float(requirement.minimum_probability)
        ),
        "resolved_column": resolved_column,
        "site_count": site_count,
        "known_localisation_count": known_count,
        "missing_or_unknown_localisation_count": missing_count,
    }


def is_permissive_localisation_requirement(
    requirement: LocalisationRequirement,
) -> bool:
    """Return whether localisation policy permits unknown or low-confidence sites."""

    return requirement.minimum_probability is None


def build_reference_context_compatibility_caveat(
    warning: ReferenceContextCompatibilityWarning,
    *,
    policy: object,
    workflow_scope: str,
) -> ResultCaveat:
    """Convert a reference-context compatibility warning into a result caveat."""

    details = warning.to_payload()
    details["policy"] = str(policy)
    details["workflow_scope"] = workflow_scope
    return ResultCaveat(
        code=warning.code,
        severity="warning",
        message=warning.message,
        details=details,
    )


def build_direct_trusted_dataset_construction_caveat(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    code: str,
    workflow_scope: str,
    workflow_label: str,
) -> ResultCaveat | None:
    """Build a workflow caveat for direct trusted dataset construction."""

    details = direct_trusted_dataset_construction_details(
        dataset=dataset,
        workflow_scope=workflow_scope,
    )
    if details is None:
        return None
    severity = direct_trusted_dataset_caveat_severity(details)
    message = (
        "Input dataset was directly constructed as trusted analysis-ready state; "
        f"{workflow_label} did not rerun dataset-building validation."
    )
    if severity == "warning":
        message = (
            message
            + " Trusted construction assertion metadata is absent or incomplete."
        )
    return ResultCaveat(
        code=code,
        severity=severity,
        message=message,
        details=details,
    )


def direct_trusted_dataset_construction_details(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    workflow_scope: str,
) -> dict[str, object] | None:
    """Return direct trusted-construction details for result caveats."""

    provenance = dataset.provenance
    construction = (
        {}
        if provenance is None
        else _construction_payload(provenance.workflow_parameters)
    )
    assertions = dataset.trusted_construction_assertions
    source = construction.get("source")
    if (
        assertions is None
        and provenance is not None
        and provenance.workflow_name != DIRECT_CONSTRUCTION_WORKFLOW_NAME
        and source != DIRECT_CONSTRUCTION_SOURCE
    ):
        return None

    assertion_payload = _trusted_assertion_payload(
        assertions=assertions,
        construction=construction,
    )
    missing_assertions = _missing_trusted_assertions(assertion_payload)
    metadata_provided = assertion_payload.get("assertion_metadata_provided")
    metadata_provided_bool = (
        bool(metadata_provided)
        if isinstance(
            metadata_provided,
            bool,
        )
        else False
    )
    details: dict[str, object] = {
        "workflow_scope": workflow_scope,
        "dataset_workflow_name": None
        if provenance is None
        else provenance.workflow_name,
        "construction_source": (
            DIRECT_CONSTRUCTION_SOURCE if source is None else str(source)
        ),
        "trusted_assertion_metadata_provided": metadata_provided_bool,
        "trusted_construction_assertions": assertion_payload,
        "missing_trusted_assertions": list(missing_assertions),
    }
    builder_used = construction.get("builder_used")
    if isinstance(builder_used, bool):
        details["builder_used"] = builder_used
    warning = construction.get("warning")
    if isinstance(warning, str) and warning.strip():
        details["construction_warning"] = warning.strip()
    assertion_warning = construction.get("assertion_warning")
    if isinstance(assertion_warning, str) and assertion_warning.strip():
        details["assertion_warning"] = assertion_warning.strip()
    return details


def direct_trusted_dataset_caveat_severity(
    details: Mapping[str, object],
) -> ResultCaveatSeverity:
    """Return warning when trusted construction assertion metadata is incomplete."""

    missing = details.get("missing_trusted_assertions")
    metadata_provided = details.get("trusted_assertion_metadata_provided")
    if metadata_provided is not True:
        return "warning"
    if isinstance(missing, list) and missing:
        return "warning"
    return "info"


def deduplicate_caveats(
    caveats: Iterable[ResultCaveat],
) -> tuple[ResultCaveat, ...]:
    """Return caveats with duplicate code/details pairs removed, preserving order."""

    observed: set[tuple[str, tuple[tuple[str, str], ...]]] = set()
    unique: list[ResultCaveat] = []
    for caveat in caveats:
        detail_key = tuple(
            sorted((str(key), repr(value)) for key, value in caveat.details.items())
        )
        key = (caveat.code, detail_key)
        if key in observed:
            continue
        observed.add(key)
        unique.append(caveat)
    return tuple(unique)


def _resolve_localisation_column_name(site_metadata: pd.DataFrame) -> str | None:
    if "localisation_confidence" in site_metadata.columns:
        return "localisation_confidence"
    if "localisation_probability" in site_metadata.columns:
        return "localisation_probability"
    return None


def _construction_payload(parameters: Mapping[str, object]) -> Mapping[str, object]:
    construction = parameters.get("construction")
    if isinstance(construction, Mapping):
        return construction
    return {}


def _trusted_assertion_payload(
    *,
    assertions: TrustedDatasetConstructionAssertions | None,
    construction: Mapping[str, object],
) -> dict[str, object]:
    if assertions is not None:
        return assertions.to_payload()
    raw_payload = construction.get("trusted_construction_assertions")
    if isinstance(raw_payload, Mapping):
        return {str(key): value for key, value in raw_payload.items()}
    return TrustedDatasetConstructionAssertions.missing().to_payload()


def _missing_trusted_assertions(
    assertion_payload: Mapping[str, object],
) -> tuple[str, ...]:
    raw_missing = assertion_payload.get("missing_assertions")
    if isinstance(raw_missing, list):
        return tuple(str(item) for item in raw_missing)
    return tuple(
        field_name
        for field_name in _TRUSTED_CONSTRUCTION_ASSERTION_FIELDS
        if assertion_payload.get(field_name) is not True
    )


__all__ = [
    "build_direct_trusted_dataset_construction_caveat",
    "build_reference_context_compatibility_caveat",
    "build_localisation_policy_details",
    "deduplicate_caveats",
    "direct_trusted_dataset_caveat_severity",
    "direct_trusted_dataset_construction_details",
    "is_permissive_localisation_requirement",
]
