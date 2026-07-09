"""Shared helpers for workflow result caveat assembly."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from phospy.contracts.configs.localisation import LocalisationRequirement
from phospy.contracts.result_caveats import ResultCaveat
from phospy.validation.identity_contracts import ReferenceContextCompatibilityWarning


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
    workflow_scope: str,
) -> ResultCaveat:
    """Convert a reference-context compatibility warning into a result caveat."""

    details = warning.to_payload()
    details["workflow_scope"] = workflow_scope
    return ResultCaveat(
        code=warning.code,
        severity="warning",
        message=warning.message,
        details=details,
    )


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


__all__ = [
    "build_reference_context_compatibility_caveat",
    "build_localisation_policy_details",
    "deduplicate_caveats",
    "is_permissive_localisation_requirement",
]
