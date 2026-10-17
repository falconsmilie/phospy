"""Localisation-evidence details used by workflow caveat assembly."""

from __future__ import annotations

import pandas as pd

from phospy.evidence.config import LocalisationRequirement


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
    """Return whether policy permits unknown or low-confidence sites."""

    return requirement.minimum_probability is None


def _resolve_localisation_column_name(site_metadata: pd.DataFrame) -> str | None:
    if "localisation_confidence" in site_metadata.columns:
        return "localisation_confidence"
    if "localisation_probability" in site_metadata.columns:
        return "localisation_probability"
    return None


__all__ = [
    "build_localisation_policy_details",
    "is_permissive_localisation_requirement",
]
