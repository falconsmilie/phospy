"""Internal kinase-substrate contribution collection."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import SupportsFloat, cast

import pandas as pd

from phospy.science.references.models import ReferenceBundle
from phospy.science.tables.kinase import (
    KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS,
    KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_BELOW_MIN_SUBSTRATES,
    KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_MISSING_SCORE_VALUE,
    KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NO_SCORE_COLUMN,
    KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_IN_PROFILE_SUPPORT,
    KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_QUANTIFIED,
    KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED,
    KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED,
)


@dataclass(frozen=True, slots=True)
class KinaseSubstrateContributionReferenceSource:
    """Reference metadata already carried by the resolved bundle."""

    source_name: str | None = None
    source_version: str | None = None
    bundle_id: str | None = None
    identifier_namespace: str | None = None


def reference_source_from_bundle(
    references: ReferenceBundle,
) -> KinaseSubstrateContributionReferenceSource:
    """Extract nullable reference source fields for contribution records."""

    manifest = references.manifest
    provenance = references.provenance
    return KinaseSubstrateContributionReferenceSource(
        source_name=(
            manifest.source_name
            if manifest is not None
            else None
            if provenance is None
            else provenance.source_name
        ),
        source_version=(
            manifest.source_version
            if manifest is not None
            else None
            if provenance is None
            else provenance.source_version
        ),
        bundle_id=(
            manifest.bundle_id
            if manifest is not None
            else None
            if provenance is None
            else provenance.bundle_id
        ),
        identifier_namespace=(
            manifest.identifier_namespace
            if manifest is not None
            else None
            if provenance is None
            else provenance.identifier_namespace
        ),
    )


def build_kinase_substrate_contribution_table(
    *,
    kinase_substrate_map: pd.DataFrame,
    scoring_values: pd.DataFrame,
    score_component: str,
    quantified_substrates: Mapping[str, Sequence[str]],
    substrate_counts: pd.Series,
    min_substrates: int,
    score_source_matrix: pd.DataFrame | None = None,
    reference_source: KinaseSubstrateContributionReferenceSource | None = None,
    display_reference_matching: Mapping[str, object] | None = None,
) -> pd.DataFrame:
    """Build internal substrate-level evidence rows for kinase scoring."""

    source = reference_source or KinaseSubstrateContributionReferenceSource()
    contribution_map = _deduplicate_contribution_map(kinase_substrate_map)
    ambiguous_site_keys, ambiguous_display_ids = _collect_ambiguous_reference_ids(
        kinase_substrate_map=contribution_map,
        display_reference_matching=display_reference_matching,
    )
    rows: list[dict[str, object]] = []
    scoring_index = {str(value) for value in scoring_values.index.tolist()}
    scoring_columns = {str(value) for value in scoring_values.columns.tolist()}
    support_sites_by_kinase = {
        str(kinase): {str(site) for site in sites}
        for kinase, sites in quantified_substrates.items()
    }
    fallback_counts = _fallback_substrate_counts(
        kinase_substrate_map=contribution_map,
        scoring_index=scoring_index,
    )

    for record in contribution_map.to_dict(orient="records"):
        kinase = str(record["kinase"])
        substrate_site = str(record["substrate_site"])
        substrate_identifier = _optional_text(record.get("display_id"))
        score_value = _lookup_float(scoring_values, row=substrate_site, column=kinase)
        score_source = _lookup_text(
            score_source_matrix,
            row=substrate_site,
            column=kinase,
        )
        if score_source is None and kinase in scoring_columns:
            score_source = str(score_component)
        exclusion_reason = _resolve_exclusion_reason(
            kinase=kinase,
            substrate_site=substrate_site,
            score_value=score_value,
            scoring_index=scoring_index,
            scoring_columns=scoring_columns,
            support_sites_by_kinase=support_sites_by_kinase,
            substrate_counts=substrate_counts,
            fallback_counts=fallback_counts,
            min_substrates=int(min_substrates),
        )
        status = (
            KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED
            if exclusion_reason is None
            else KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED
        )
        rows.append(
            {
                "kinase": kinase,
                "substrate_site": substrate_site,
                "substrate_identifier": substrate_identifier,
                "value_used_in_scoring": score_value,
                "score_component": str(score_component),
                "score_source": score_source,
                "reference_source_name": source.source_name,
                "reference_source_version": source.source_version,
                "reference_bundle_id": source.bundle_id,
                "reference_identifier_namespace": source.identifier_namespace,
                "status": status,
                "exclusion_reason": exclusion_reason,
                "ambiguous": bool(
                    substrate_site in ambiguous_site_keys
                    or (
                        substrate_identifier is not None
                        and substrate_identifier in ambiguous_display_ids
                    )
                ),
            }
        )
    return pd.DataFrame.from_records(
        rows,
        columns=pd.Index(KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS),
    )


def _deduplicate_contribution_map(kinase_substrate_map: pd.DataFrame) -> pd.DataFrame:
    columns = ["kinase", "substrate_site"]
    if "display_id" in kinase_substrate_map.columns:
        columns.append("display_id")
    return (
        kinase_substrate_map.loc[:, columns]
        .copy(deep=True)
        .drop_duplicates(ignore_index=True)
    )


def _fallback_substrate_counts(
    *,
    kinase_substrate_map: pd.DataFrame,
    scoring_index: set[str],
) -> dict[str, int]:
    counts: dict[str, int] = {}
    for kinase, grouped in kinase_substrate_map.groupby("kinase", sort=False):
        unique_sites = dict.fromkeys(
            str(site) for site in grouped.loc[:, "substrate_site"].tolist()
        )
        counts[str(kinase)] = sum(1 for site in unique_sites if site in scoring_index)
    return counts


def _resolve_exclusion_reason(
    *,
    kinase: str,
    substrate_site: str,
    score_value: float | None,
    scoring_index: set[str],
    scoring_columns: set[str],
    support_sites_by_kinase: Mapping[str, set[str]],
    substrate_counts: pd.Series,
    fallback_counts: Mapping[str, int],
    min_substrates: int,
) -> str | None:
    if substrate_site not in scoring_index:
        return KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_QUANTIFIED
    substrate_count = _substrate_count(
        substrate_counts=substrate_counts,
        fallback_counts=fallback_counts,
        kinase=kinase,
    )
    if substrate_count < min_substrates:
        return KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_BELOW_MIN_SUBSTRATES
    if kinase not in scoring_columns:
        return KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NO_SCORE_COLUMN
    if substrate_site not in support_sites_by_kinase.get(kinase, set()):
        return KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_IN_PROFILE_SUPPORT
    if score_value is None or pd.isna(score_value):
        return KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_MISSING_SCORE_VALUE
    return None


def _substrate_count(
    *,
    substrate_counts: pd.Series,
    fallback_counts: Mapping[str, int],
    kinase: str,
) -> int:
    if kinase in substrate_counts.index:
        return int(substrate_counts.loc[kinase])
    return int(fallback_counts.get(kinase, 0))


def _lookup_float(
    frame: pd.DataFrame,
    *,
    row: str,
    column: str,
) -> float | None:
    if row not in frame.index or column not in frame.columns:
        return None
    value = frame.at[row, column]
    if pd.isna(value):
        return float("nan")
    return float(cast(SupportsFloat, value))


def _lookup_text(
    frame: pd.DataFrame | None,
    *,
    row: str,
    column: str,
) -> str | None:
    if frame is None or row not in frame.index or column not in frame.columns:
        return None
    value = frame.at[row, column]
    if pd.isna(value):
        return None
    return str(value)


def _optional_text(value: object | None) -> str | None:
    if value is None or value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    text = str(value).strip()
    return text or None


def _collect_ambiguous_reference_ids(
    *,
    kinase_substrate_map: pd.DataFrame,
    display_reference_matching: Mapping[str, object] | None,
) -> tuple[set[str], set[str]]:
    site_keys: set[str] = set()
    display_ids: set[str] = set()
    if isinstance(display_reference_matching, Mapping):
        matches = display_reference_matching.get(
            "one_to_many_display_reference_matches",
        )
        if isinstance(matches, list):
            for match in matches:
                if not isinstance(match, Mapping):
                    continue
                display_id = _optional_text(match.get("display_id"))
                if display_id is not None:
                    display_ids.add(display_id)
                matched_site_keys = match.get("site_keys")
                if isinstance(matched_site_keys, list):
                    site_keys.update(
                        str(site_key)
                        for site_key in matched_site_keys
                        if _optional_text(site_key) is not None
                    )
    if "display_id" in kinase_substrate_map.columns:
        for display_id, grouped in kinase_substrate_map.groupby(
            "display_id",
            sort=False,
            dropna=True,
        ):
            sites = {str(site) for site in grouped.loc[:, "substrate_site"].tolist()}
            if len(sites) < 2:
                continue
            text = _optional_text(display_id)
            if text is not None:
                display_ids.add(text)
            site_keys.update(sites)
    return site_keys, display_ids


__all__ = [
    "KINASE_SUBSTRATE_CONTRIBUTION_COLUMNS",
    "KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_BELOW_MIN_SUBSTRATES",
    "KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_MISSING_SCORE_VALUE",
    "KINASE_SUBSTRATE_CONTRIBUTION_EXCLUDED_NOT_QUANTIFIED",
    "KINASE_SUBSTRATE_CONTRIBUTION_STATUS_EXCLUDED",
    "KINASE_SUBSTRATE_CONTRIBUTION_STATUS_INCLUDED",
    "KinaseSubstrateContributionReferenceSource",
    "build_kinase_substrate_contribution_table",
    "reference_source_from_bundle",
]
