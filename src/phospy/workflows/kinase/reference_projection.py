"""Kinase-substrate reference projection onto dataset site_key rows."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import NoReturn

import pandas as pd

from phospy.contracts.configs.kinase import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES,
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS,
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR,
    KinaseReferenceDisplayAmbiguityPolicy,
)
from phospy.contracts.kinase_reference_projection import (
    KINASE_REFERENCE_PROJECTION_INTERPRETER_VERSION as _INTERPRETER_VERSION,
)
from phospy.contracts.kinase_reference_projection import (
    KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_DISPLAY_ID as _SOURCE_IDENTIFIER_KIND_DATASET_DISPLAY_ID,
)
from phospy.contracts.kinase_reference_projection import (
    KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_DATASET_SITE_KEY as _SOURCE_IDENTIFIER_KIND_DATASET_SITE_KEY,
)
from phospy.contracts.kinase_reference_projection import (
    KINASE_REFERENCE_PROJECTION_SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE as _SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE,
)
from phospy.contracts.kinase_reference_projection import (
    KinaseReferenceProjectionSummary,
)
from phospy.errors.workflows import WorkflowBoundaryError

_PROJECTED_COLUMNS = ["kinase", "substrate_site", "display_id"]


@dataclass(slots=True)
class _AmbiguousDisplayReferenceRows:
    site_keys: tuple[str, ...]
    reference_rows: list[dict[str, object]]


@dataclass(frozen=True, slots=True)
class KinaseReferenceDisplayAmbiguityDiagnostic:
    """One ambiguous display-level kinase-substrate reference projection."""

    display_id: str
    site_keys: tuple[str, ...]
    reference_rows: tuple[Mapping[str, object], ...]
    projected_rows: int
    interpreter_version: str = _INTERPRETER_VERSION

    @property
    def matched_row_count(self) -> int:
        return len(self.site_keys)

    @property
    def reference_row_count(self) -> int:
        return len(self.reference_rows)

    def to_payload(self) -> dict[str, object]:
        return {
            "display_id": self.display_id,
            "site_keys": list(self.site_keys),
            "matched_row_count": int(self.matched_row_count),
            "reference_row_count": int(self.reference_row_count),
            "reference_rows": [dict(row) for row in self.reference_rows],
            "projected_rows": int(self.projected_rows),
            "interpreter_version": self.interpreter_version,
        }


@dataclass(frozen=True, slots=True)
class KinaseReferenceProjectionResult:
    """Projected kinase-substrate map plus projection diagnostics."""

    kinase_substrate_map: pd.DataFrame
    ambiguity_policy: KinaseReferenceDisplayAmbiguityPolicy
    ambiguity_diagnostics: tuple[KinaseReferenceDisplayAmbiguityDiagnostic, ...]
    projection_summary: KinaseReferenceProjectionSummary

    def display_reference_matching_payload(self) -> dict[str, object]:
        matches = [item.to_payload() for item in self.ambiguity_diagnostics]
        return {
            "reference_key": "display_id",
            "dataset_row_identity": "site_key",
            "ambiguity_policy": self.ambiguity_policy,
            "one_to_many_display_reference_match_count": len(matches),
            "one_to_many_display_reference_site_key_rows": sum(
                item.matched_row_count for item in self.ambiguity_diagnostics
            ),
            "one_to_many_display_reference_matches": matches,
            "interpreter_version": _INTERPRETER_VERSION,
        }


class KinaseReferenceProjector:
    """Project kinase-substrate reference rows to dataset site_key identity."""

    def run(
        self,
        *,
        reference_kinase_substrate_map: pd.DataFrame,
        site_identity_map: pd.DataFrame,
        ambiguity_policy: KinaseReferenceDisplayAmbiguityPolicy,
    ) -> KinaseReferenceProjectionResult:
        self._validate_ambiguity_policy(ambiguity_policy)
        display_lookup, site_key_to_display_id = self._build_identity_lookups(
            site_identity_map=site_identity_map
        )
        ambiguity_rows = self._collect_display_ambiguities(
            reference_kinase_substrate_map=reference_kinase_substrate_map,
            display_lookup=display_lookup,
        )
        preliminary_diagnostics = self._build_ambiguity_diagnostics(
            ambiguity_rows=ambiguity_rows,
            projected_map=None,
        )
        if (
            preliminary_diagnostics
            and ambiguity_policy == KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR
        ):
            self._raise_display_ambiguity_error(preliminary_diagnostics)
        projected = self._project_rows(
            reference_kinase_substrate_map=reference_kinase_substrate_map,
            display_lookup=display_lookup,
            site_key_to_display_id=site_key_to_display_id,
        )
        diagnostics = self._build_ambiguity_diagnostics(
            ambiguity_rows=ambiguity_rows,
            projected_map=projected,
        )
        projection_summary = self._build_projection_summary(
            reference_kinase_substrate_map=reference_kinase_substrate_map,
            display_lookup=display_lookup,
            site_key_to_display_id=site_key_to_display_id,
            projected_map=projected,
            ambiguity_diagnostics=diagnostics,
        )
        return KinaseReferenceProjectionResult(
            kinase_substrate_map=projected,
            ambiguity_policy=ambiguity_policy,
            ambiguity_diagnostics=diagnostics,
            projection_summary=projection_summary,
        )

    @staticmethod
    def _validate_ambiguity_policy(
        ambiguity_policy: KinaseReferenceDisplayAmbiguityPolicy,
    ) -> None:
        if ambiguity_policy in KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES:
            return
        supported = ", ".join(sorted(KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES))
        raise WorkflowBoundaryError(
            seam="kinase.interpreter.reference_display_ambiguity_policy",
            next_action=(
                "use reference_display_ambiguity_policy='error' or "
                "'allow_with_diagnostics' on KinaseWorkflowRequest"
            ),
            details={
                "reference_display_ambiguity_policy": str(ambiguity_policy),
                "supported_policies": supported,
            },
            message_prefix="kinase workflow boundary validation failed",
        )

    @staticmethod
    def _build_identity_lookups(
        *,
        site_identity_map: pd.DataFrame,
    ) -> tuple[dict[str, list[str]], dict[str, str]]:
        display_lookup: dict[str, list[str]] = {}
        site_key_to_display_id: dict[str, str] = {}
        for site_key, display_id in site_identity_map.loc[
            :, ["site_key", "display_id"]
        ].itertuples(index=False):
            site_key_value = str(site_key)
            display_value = str(display_id)
            display_lookup.setdefault(display_value, []).append(site_key_value)
            site_key_to_display_id[site_key_value] = display_value
        return display_lookup, site_key_to_display_id

    @staticmethod
    def _collect_display_ambiguities(
        *,
        reference_kinase_substrate_map: pd.DataFrame,
        display_lookup: dict[str, list[str]],
    ) -> dict[str, _AmbiguousDisplayReferenceRows]:
        ambiguity_rows: dict[str, _AmbiguousDisplayReferenceRows] = {}
        for row_position, row in enumerate(
            reference_kinase_substrate_map.loc[
                :, ["kinase", "substrate_site"]
            ].itertuples(index=True, name=None)
        ):
            row_index, kinase, substrate_site = row
            display_id = str(substrate_site)
            site_keys = tuple(display_lookup.get(display_id, ()))
            if len(site_keys) < 2:
                continue
            entry = ambiguity_rows.setdefault(
                display_id,
                _AmbiguousDisplayReferenceRows(
                    site_keys=site_keys,
                    reference_rows=[],
                ),
            )
            entry.reference_rows.append(
                {
                    "row_position": int(row_position),
                    "row_index": str(row_index),
                    "kinase": str(kinase),
                    "substrate_site": display_id,
                }
            )
        return ambiguity_rows

    @staticmethod
    def _project_rows(
        *,
        reference_kinase_substrate_map: pd.DataFrame,
        display_lookup: dict[str, list[str]],
        site_key_to_display_id: dict[str, str],
    ) -> pd.DataFrame:
        rows: list[dict[str, str]] = []
        for kinase, substrate_site in reference_kinase_substrate_map.loc[
            :, ["kinase", "substrate_site"]
        ].itertuples(index=False):
            substrate_value = str(substrate_site)
            if substrate_value in site_key_to_display_id:
                rows.append(
                    {
                        "kinase": str(kinase),
                        "substrate_site": substrate_value,
                        "display_id": site_key_to_display_id[substrate_value],
                    }
                )
                continue
            matched_site_keys = display_lookup.get(substrate_value, [])
            for site_key in matched_site_keys:
                rows.append(
                    {
                        "kinase": str(kinase),
                        "substrate_site": str(site_key),
                        "display_id": substrate_value,
                    }
                )
        if not rows:
            return pd.DataFrame(columns=pd.Index(_PROJECTED_COLUMNS))
        return pd.DataFrame.from_records(
            rows,
            columns=_PROJECTED_COLUMNS,
        ).drop_duplicates(ignore_index=True)

    @staticmethod
    def _build_projection_summary(
        *,
        reference_kinase_substrate_map: pd.DataFrame,
        display_lookup: dict[str, list[str]],
        site_key_to_display_id: dict[str, str],
        projected_map: pd.DataFrame,
        ambiguity_diagnostics: tuple[KinaseReferenceDisplayAmbiguityDiagnostic, ...],
    ) -> KinaseReferenceProjectionSummary:
        source_identifiers = _source_substrate_identifiers(
            reference_kinase_substrate_map
        )
        matched_source_identifiers: list[str] = []
        unmatched_source_identifiers: list[str] = []
        source_identifier_kinds: list[str] = []
        for substrate_identifier in source_identifiers:
            if substrate_identifier in site_key_to_display_id:
                matched_source_identifiers.append(substrate_identifier)
                source_identifier_kinds.append(_SOURCE_IDENTIFIER_KIND_DATASET_SITE_KEY)
                continue
            if display_lookup.get(substrate_identifier):
                matched_source_identifiers.append(substrate_identifier)
                source_identifier_kinds.append(
                    _SOURCE_IDENTIFIER_KIND_DATASET_DISPLAY_ID
                )
                continue
            unmatched_source_identifiers.append(substrate_identifier)
            source_identifier_kinds.append(
                _SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE
            )
        projected_site_keys: set[str] = set()
        if "substrate_site" in projected_map.columns:
            projected_site_keys = set(
                str(value) for value in projected_map.loc[:, "substrate_site"].tolist()
            )
        return KinaseReferenceProjectionSummary(
            source_reference_row_count=int(reference_kinase_substrate_map.shape[0]),
            matched_source_substrate_identifiers=tuple(
                sorted(matched_source_identifiers)
            ),
            unmatched_source_substrate_identifiers=tuple(
                sorted(unmatched_source_identifiers)
            ),
            projected_dataset_site_key_count=len(projected_site_keys),
            source_identifier_kinds=tuple(sorted(set(source_identifier_kinds))),
            one_to_many_display_reference_match_count=len(ambiguity_diagnostics),
            one_to_many_display_reference_site_key_rows=sum(
                item.matched_row_count for item in ambiguity_diagnostics
            ),
        )

    @staticmethod
    def _build_ambiguity_diagnostics(
        *,
        ambiguity_rows: dict[str, _AmbiguousDisplayReferenceRows],
        projected_map: pd.DataFrame | None,
    ) -> tuple[KinaseReferenceDisplayAmbiguityDiagnostic, ...]:
        diagnostics: list[KinaseReferenceDisplayAmbiguityDiagnostic] = []
        for display_id in sorted(ambiguity_rows):
            entry = ambiguity_rows[display_id]
            site_keys = tuple(str(value) for value in entry.site_keys)
            reference_rows = tuple(dict(row) for row in entry.reference_rows)
            projected_rows = 0
            if projected_map is not None and "display_id" in projected_map.columns:
                projected_rows = int(
                    (
                        projected_map.loc[:, "display_id"].astype(str)
                        == str(display_id)
                    ).sum()
                )
            diagnostics.append(
                KinaseReferenceDisplayAmbiguityDiagnostic(
                    display_id=str(display_id),
                    site_keys=site_keys,
                    reference_rows=reference_rows,
                    projected_rows=projected_rows,
                )
            )
        return tuple(diagnostics)

    @staticmethod
    def _raise_display_ambiguity_error(
        diagnostics: tuple[KinaseReferenceDisplayAmbiguityDiagnostic, ...],
    ) -> NoReturn:
        payload = [item.to_payload() for item in diagnostics]
        raise WorkflowBoundaryError(
            seam="kinase.interpreter.reference_display_ambiguity",
            next_action=(
                "provide protein-scoped site_key reference identity for ambiguous "
                "kinase-substrate evidence, or explicitly set "
                "reference_display_ambiguity_policy='allow_with_diagnostics' on "
                "KinaseWorkflowRequest if one-to-many display projection is intended"
            ),
            details={
                "ambiguous_display_ids": [item.display_id for item in diagnostics],
                "ambiguity_count": len(diagnostics),
                "ambiguity_diagnostics": payload,
            },
            message_prefix="kinase workflow boundary validation failed",
        )


__all__ = [
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ALLOW_WITH_DIAGNOSTICS",
    "KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICY_ERROR",
    "KinaseReferenceDisplayAmbiguityDiagnostic",
    "KinaseReferenceDisplayAmbiguityPolicy",
    "KinaseReferenceProjectionResult",
    "KinaseReferenceProjectionSummary",
    "KinaseReferenceProjector",
]


def _source_substrate_identifiers(
    reference_kinase_substrate_map: pd.DataFrame,
) -> tuple[str, ...]:
    if "substrate_site" not in reference_kinase_substrate_map.columns:
        return ()
    return tuple(
        sorted(
            {
                str(value)
                for value in reference_kinase_substrate_map.loc[
                    :, "substrate_site"
                ].tolist()
            }
        )
    )
