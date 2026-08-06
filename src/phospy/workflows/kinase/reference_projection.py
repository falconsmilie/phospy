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
from phospy.errors.workflows import WorkflowBoundaryError

_INTERPRETER_VERSION = "phospy.workflows.kinase.reference_projector.v1"
_PROJECTION_SUMMARY_SCHEMA_VERSION = 1
_PROJECTED_COLUMNS = ["kinase", "substrate_site", "display_id"]
_UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT = 5
_SOURCE_IDENTIFIER_NAMESPACE = "references.kinase_substrate_map.substrate_site"
_OUTPUT_IDENTIFIER_NAMESPACE = "dataset.site_key"
_SOURCE_IDENTITY_SEMANTICS = (
    "source reference substrate identifiers before dataset projection; values may "
    "be dataset site_key values, dataset display_id values, or unmatched reference "
    "substrate identifiers"
)
_OUTPUT_IDENTITY_SEMANTICS = (
    "dataset site_key rows produced by projecting source reference substrate "
    "identifiers through dataset site_key/display_id identity"
)
_SOURCE_IDENTIFIER_KIND_DATASET_SITE_KEY = "dataset_site_key"
_SOURCE_IDENTIFIER_KIND_DATASET_DISPLAY_ID = "dataset_display_id"
_SOURCE_IDENTIFIER_KIND_UNMATCHED_REFERENCE_SUBSTRATE = (
    "unmatched_reference_substrate_identifier"
)


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
class KinaseReferenceProjectionSummary:
    """Typed reference-substrate projection summary before unmatched rows vanish."""

    source_reference_row_count: int
    matched_source_substrate_identifiers: tuple[str, ...]
    unmatched_source_substrate_identifiers: tuple[str, ...]
    projected_dataset_site_key_count: int
    source_identifier_kinds: tuple[str, ...]
    one_to_many_display_reference_match_count: int
    one_to_many_display_reference_site_key_rows: int
    source_identifier_namespace: str = _SOURCE_IDENTIFIER_NAMESPACE
    output_identifier_namespace: str = _OUTPUT_IDENTIFIER_NAMESPACE
    source_identity_semantics: str = _SOURCE_IDENTITY_SEMANTICS
    output_identity_semantics: str = _OUTPUT_IDENTITY_SEMANTICS
    unmatched_identifier_example_limit: int = _UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT
    schema_version: int = _PROJECTION_SUMMARY_SCHEMA_VERSION
    interpreter_version: str = _INTERPRETER_VERSION

    def __post_init__(self) -> None:
        source_reference_row_count = _require_non_negative_int(
            self.source_reference_row_count,
            field_name="kinase_reference_projection.source_reference_row_count",
        )
        projected_dataset_site_key_count = _require_non_negative_int(
            self.projected_dataset_site_key_count,
            field_name="kinase_reference_projection.projected_dataset_site_key_count",
        )
        one_to_many_display_reference_match_count = _require_non_negative_int(
            self.one_to_many_display_reference_match_count,
            field_name=(
                "kinase_reference_projection.one_to_many_display_reference_match_count"
            ),
        )
        one_to_many_display_reference_site_key_rows = _require_non_negative_int(
            self.one_to_many_display_reference_site_key_rows,
            field_name=(
                "kinase_reference_projection."
                "one_to_many_display_reference_site_key_rows"
            ),
        )
        example_limit = _require_non_negative_int(
            self.unmatched_identifier_example_limit,
            field_name=(
                "kinase_reference_projection.unmatched_identifier_example_limit"
            ),
        )
        schema_version = _require_non_negative_int(
            self.schema_version,
            field_name="kinase_reference_projection.schema_version",
        )
        matched = _unique_sorted_text_tuple(
            self.matched_source_substrate_identifiers,
            field_name=(
                "kinase_reference_projection.matched_source_substrate_identifiers"
            ),
        )
        unmatched = _unique_sorted_text_tuple(
            self.unmatched_source_substrate_identifiers,
            field_name=(
                "kinase_reference_projection.unmatched_source_substrate_identifiers"
            ),
        )
        overlap = set(matched).intersection(unmatched)
        if overlap:
            examples = sorted(overlap)[:_UNMATCHED_IDENTIFIER_EXAMPLE_LIMIT]
            raise WorkflowBoundaryError(
                "kinase reference projection summary matched and unmatched "
                f"identifier sets must be disjoint; overlap_examples={examples}"
            )
        source_identifier_kinds = _unique_sorted_text_tuple(
            self.source_identifier_kinds,
            field_name="kinase_reference_projection.source_identifier_kinds",
        )
        object.__setattr__(
            self,
            "source_reference_row_count",
            source_reference_row_count,
        )
        object.__setattr__(self, "matched_source_substrate_identifiers", matched)
        object.__setattr__(self, "unmatched_source_substrate_identifiers", unmatched)
        object.__setattr__(
            self,
            "projected_dataset_site_key_count",
            projected_dataset_site_key_count,
        )
        object.__setattr__(self, "source_identifier_kinds", source_identifier_kinds)
        object.__setattr__(
            self,
            "one_to_many_display_reference_match_count",
            one_to_many_display_reference_match_count,
        )
        object.__setattr__(
            self,
            "one_to_many_display_reference_site_key_rows",
            one_to_many_display_reference_site_key_rows,
        )
        object.__setattr__(
            self,
            "source_identifier_namespace",
            _require_non_empty_text(
                self.source_identifier_namespace,
                field_name="kinase_reference_projection.source_identifier_namespace",
            ),
        )
        object.__setattr__(
            self,
            "output_identifier_namespace",
            _require_non_empty_text(
                self.output_identifier_namespace,
                field_name="kinase_reference_projection.output_identifier_namespace",
            ),
        )
        object.__setattr__(
            self,
            "source_identity_semantics",
            _require_non_empty_text(
                self.source_identity_semantics,
                field_name="kinase_reference_projection.source_identity_semantics",
            ),
        )
        object.__setattr__(
            self,
            "output_identity_semantics",
            _require_non_empty_text(
                self.output_identity_semantics,
                field_name="kinase_reference_projection.output_identity_semantics",
            ),
        )
        object.__setattr__(
            self,
            "unmatched_identifier_example_limit",
            example_limit,
        )
        object.__setattr__(self, "schema_version", schema_version)
        object.__setattr__(
            self,
            "interpreter_version",
            _require_non_empty_text(
                self.interpreter_version,
                field_name="kinase_reference_projection.interpreter_version",
            ),
        )

    @property
    def unique_source_substrate_identifier_count(self) -> int:
        return int(
            len(self.matched_source_substrate_identifiers)
            + len(self.unmatched_source_substrate_identifiers)
        )

    @property
    def matched_source_substrate_identifier_count(self) -> int:
        return len(self.matched_source_substrate_identifiers)

    @property
    def unmatched_source_substrate_identifier_count(self) -> int:
        return len(self.unmatched_source_substrate_identifiers)

    @property
    def unmatched_source_substrate_identifier_examples(self) -> tuple[str, ...]:
        return self.unmatched_source_substrate_identifiers[
            : self.unmatched_identifier_example_limit
        ]

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": int(self.schema_version),
            "source_reference_row_count": int(self.source_reference_row_count),
            "unique_source_substrate_identifier_count": int(
                self.unique_source_substrate_identifier_count
            ),
            "matched_source_substrate_identifier_count": int(
                self.matched_source_substrate_identifier_count
            ),
            "unmatched_source_substrate_identifier_count": int(
                self.unmatched_source_substrate_identifier_count
            ),
            "matched_source_substrate_identifiers": list(
                self.matched_source_substrate_identifiers
            ),
            "unmatched_source_substrate_identifiers": list(
                self.unmatched_source_substrate_identifiers
            ),
            "unmatched_source_substrate_identifier_examples": list(
                self.unmatched_source_substrate_identifier_examples
            ),
            "unmatched_identifier_example_limit": int(
                self.unmatched_identifier_example_limit
            ),
            "projected_dataset_site_key_count": int(
                self.projected_dataset_site_key_count
            ),
            "source_identifier_namespace": self.source_identifier_namespace,
            "output_identifier_namespace": self.output_identifier_namespace,
            "source_identifier_kinds": list(self.source_identifier_kinds),
            "source_identity_semantics": self.source_identity_semantics,
            "output_identity_semantics": self.output_identity_semantics,
            "one_to_many_display_reference_match_count": int(
                self.one_to_many_display_reference_match_count
            ),
            "one_to_many_display_reference_site_key_rows": int(
                self.one_to_many_display_reference_site_key_rows
            ),
            "one_to_many_projection_diagnostics": (
                "display_reference_matching.one_to_many_display_reference_matches"
            ),
            "interpreter_version": self.interpreter_version,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> KinaseReferenceProjectionSummary:
        return cls(
            source_reference_row_count=_payload_int(
                payload,
                "source_reference_row_count",
            ),
            matched_source_substrate_identifiers=_payload_text_tuple(
                payload,
                "matched_source_substrate_identifiers",
            ),
            unmatched_source_substrate_identifiers=_payload_text_tuple(
                payload,
                "unmatched_source_substrate_identifiers",
            ),
            projected_dataset_site_key_count=_payload_int(
                payload,
                "projected_dataset_site_key_count",
            ),
            source_identifier_kinds=_payload_text_tuple(
                payload,
                "source_identifier_kinds",
            ),
            one_to_many_display_reference_match_count=_payload_int(
                payload,
                "one_to_many_display_reference_match_count",
            ),
            one_to_many_display_reference_site_key_rows=_payload_int(
                payload,
                "one_to_many_display_reference_site_key_rows",
            ),
            source_identifier_namespace=_payload_text(
                payload,
                "source_identifier_namespace",
            ),
            output_identifier_namespace=_payload_text(
                payload,
                "output_identifier_namespace",
            ),
            source_identity_semantics=_payload_text(
                payload,
                "source_identity_semantics",
            ),
            output_identity_semantics=_payload_text(
                payload,
                "output_identity_semantics",
            ),
            unmatched_identifier_example_limit=_payload_int(
                payload,
                "unmatched_identifier_example_limit",
            ),
            schema_version=_payload_int(payload, "schema_version"),
            interpreter_version=_payload_text(payload, "interpreter_version"),
        )


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


def _unique_sorted_text_tuple(
    values: tuple[str, ...],
    *,
    field_name: str,
) -> tuple[str, ...]:
    return tuple(
        sorted(
            {_require_non_empty_text(value, field_name=field_name) for value in values}
        )
    )


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        raise WorkflowBoundaryError(f"{field_name} must be a non-empty string")
    return text


def _require_non_negative_int(value: object, *, field_name: str) -> int:
    if isinstance(value, bool):
        raise WorkflowBoundaryError(f"{field_name} must be a non-negative integer")
    try:
        normalized = int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise WorkflowBoundaryError(
            f"{field_name} must be a non-negative integer"
        ) from exc
    if normalized < 0:
        raise WorkflowBoundaryError(f"{field_name} must be a non-negative integer")
    return normalized


def _payload_text(payload: Mapping[str, object], key: str) -> str:
    if key not in payload:
        raise WorkflowBoundaryError(
            f"kinase reference projection payload is missing required key: {key}"
        )
    return _require_non_empty_text(
        payload[key],
        field_name=f"kinase_reference_projection.{key}",
    )


def _payload_int(payload: Mapping[str, object], key: str) -> int:
    if key not in payload:
        raise WorkflowBoundaryError(
            f"kinase reference projection payload is missing required key: {key}"
        )
    return _require_non_negative_int(
        payload[key],
        field_name=f"kinase_reference_projection.{key}",
    )


def _payload_text_tuple(payload: Mapping[str, object], key: str) -> tuple[str, ...]:
    if key not in payload:
        raise WorkflowBoundaryError(
            f"kinase reference projection payload is missing required key: {key}"
        )
    value = payload[key]
    if isinstance(value, str) or not isinstance(value, (list, tuple)):
        raise WorkflowBoundaryError(
            f"kinase_reference_projection.{key} must be a sequence of strings"
        )
    return tuple(
        _require_non_empty_text(
            item,
            field_name=f"kinase_reference_projection.{key}[]",
        )
        for item in value
    )
