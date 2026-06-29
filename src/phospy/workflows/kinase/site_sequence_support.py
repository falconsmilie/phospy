"""Site-sequence support assembly for kinase workflow interpretation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal, cast

import pandas as pd

from phospy.contracts.configs.kinase import (
    KINASE_SITE_SEQUENCE_CONFLICT_POLICIES,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    KinaseSiteSequenceConflictPolicy,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.workflows.kinase.site_sequence_policy import (
    resolve_site_sequence_conflict_policy,
)

_SITE_SEQUENCE_COLUMN = "site_sequence"
_INTERPRETER_VERSION = "phospy.workflows.kinase.site_sequence_support_builder.v1"
_CONFLICT_NEXT_ACTION = (
    "fix dataset site_sequence values for conflicting sites or use "
    "site_sequence_conflict_policy='prefer_reference' or 'prefer_dataset' on "
    "KinaseWorkflowRequest"
)


@dataclass(frozen=True, slots=True)
class KinaseSiteSequenceConflictDiagnostic:
    site_key: str
    display_id: str
    dataset_sequence: str
    reference_sequence: str
    selected_policy: KinaseSiteSequenceConflictPolicy
    selected_sequence: str | None
    selected_sequence_source: Literal["reference", "dataset", "unresolved"]
    diagnostic: str
    interpreter_version: str = _INTERPRETER_VERSION

    def to_payload(self) -> dict[str, object]:
        return {
            "site_key": self.site_key,
            "display_id": self.display_id,
            "dataset_sequence": self.dataset_sequence,
            "reference_sequence": self.reference_sequence,
            "policy": self.selected_policy,
            "selected_policy": self.selected_policy,
            "selected_sequence": self.selected_sequence,
            "selected_sequence_source": self.selected_sequence_source,
            "diagnostic": self.diagnostic,
            "interpreter_version": self.interpreter_version,
        }


@dataclass(frozen=True, slots=True)
class KinaseSiteSequenceSelectionDiagnostic:
    site_key: str
    display_id: str
    selected_sequence: str
    selected_sequence_source: Literal["reference", "dataset"]
    policy: KinaseSiteSequenceConflictPolicy
    dataset_sequence: str | None
    reference_sequence: str | None
    diagnostic: str
    interpreter_version: str = _INTERPRETER_VERSION

    def to_payload(self) -> dict[str, object]:
        return {
            "site_key": self.site_key,
            "display_id": self.display_id,
            "selected_sequence": self.selected_sequence,
            "selected_sequence_source": self.selected_sequence_source,
            "policy": self.policy,
            "dataset_sequence": self.dataset_sequence,
            "reference_sequence": self.reference_sequence,
            "diagnostic": self.diagnostic,
            "interpreter_version": self.interpreter_version,
        }


@dataclass(frozen=True, slots=True)
class KinaseSiteSequenceSupportResult:
    site_sequences: pd.DataFrame
    dataset_sequences_added: int
    dataset_sequences_missing: int
    dataset_sequences_available: int
    conflict_policy: KinaseSiteSequenceConflictPolicy
    conflicts: tuple[KinaseSiteSequenceConflictDiagnostic, ...]
    display_reference_multi_matches: tuple[dict[str, object], ...]
    sequence_source_records: tuple[KinaseSiteSequenceSelectionDiagnostic, ...] = ()
    interpreter_version: str = _INTERPRETER_VERSION

    @property
    def dataset_reference_conflict_count(self) -> int:
        return len(self.conflicts)

    def diagnostics_payload(self) -> dict[str, object]:
        return {
            "dataset_sequences_added": int(self.dataset_sequences_added),
            "dataset_reference_conflict_count": int(
                self.dataset_reference_conflict_count
            ),
            "dataset_sequences_missing": int(self.dataset_sequences_missing),
            "dataset_sequences_available": int(self.dataset_sequences_available),
            "display_reference_multi_match_count": int(
                len(self.display_reference_multi_matches)
            ),
            "display_reference_multi_matches": [
                {
                    "display_id": str(item["display_id"]),
                    "site_keys": tuple(
                        str(value)
                        for value in cast(Sequence[object], item["site_keys"])
                    ),
                }
                for item in self.display_reference_multi_matches
            ],
            "conflict_policy": self.conflict_policy,
            "conflict_diagnostics": [item.to_payload() for item in self.conflicts],
            "selected_sequence_sources": [
                item.to_payload() for item in self.sequence_source_records
            ],
            "interpreter_version": self.interpreter_version,
        }


class KinaseSiteSequenceConflictError(WorkflowBoundaryError):
    """Dataset and reference sequence support disagree under the error policy."""


class KinaseSiteSequenceSupportBuilder:
    """Build execution-time site-sequence support for kinase workflow scoring."""

    def run(
        self,
        *,
        dataset: pd.DataFrame,
        site_metadata: pd.DataFrame,
        reference_site_sequences: pd.DataFrame,
        conflict_policy: KinaseSiteSequenceConflictPolicy,
    ) -> KinaseSiteSequenceSupportResult:
        resolved_policy = resolve_site_sequence_conflict_policy(
            conflict_policy,
            field_name="site_sequence_conflict_policy",
            error_type=WorkflowBoundaryError,
        )
        merged = pd.DataFrame(columns=[_SITE_SEQUENCE_COLUMN, "display_id"]).astype(
            "object"
        )
        merged.index = pd.Index([], dtype="object", name=dataset.index.name)
        dataset_sequences_available = 0
        dataset_sequences_missing = 0
        dataset_sequences_added = 0
        conflicts: list[KinaseSiteSequenceConflictDiagnostic] = []
        sequence_source_records: list[KinaseSiteSequenceSelectionDiagnostic] = []
        display_reference_multi_matches: list[dict[str, object]] = []
        if _SITE_SEQUENCE_COLUMN not in site_metadata.columns:
            return KinaseSiteSequenceSupportResult(
                site_sequences=merged,
                dataset_sequences_added=0,
                dataset_sequences_missing=int(dataset.shape[0]),
                dataset_sequences_available=0,
                conflict_policy=resolved_policy,
                conflicts=(),
                sequence_source_records=(),
                display_reference_multi_matches=(),
            )
        display_series = (
            site_metadata.reindex(dataset.index)
            .loc[:, "display_id"]
            .astype("string")
            .str.strip()
            if "display_id" in site_metadata.columns
            else pd.Series(
                dataset.index.astype(str), index=dataset.index, dtype="string"
            )
        )
        reference_sequence_by_display_id = _reference_sequence_by_display_id(
            reference_site_sequences
        )
        reference_sequence_by_site_key: dict[str, str] = {}
        display_id_by_site_key: dict[str, str] = {}
        display_to_site_keys: dict[str, list[str]] = {}
        for site_id in dataset.index.tolist():
            site_key = str(site_id)
            display_value = display_series.loc[site_id]
            display_id = (
                site_key
                if not bool(pd.notna(display_value)) or str(display_value) == ""
                else str(display_value)
            )
            display_id_by_site_key[site_key] = display_id
            display_to_site_keys.setdefault(display_id, []).append(site_key)
            reference_sequence = reference_sequence_by_display_id.get(display_id)
            if reference_sequence is not None:
                reference_sequence_by_site_key[site_key] = reference_sequence
                merged.loc[site_key, _SITE_SEQUENCE_COLUMN] = reference_sequence
                merged.loc[site_key, "display_id"] = display_id
        for display_id, site_keys in display_to_site_keys.items():
            if display_id not in reference_sequence_by_display_id or len(site_keys) < 2:
                continue
            display_reference_multi_matches.append(
                {
                    "display_id": display_id,
                    "site_keys": tuple(site_keys),
                }
            )
        dataset_sequence_series = (
            site_metadata.reindex(dataset.index)
            .loc[:, _SITE_SEQUENCE_COLUMN]
            .astype("string")
            .str.strip()
        )
        for site_id in dataset.index.tolist():
            site_key = str(site_id)
            display_id = display_id_by_site_key[site_key]
            sequence_value = dataset_sequence_series.loc[site_id]
            dataset_sequence = _normalise_sequence_value(sequence_value)
            reference_sequence = reference_sequence_by_site_key.get(site_key)
            has_sequence = dataset_sequence is not None
            if not has_sequence:
                dataset_sequences_missing += 1
                if reference_sequence is not None:
                    sequence_source_records.append(
                        KinaseSiteSequenceSelectionDiagnostic(
                            site_key=site_key,
                            display_id=display_id,
                            selected_sequence=reference_sequence,
                            selected_sequence_source="reference",
                            policy=resolved_policy,
                            dataset_sequence=None,
                            reference_sequence=reference_sequence,
                            diagnostic=(
                                "reference sequence selected because dataset "
                                "sequence is missing"
                            ),
                        )
                    )
                continue
            dataset_sequences_available += 1
            if reference_sequence is None:
                merged.loc[site_key, _SITE_SEQUENCE_COLUMN] = dataset_sequence
                merged.loc[site_key, "display_id"] = display_id
                dataset_sequences_added += 1
                sequence_source_records.append(
                    KinaseSiteSequenceSelectionDiagnostic(
                        site_key=site_key,
                        display_id=display_id,
                        selected_sequence=dataset_sequence,
                        selected_sequence_source="dataset",
                        policy=resolved_policy,
                        dataset_sequence=dataset_sequence,
                        reference_sequence=None,
                        diagnostic=(
                            "dataset sequence selected because reference "
                            "sequence is missing"
                        ),
                    )
                )
                continue
            if reference_sequence == dataset_sequence:
                sequence_source_records.append(
                    KinaseSiteSequenceSelectionDiagnostic(
                        site_key=site_key,
                        display_id=display_id,
                        selected_sequence=reference_sequence,
                        selected_sequence_source="reference",
                        policy=resolved_policy,
                        dataset_sequence=dataset_sequence,
                        reference_sequence=reference_sequence,
                        diagnostic=(
                            "reference sequence selected; dataset sequence "
                            "matches reference sequence"
                        ),
                    )
                )
                continue
            selected_sequence: str | None = None
            selected_source: Literal["reference", "dataset", "unresolved"] = (
                "unresolved"
            )
            diagnostic = (
                "dataset/reference sequence conflict is unresolved under "
                "site_sequence_conflict_policy='error'"
            )
            if resolved_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET:
                selected_sequence = dataset_sequence
                selected_source = "dataset"
                merged.loc[site_key, _SITE_SEQUENCE_COLUMN] = dataset_sequence
                diagnostic = (
                    "dataset sequence selected for dataset/reference sequence "
                    "conflict under site_sequence_conflict_policy='prefer_dataset'"
                )
                sequence_source_records.append(
                    KinaseSiteSequenceSelectionDiagnostic(
                        site_key=site_key,
                        display_id=display_id,
                        selected_sequence=dataset_sequence,
                        selected_sequence_source="dataset",
                        policy=resolved_policy,
                        dataset_sequence=dataset_sequence,
                        reference_sequence=reference_sequence,
                        diagnostic=diagnostic,
                    )
                )
            elif (
                resolved_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE
            ):
                selected_sequence = reference_sequence
                selected_source = "reference"
                diagnostic = (
                    "reference sequence selected for dataset/reference sequence "
                    "conflict under site_sequence_conflict_policy='prefer_reference'"
                )
                sequence_source_records.append(
                    KinaseSiteSequenceSelectionDiagnostic(
                        site_key=site_key,
                        display_id=display_id,
                        selected_sequence=reference_sequence,
                        selected_sequence_source="reference",
                        policy=resolved_policy,
                        dataset_sequence=dataset_sequence,
                        reference_sequence=reference_sequence,
                        diagnostic=diagnostic,
                    )
                )
            conflicts.append(
                KinaseSiteSequenceConflictDiagnostic(
                    site_key=site_key,
                    display_id=display_id,
                    dataset_sequence=dataset_sequence,
                    reference_sequence=reference_sequence,
                    selected_policy=resolved_policy,
                    selected_sequence=selected_sequence,
                    selected_sequence_source=selected_source,
                    diagnostic=diagnostic,
                )
            )
        if resolved_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR and conflicts:
            _raise_site_sequence_conflict_error(
                conflict_policy=resolved_policy,
                conflicts=tuple(conflicts),
            )
        merged.index.name = dataset.index.name
        return KinaseSiteSequenceSupportResult(
            site_sequences=merged,
            dataset_sequences_added=dataset_sequences_added,
            dataset_sequences_missing=dataset_sequences_missing,
            dataset_sequences_available=dataset_sequences_available,
            conflict_policy=resolved_policy,
            conflicts=tuple(conflicts),
            sequence_source_records=tuple(sequence_source_records),
            display_reference_multi_matches=tuple(display_reference_multi_matches),
        )


def _reference_sequence_by_display_id(
    reference_site_sequences: pd.DataFrame,
) -> dict[str, str]:
    if _SITE_SEQUENCE_COLUMN not in reference_site_sequences.columns:
        return {}
    reference_sequences: dict[str, str] = {}
    for display_id, value in reference_site_sequences.loc[
        :, _SITE_SEQUENCE_COLUMN
    ].items():
        display_key = str(display_id).strip()
        sequence = _normalise_sequence_value(value)
        if display_key == "" or sequence is None:
            continue
        reference_sequences[display_key] = sequence
    return reference_sequences


def _normalise_sequence_value(value: object) -> str | None:
    if bool(pd.Series((value,), dtype="object").isna().iat[0]):
        return None
    text = str(value).strip()
    return text or None


def _raise_site_sequence_conflict_error(
    *,
    conflict_policy: KinaseSiteSequenceConflictPolicy,
    conflicts: tuple[KinaseSiteSequenceConflictDiagnostic, ...],
) -> None:
    first = conflicts[0]
    raise KinaseSiteSequenceConflictError(
        seam="kinase.interpreter.site_sequence_conflict",
        next_action=_CONFLICT_NEXT_ACTION,
        details={
            "conflict_policy": conflict_policy,
            "dataset_reference_conflict_count": int(len(conflicts)),
            "site_key": first.site_key,
            "display_id": first.display_id,
            "dataset_sequence": first.dataset_sequence,
            "reference_sequence": first.reference_sequence,
            "conflict_diagnostics": [item.to_payload() for item in conflicts],
        },
        message_prefix="kinase workflow boundary validation failed",
    )


__all__ = [
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICIES",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE",
    "KinaseSiteSequenceConflictError",
    "KinaseSiteSequenceConflictDiagnostic",
    "KinaseSiteSequenceConflictPolicy",
    "KinaseSiteSequenceSelectionDiagnostic",
    "KinaseSiteSequenceSupportBuilder",
    "KinaseSiteSequenceSupportResult",
]
