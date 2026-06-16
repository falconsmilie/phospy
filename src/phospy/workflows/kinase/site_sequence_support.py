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

_SITE_SEQUENCE_COLUMN = "site_sequence"
_INTERPRETER_VERSION = "phospy.workflows.kinase.site_sequence_support_builder.v1"


@dataclass(frozen=True, slots=True)
class KinaseSiteSequenceConflictDiagnostic:
    site_key: str
    display_id: str
    dataset_sequence: str
    reference_sequence: str
    selected_policy: KinaseSiteSequenceConflictPolicy
    selected_sequence: str
    selected_sequence_source: Literal["reference", "dataset"]
    interpreter_version: str = _INTERPRETER_VERSION

    def to_payload(self) -> dict[str, str]:
        return {
            "site_key": self.site_key,
            "display_id": self.display_id,
            "dataset_sequence": self.dataset_sequence,
            "reference_sequence": self.reference_sequence,
            "selected_policy": self.selected_policy,
            "selected_sequence": self.selected_sequence,
            "selected_sequence_source": self.selected_sequence_source,
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
            "interpreter_version": self.interpreter_version,
        }


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
        if conflict_policy not in KINASE_SITE_SEQUENCE_CONFLICT_POLICIES:
            supported = ", ".join(sorted(KINASE_SITE_SEQUENCE_CONFLICT_POLICIES))
            raise WorkflowBoundaryError(
                "kinase workflow boundary validation failed at seam="
                "kinase.interpreter.site_sequence_conflict_policy; "
                f"site_sequence_conflict_policy must be one of: {supported}; "
                "next_action=use a supported site-sequence conflict policy"
            )
        merged = pd.DataFrame(columns=[_SITE_SEQUENCE_COLUMN, "display_id"]).astype(
            "object"
        )
        merged.index = pd.Index([], dtype="object", name=dataset.index.name)
        dataset_sequences_available = 0
        dataset_sequences_missing = 0
        dataset_sequences_added = 0
        conflicts: list[KinaseSiteSequenceConflictDiagnostic] = []
        display_reference_multi_matches: list[dict[str, object]] = []
        if _SITE_SEQUENCE_COLUMN not in site_metadata.columns:
            return KinaseSiteSequenceSupportResult(
                site_sequences=merged,
                dataset_sequences_added=0,
                dataset_sequences_missing=int(dataset.shape[0]),
                dataset_sequences_available=0,
                conflict_policy=conflict_policy,
                conflicts=(),
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
        reference_sequences = reference_site_sequences.loc[:, _SITE_SEQUENCE_COLUMN]
        display_to_site_keys: dict[str, list[str]] = {}
        for site_id in dataset.index.tolist():
            site_key = str(site_id)
            display_value = display_series.loc[site_id]
            display_id = (
                site_key
                if not bool(pd.notna(display_value)) or str(display_value) == ""
                else str(display_value)
            )
            display_to_site_keys.setdefault(display_id, []).append(site_key)
            if display_id in reference_sequences.index:
                merged.loc[site_key, _SITE_SEQUENCE_COLUMN] = str(
                    reference_sequences.at[display_id]
                )
                merged.loc[site_key, "display_id"] = display_id
        for display_id, site_keys in display_to_site_keys.items():
            if display_id not in reference_sequences.index or len(site_keys) < 2:
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
            display_value = display_series.loc[site_id]
            display_id = (
                site_key
                if not bool(pd.notna(display_value)) or str(display_value) == ""
                else str(display_value)
            )
            sequence_value = dataset_sequence_series.loc[site_id]
            has_sequence = bool(pd.notna(sequence_value)) and str(sequence_value) != ""
            if not has_sequence:
                dataset_sequences_missing += 1
                continue
            dataset_sequences_available += 1
            dataset_sequence = str(sequence_value)
            if site_key not in merged.index:
                merged.loc[site_key, _SITE_SEQUENCE_COLUMN] = dataset_sequence
                merged.loc[site_key, "display_id"] = display_id
                dataset_sequences_added += 1
                continue
            reference_sequence = str(merged.at[site_key, _SITE_SEQUENCE_COLUMN])
            if reference_sequence == dataset_sequence:
                continue
            selected_sequence = reference_sequence
            selected_source: Literal["reference", "dataset"] = "reference"
            if conflict_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR:
                selected_sequence = reference_sequence
                selected_source = "reference"
            elif conflict_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET:
                selected_sequence = dataset_sequence
                selected_source = "dataset"
                merged.loc[site_key, _SITE_SEQUENCE_COLUMN] = dataset_sequence
            elif (
                conflict_policy == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE
            ):
                selected_sequence = reference_sequence
                selected_source = "reference"
            conflicts.append(
                KinaseSiteSequenceConflictDiagnostic(
                    site_key=site_key,
                    display_id=display_id,
                    dataset_sequence=dataset_sequence,
                    reference_sequence=reference_sequence,
                    selected_policy=conflict_policy,
                    selected_sequence=selected_sequence,
                    selected_sequence_source=selected_source,
                )
            )
        merged.index.name = dataset.index.name
        return KinaseSiteSequenceSupportResult(
            site_sequences=merged,
            dataset_sequences_added=dataset_sequences_added,
            dataset_sequences_missing=dataset_sequences_missing,
            dataset_sequences_available=dataset_sequences_available,
            conflict_policy=conflict_policy,
            conflicts=tuple(conflicts),
            display_reference_multi_matches=tuple(display_reference_multi_matches),
        )


__all__ = [
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICIES",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET",
    "KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE",
    "KinaseSiteSequenceConflictPolicy",
    "KinaseSiteSequenceSupportBuilder",
    "KinaseSiteSequenceSupportResult",
]
