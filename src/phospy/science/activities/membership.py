"""Typed substrate-membership provenance for activity inference."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Final

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance.immutability import (
    FrozenJsonMapping,
    freeze_json_mapping_with_error_type,
    thaw_json_mapping,
)
from phospy.provenance.models import TableFingerprint
from phospy.provenance.serialization import (
    table_fingerprint_from_payload,
    table_fingerprint_to_payload,
)

ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE: Final = "fixed_external_reference"
ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF: Final = "sequence_only_motif"
ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED: Final = "profile_derived"
ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF: Final = "fused_profile_motif"
ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED: Final = "prediction_selected"
ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN: Final = "unknown"

ACTIVITY_MEMBERSHIP_SOURCE_CATEGORIES: Final = frozenset(
    {
        ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE,
        ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF,
        ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED,
        ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF,
        ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED,
        ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN,
    }
)

ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION: Final = "1"

KSEA_MEMBERSHIP_INFERENCE_STATUS_ELIGIBLE: Final = "ordinary_p_q_eligible"
KSEA_MEMBERSHIP_INFERENCE_STATUS_UNAVAILABLE: Final = "ordinary_p_q_unavailable"
KSEA_MEMBERSHIP_ELIGIBLE_REASON: Final = (
    "membership_selection_independent_of_tested_quantitative_matrix"
)
KSEA_MEMBERSHIP_CONSUMED_TESTED_MATRIX_REASON: Final = (
    "membership_selection_consumed_tested_quantitative_matrix"
)
KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON: Final = (
    "missing_membership_selection_provenance"
)


@dataclass(frozen=True, slots=True)
class ActivityMembershipSelection:
    """Provenance and inferential eligibility for activity substrate membership.

    The activity science domain owns whether ordinary inferential p/q values are
    eligible. Workflow code may construct this record from resolved scoring and
    prediction provenance, but KSEA uses only this typed contract and never
    infers independence from raw matrices or method-name strings.
    """

    source_category: str
    selection_method: str
    selection_method_version: str
    score_source: str
    threshold_top_k_policy: Mapping[str, object] = field(default_factory=dict)
    source_reference_fingerprints: tuple[TableFingerprint, ...] = ()
    quantitative_dataset_fingerprint: TableFingerprint | None = None
    consumed_tested_matrix: bool = False
    selected_kinase_universe: tuple[str, ...] = ()
    selected_substrate_universe: tuple[str, ...] = ()
    inferential_eligible: bool = False
    inferential_eligibility_reason: str = KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON

    def __post_init__(self) -> None:
        source_category = str(self.source_category).strip()
        if source_category not in ACTIVITY_MEMBERSHIP_SOURCE_CATEGORIES:
            allowed = ", ".join(sorted(ACTIVITY_MEMBERSHIP_SOURCE_CATEGORIES))
            raise WorkflowBoundaryError(
                "activity membership source_category must be one of: "
                f"{allowed}; got {source_category!r}"
            )
        consumed_tested_matrix = _require_bool(
            self.consumed_tested_matrix,
            field_name="activity_membership_selection.consumed_tested_matrix",
        )
        inferential_eligible = _require_bool(
            self.inferential_eligible,
            field_name="activity_membership_selection.inferential_eligible",
        )
        if inferential_eligible and consumed_tested_matrix:
            raise WorkflowBoundaryError(
                "activity membership selection cannot be inferentially eligible "
                "when consumed_tested_matrix is true"
            )
        if consumed_tested_matrix and self.quantitative_dataset_fingerprint is None:
            raise WorkflowBoundaryError(
                "activity membership selection that consumed the tested matrix "
                "must record quantitative_dataset_fingerprint"
            )
        source_reference_fingerprints = _coerce_table_fingerprints(
            self.source_reference_fingerprints,
            field_name="activity_membership_selection.source_reference_fingerprints",
        )
        quantitative_fingerprint = self.quantitative_dataset_fingerprint
        if quantitative_fingerprint is not None and not isinstance(
            quantitative_fingerprint,
            TableFingerprint,
        ):
            raise WorkflowBoundaryError(
                "activity_membership_selection.quantitative_dataset_fingerprint "
                "must be TableFingerprint or None"
            )
        object.__setattr__(self, "source_category", source_category)
        object.__setattr__(
            self,
            "selection_method",
            _require_non_empty_text(
                self.selection_method,
                field_name="activity_membership_selection.selection_method",
            ),
        )
        object.__setattr__(
            self,
            "selection_method_version",
            _require_non_empty_text(
                self.selection_method_version,
                field_name="activity_membership_selection.selection_method_version",
            ),
        )
        object.__setattr__(
            self,
            "score_source",
            _require_non_empty_text(
                self.score_source,
                field_name="activity_membership_selection.score_source",
            ),
        )
        object.__setattr__(
            self,
            "threshold_top_k_policy",
            _freeze_policy_mapping(self.threshold_top_k_policy),
        )
        object.__setattr__(
            self,
            "source_reference_fingerprints",
            source_reference_fingerprints,
        )
        object.__setattr__(
            self,
            "selected_kinase_universe",
            _string_tuple(
                self.selected_kinase_universe,
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
        )
        object.__setattr__(
            self,
            "selected_substrate_universe",
            _string_tuple(
                self.selected_substrate_universe,
                field_name="activity_membership_selection.selected_substrate_universe",
            ),
        )
        object.__setattr__(self, "consumed_tested_matrix", consumed_tested_matrix)
        object.__setattr__(self, "inferential_eligible", inferential_eligible)
        object.__setattr__(
            self,
            "inferential_eligibility_reason",
            _require_non_empty_text(
                self.inferential_eligibility_reason,
                field_name=(
                    "activity_membership_selection.inferential_eligibility_reason"
                ),
            ),
        )

    @property
    def inferential_status(self) -> str:
        """Return a compact status token for p/q-value availability."""

        if self.inferential_eligible:
            return KSEA_MEMBERSHIP_INFERENCE_STATUS_ELIGIBLE
        return KSEA_MEMBERSHIP_INFERENCE_STATUS_UNAVAILABLE

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible membership-selection payload."""

        source_reference_fingerprints = [
            table_fingerprint_to_payload(item)
            for item in self.source_reference_fingerprints
        ]
        return {
            "source_category": self.source_category,
            "selection_method": self.selection_method,
            "selection_method_version": self.selection_method_version,
            "score_source": self.score_source,
            "threshold_top_k_policy": thaw_json_mapping(
                self.threshold_top_k_policy,
                field_name="activity_membership_selection.threshold_top_k_policy",
            ),
            "source_reference_fingerprint": (
                None
                if not source_reference_fingerprints
                else source_reference_fingerprints[0]
            ),
            "source_reference_fingerprints": source_reference_fingerprints,
            "quantitative_dataset_fingerprint": (
                None
                if self.quantitative_dataset_fingerprint is None
                else table_fingerprint_to_payload(self.quantitative_dataset_fingerprint)
            ),
            "consumed_tested_matrix": bool(self.consumed_tested_matrix),
            "selected_kinase_universe": list(self.selected_kinase_universe),
            "selected_substrate_universe": list(self.selected_substrate_universe),
            "inferential_eligible": bool(self.inferential_eligible),
            "inferential_status": self.inferential_status,
            "inferential_eligibility_reason": self.inferential_eligibility_reason,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMembershipSelection:
        """Reconstruct a membership-selection record from a payload."""

        if not isinstance(payload, Mapping):
            raise WorkflowBoundaryError(
                "activity membership selection payload must be a mapping"
            )
        source_reference_fingerprints = _fingerprints_from_payload(payload)
        quantitative_payload = payload.get("quantitative_dataset_fingerprint")
        quantitative_fingerprint = (
            None
            if quantitative_payload is None
            else table_fingerprint_from_payload(
                _require_mapping(
                    quantitative_payload,
                    field_name=(
                        "activity_membership_selection.quantitative_dataset_fingerprint"
                    ),
                )
            )
        )
        return cls(
            source_category=str(payload.get("source_category", "")),
            selection_method=str(payload.get("selection_method", "")),
            selection_method_version=str(payload.get("selection_method_version", "")),
            score_source=str(payload.get("score_source", "")),
            threshold_top_k_policy=_require_mapping(
                payload.get("threshold_top_k_policy", {}),
                field_name="activity_membership_selection.threshold_top_k_policy",
            ),
            source_reference_fingerprints=source_reference_fingerprints,
            quantitative_dataset_fingerprint=quantitative_fingerprint,
            consumed_tested_matrix=_require_bool(
                payload.get("consumed_tested_matrix", False),
                field_name=("activity_membership_selection.consumed_tested_matrix"),
            ),
            selected_kinase_universe=_string_tuple(
                payload.get("selected_kinase_universe", ()),
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
            selected_substrate_universe=_string_tuple(
                payload.get("selected_substrate_universe", ()),
                field_name=(
                    "activity_membership_selection.selected_substrate_universe"
                ),
            ),
            inferential_eligible=_require_bool(
                payload.get("inferential_eligible", False),
                field_name="activity_membership_selection.inferential_eligible",
            ),
            inferential_eligibility_reason=str(
                payload.get(
                    "inferential_eligibility_reason",
                    KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON,
                )
            ),
        )

    @classmethod
    def missing(
        cls,
        *,
        selected_kinase_universe: Iterable[object] = (),
        selected_substrate_universe: Iterable[object] = (),
    ) -> ActivityMembershipSelection:
        """Return an explicit ineligible record for missing provenance."""

        return cls(
            source_category=ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN,
            selection_method="missing_membership_selection_provenance",
            selection_method_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
            score_source="unknown",
            threshold_top_k_policy={},
            source_reference_fingerprints=(),
            quantitative_dataset_fingerprint=None,
            consumed_tested_matrix=False,
            selected_kinase_universe=_string_tuple(
                selected_kinase_universe,
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
            selected_substrate_universe=_string_tuple(
                selected_substrate_universe,
                field_name=(
                    "activity_membership_selection.selected_substrate_universe"
                ),
            ),
            inferential_eligible=False,
            inferential_eligibility_reason=KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON,
        )


def selected_substrate_universe_from_prediction_matrix(
    pred_mat: pd.DataFrame,
    *,
    threshold: float,
) -> tuple[str, ...]:
    """Return sites selected by finite prediction support at a threshold."""

    if pred_mat.empty:
        return ()
    values = pred_mat.to_numpy(dtype=float, copy=False)
    selected = np.isfinite(values) & (values >= float(threshold))
    if not bool(selected.any()):
        return ()
    site_mask = np.asarray(selected.any(axis=1), dtype=bool).reshape(-1)
    return tuple(
        str(site_id)
        for position, site_id in enumerate(pred_mat.index.tolist())
        if bool(site_mask[position])
    )


def _fingerprints_from_payload(
    payload: Mapping[str, object],
) -> tuple[TableFingerprint, ...]:
    raw_fingerprints = payload.get("source_reference_fingerprints")
    if raw_fingerprints is None:
        raw_fingerprint = payload.get("source_reference_fingerprint")
        if raw_fingerprint is None:
            return ()
        raw_fingerprints = [raw_fingerprint]
    if isinstance(raw_fingerprints, (str, bytes, bytearray)):
        raise WorkflowBoundaryError(
            "activity_membership_selection.source_reference_fingerprints "
            "must be a sequence"
        )
    try:
        values = tuple(raw_fingerprints)  # type: ignore[arg-type]
    except TypeError as exc:
        raise WorkflowBoundaryError(
            "activity_membership_selection.source_reference_fingerprints "
            "must be a sequence"
        ) from exc
    fingerprints: list[TableFingerprint] = []
    for position, value in enumerate(values):
        fingerprints.append(
            table_fingerprint_from_payload(
                _require_mapping(
                    value,
                    field_name=(
                        "activity_membership_selection."
                        f"source_reference_fingerprints[{position}]"
                    ),
                )
            )
        )
    return tuple(fingerprints)


def _coerce_table_fingerprints(
    value: object,
    *,
    field_name: str,
) -> tuple[TableFingerprint, ...]:
    if isinstance(value, (str, bytes, bytearray)):
        raise WorkflowBoundaryError(f"{field_name} must be a sequence")
    try:
        fingerprints = tuple(value)  # type: ignore[arg-type]
    except TypeError as exc:
        raise WorkflowBoundaryError(f"{field_name} must be a sequence") from exc
    for fingerprint in fingerprints:
        if not isinstance(fingerprint, TableFingerprint):
            raise WorkflowBoundaryError(
                f"{field_name} must contain only TableFingerprint values"
            )
    return fingerprints


def _freeze_policy_mapping(value: object) -> FrozenJsonMapping:
    if not isinstance(value, Mapping):
        raise WorkflowBoundaryError(
            "activity_membership_selection.threshold_top_k_policy must be a mapping"
        )
    return freeze_json_mapping_with_error_type(
        value,
        field_name="activity_membership_selection.threshold_top_k_policy",
        error_type=WorkflowBoundaryError,
    )


def _string_tuple(value: object, *, field_name: str) -> tuple[str, ...]:
    if isinstance(value, (str, bytes, bytearray)):
        raise WorkflowBoundaryError(f"{field_name} must be a sequence of strings")
    try:
        values = tuple(value)  # type: ignore[arg-type]
    except TypeError as exc:
        raise WorkflowBoundaryError(
            f"{field_name} must be a sequence of strings"
        ) from exc
    result: list[str] = []
    for item in values:
        text = str(item).strip()
        if not text:
            raise WorkflowBoundaryError(f"{field_name} must not contain blank labels")
        result.append(text)
    return tuple(result)


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        raise WorkflowBoundaryError(f"{field_name} must be a non-empty string")
    return text


def _require_bool(value: object, *, field_name: str) -> bool:
    if not isinstance(value, bool):
        raise WorkflowBoundaryError(f"{field_name} must be a bool")
    return value


def _require_mapping(value: object, *, field_name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise WorkflowBoundaryError(f"{field_name} must be a mapping")
    return {str(key): item for key, item in value.items()}


__all__ = [
    "ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION",
    "ACTIVITY_MEMBERSHIP_SOURCE_CATEGORIES",
    "ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE",
    "ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF",
    "ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED",
    "ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED",
    "ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF",
    "ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN",
    "ActivityMembershipSelection",
    "KSEA_MEMBERSHIP_CONSUMED_TESTED_MATRIX_REASON",
    "KSEA_MEMBERSHIP_ELIGIBLE_REASON",
    "KSEA_MEMBERSHIP_INFERENCE_STATUS_ELIGIBLE",
    "KSEA_MEMBERSHIP_INFERENCE_STATUS_UNAVAILABLE",
    "KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON",
    "selected_substrate_universe_from_prediction_matrix",
]
