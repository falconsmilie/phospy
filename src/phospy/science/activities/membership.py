"""Typed substrate-membership provenance for activity inference."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Final

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance.hashing import fingerprint_table_normalized_axes
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

ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION: Final = "2"
KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION: Final = "2"
KSEA_TESTED_QUANTITATIVE_MATRIX_FINGERPRINT_NAME: Final = (
    "dataset.ksea_background_phospho_matrix"
)
KSEA_SELECTION_QUANTITATIVE_MATRIX_FINGERPRINT_NAME: Final = (
    "dataset.scoring_phospho_matrix"
)
KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION: Final = "1"
KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_FIXED_EXTERNAL_REFERENCE: Final = (
    "fixed_external_reference_membership_independent_of_tested_matrix"
)
KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_SEQUENCE_ONLY_MOTIF: Final = (
    "sequence_only_motif_membership_independent_of_tested_matrix"
)

KSEA_MEMBERSHIP_INFERENCE_STATUS_ELIGIBLE: Final = "ordinary_p_q_available"
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
KSEA_MEMBERSHIP_PROFILE_DERIVED_REASON: Final = (
    "membership_selection_profile_derived_from_quantitative_data"
)
KSEA_MEMBERSHIP_FUSED_PROFILE_MOTIF_REASON: Final = (
    "membership_selection_fused_profile_motif_uses_quantitative_data"
)
KSEA_MEMBERSHIP_PREDICTION_SELECTED_REASON: Final = (
    "prediction_selected_membership_not_independently_eligible"
)
KSEA_MEMBERSHIP_INCOMPLETE_INDEPENDENCE_EVIDENCE_REASON: Final = (
    "incomplete_independent_membership_evidence"
)
KSEA_MEMBERSHIP_QUANTITATIVE_SELECTION_REASON: Final = (
    "membership_selection_consumed_quantitative_data"
)

_KSEA_SEQUENCE_ONLY_MOTIF_SCORE_SOURCE: Final = "kinase_library_motif_scores"


@dataclass(frozen=True, slots=True)
class KseaMembershipInferentialDecision:
    """Science-owned ordinary KSEA p/q-value availability decision."""

    ordinary_p_q_available: bool
    status: str
    reason: str
    missing_evidence: tuple[str, ...] = ()
    policy_version: str = KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION

    def __post_init__(self) -> None:
        ordinary_p_q_available = _require_bool(
            self.ordinary_p_q_available,
            field_name="ksea_membership_inferential_decision.ordinary_p_q_available",
        )
        status = _require_non_empty_text(
            self.status,
            field_name="ksea_membership_inferential_decision.status",
        )
        allowed_statuses = {
            KSEA_MEMBERSHIP_INFERENCE_STATUS_ELIGIBLE,
            KSEA_MEMBERSHIP_INFERENCE_STATUS_UNAVAILABLE,
        }
        if status not in allowed_statuses:
            allowed = ", ".join(sorted(allowed_statuses))
            raise WorkflowBoundaryError(
                "ksea membership inferential decision status must be one of: "
                f"{allowed}; got {status!r}"
            )
        if (
            ordinary_p_q_available
            and status != KSEA_MEMBERSHIP_INFERENCE_STATUS_ELIGIBLE
        ):
            raise WorkflowBoundaryError(
                "ksea membership inferential decision status must be "
                "ordinary_p_q_available when ordinary_p_q_available is true"
            )
        if (
            not ordinary_p_q_available
            and status != KSEA_MEMBERSHIP_INFERENCE_STATUS_UNAVAILABLE
        ):
            raise WorkflowBoundaryError(
                "ksea membership inferential decision status must be "
                "ordinary_p_q_unavailable when ordinary_p_q_available is false"
            )
        object.__setattr__(
            self,
            "ordinary_p_q_available",
            ordinary_p_q_available,
        )
        object.__setattr__(self, "status", status)
        object.__setattr__(
            self,
            "reason",
            _require_non_empty_text(
                self.reason,
                field_name="ksea_membership_inferential_decision.reason",
            ),
        )
        object.__setattr__(
            self,
            "missing_evidence",
            _string_tuple(
                self.missing_evidence,
                field_name="ksea_membership_inferential_decision.missing_evidence",
            ),
        )
        object.__setattr__(
            self,
            "policy_version",
            _require_non_empty_text(
                self.policy_version,
                field_name="ksea_membership_inferential_decision.policy_version",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible decision payload."""

        return {
            "ordinary_p_q_available": bool(self.ordinary_p_q_available),
            "status": self.status,
            "reason": self.reason,
            "missing_evidence": list(self.missing_evidence),
            "policy_version": self.policy_version,
        }


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
    selection_quantitative_matrix_fingerprint: TableFingerprint | None = None
    tested_quantitative_matrix_fingerprint: TableFingerprint | None = None
    # Legacy compatibility alias for pre-v2 payloads and constructor calls. This
    # is treated only as the selection quantitative matrix fingerprint, never as
    # KSEA tested-matrix evidence.
    quantitative_dataset_fingerprint: TableFingerprint | None = None
    consumed_tested_matrix: bool = False
    selected_kinase_universe: tuple[str, ...] = ()
    selected_substrate_universe: tuple[str, ...] = ()
    inferential_eligible: bool | None = None
    inferential_eligibility_reason: str | None = None
    _inferential_decision: KseaMembershipInferentialDecision = field(
        init=False,
        repr=False,
    )

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
        source_reference_fingerprints = _coerce_table_fingerprints(
            self.source_reference_fingerprints,
            field_name="activity_membership_selection.source_reference_fingerprints",
        )
        legacy_quantitative_fingerprint = _coerce_optional_table_fingerprint(
            self.quantitative_dataset_fingerprint,
            field_name="activity_membership_selection.quantitative_dataset_fingerprint",
        )
        selection_quantitative_fingerprint = _coerce_optional_table_fingerprint(
            self.selection_quantitative_matrix_fingerprint,
            field_name=(
                "activity_membership_selection."
                "selection_quantitative_matrix_fingerprint"
            ),
        )
        if (
            legacy_quantitative_fingerprint is not None
            and selection_quantitative_fingerprint is not None
            and legacy_quantitative_fingerprint != selection_quantitative_fingerprint
        ):
            raise WorkflowBoundaryError(
                "activity_membership_selection.quantitative_dataset_fingerprint "
                "is a legacy alias for selection_quantitative_matrix_fingerprint "
                "and must match it when both are supplied"
            )
        if selection_quantitative_fingerprint is None:
            selection_quantitative_fingerprint = legacy_quantitative_fingerprint
        tested_quantitative_fingerprint = _coerce_optional_table_fingerprint(
            self.tested_quantitative_matrix_fingerprint,
            field_name=(
                "activity_membership_selection.tested_quantitative_matrix_fingerprint"
            ),
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
            "selection_quantitative_matrix_fingerprint",
            selection_quantitative_fingerprint,
        )
        object.__setattr__(
            self,
            "tested_quantitative_matrix_fingerprint",
            tested_quantitative_fingerprint,
        )
        object.__setattr__(
            self,
            "quantitative_dataset_fingerprint",
            selection_quantitative_fingerprint,
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
        decision = derive_ksea_membership_inferential_decision(self)
        _validate_requested_inferential_fields(
            requested_eligible=self.inferential_eligible,
            requested_reason=self.inferential_eligibility_reason,
            decision=decision,
        )
        object.__setattr__(self, "_inferential_decision", decision)
        object.__setattr__(
            self,
            "inferential_eligible",
            bool(decision.ordinary_p_q_available),
        )
        object.__setattr__(
            self,
            "inferential_eligibility_reason",
            decision.reason,
        )

    @property
    def inferential_decision(self) -> KseaMembershipInferentialDecision:
        """Return the science-derived KSEA ordinary p/q availability decision."""

        return self._inferential_decision

    @property
    def inferential_status(self) -> str:
        """Return a compact status token for p/q-value availability."""

        return self._inferential_decision.status

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
            "selection_quantitative_matrix_fingerprint": (
                None
                if self.selection_quantitative_matrix_fingerprint is None
                else table_fingerprint_to_payload(
                    self.selection_quantitative_matrix_fingerprint
                )
            ),
            "tested_quantitative_matrix_fingerprint": (
                None
                if self.tested_quantitative_matrix_fingerprint is None
                else table_fingerprint_to_payload(
                    self.tested_quantitative_matrix_fingerprint
                )
            ),
            "consumed_tested_matrix": bool(self.consumed_tested_matrix),
            "selected_kinase_universe": list(self.selected_kinase_universe),
            "selected_substrate_universe": list(self.selected_substrate_universe),
            "inferential_eligible": bool(self.inferential_eligible),
            "inferential_status": self.inferential_status,
            "inferential_eligibility_reason": self.inferential_eligibility_reason,
            "inferential_decision": self.inferential_decision.to_payload(),
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
        legacy_quantitative_fingerprint = _optional_fingerprint_from_payload(
            payload.get("quantitative_dataset_fingerprint"),
            field_name="activity_membership_selection.quantitative_dataset_fingerprint",
        )
        selection_quantitative_fingerprint = _optional_fingerprint_from_payload(
            payload.get("selection_quantitative_matrix_fingerprint"),
            field_name=(
                "activity_membership_selection."
                "selection_quantitative_matrix_fingerprint"
            ),
        )
        if selection_quantitative_fingerprint is None:
            selection_quantitative_fingerprint = legacy_quantitative_fingerprint
        tested_quantitative_fingerprint = _optional_fingerprint_from_payload(
            payload.get("tested_quantitative_matrix_fingerprint"),
            field_name=(
                "activity_membership_selection.tested_quantitative_matrix_fingerprint"
            ),
        )
        selection = cls(
            source_category=str(payload.get("source_category", "")),
            selection_method=str(payload.get("selection_method", "")),
            selection_method_version=str(payload.get("selection_method_version", "")),
            score_source=str(payload.get("score_source", "")),
            threshold_top_k_policy=_require_mapping(
                payload.get("threshold_top_k_policy", {}),
                field_name="activity_membership_selection.threshold_top_k_policy",
            ),
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=(
                selection_quantitative_fingerprint
            ),
            tested_quantitative_matrix_fingerprint=tested_quantitative_fingerprint,
            quantitative_dataset_fingerprint=legacy_quantitative_fingerprint,
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
            inferential_eligible=_optional_bool(
                payload.get("inferential_eligible"),
                field_name="activity_membership_selection.inferential_eligible",
            ),
            inferential_eligibility_reason=_optional_text(
                payload.get("inferential_eligibility_reason"),
                field_name=(
                    "activity_membership_selection.inferential_eligibility_reason"
                ),
            ),
        )
        _validate_serialized_inferential_fields(payload, selection=selection)
        return selection

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
            selection_quantitative_matrix_fingerprint=None,
            tested_quantitative_matrix_fingerprint=None,
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
        )


def fingerprint_ksea_tested_quantitative_matrix(
    matrix: pd.DataFrame,
) -> TableFingerprint:
    """Fingerprint the exact KSEA quantitative/background matrix under policy v2."""

    return fingerprint_table_normalized_axes(
        matrix,
        name=KSEA_TESTED_QUANTITATIVE_MATRIX_FINGERPRINT_NAME,
    )


def fingerprint_ksea_selection_quantitative_matrix(
    matrix: pd.DataFrame,
) -> TableFingerprint:
    """Fingerprint the quantitative matrix consumed by membership selection."""

    return fingerprint_table_normalized_axes(
        matrix,
        name=KSEA_SELECTION_QUANTITATIVE_MATRIX_FINGERPRINT_NAME,
    )


def derive_ksea_membership_inferential_decision(
    selection: ActivityMembershipSelection,
) -> KseaMembershipInferentialDecision:
    """Derive ordinary KSEA p/q availability from typed membership facts."""

    source_category = selection.source_category
    missing_evidence = _base_missing_evidence(selection)
    if source_category == ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN:
        return _unavailable_decision(
            KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON,
            missing_evidence=("source_category", *missing_evidence),
        )
    if source_category == ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED:
        return _unavailable_decision(
            (
                KSEA_MEMBERSHIP_CONSUMED_TESTED_MATRIX_REASON
                if selection.consumed_tested_matrix
                else KSEA_MEMBERSHIP_PROFILE_DERIVED_REASON
            ),
            missing_evidence=missing_evidence,
        )
    if source_category == ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF:
        return _unavailable_decision(
            (
                KSEA_MEMBERSHIP_CONSUMED_TESTED_MATRIX_REASON
                if selection.consumed_tested_matrix
                else KSEA_MEMBERSHIP_FUSED_PROFILE_MOTIF_REASON
            ),
            missing_evidence=missing_evidence,
        )
    if source_category == ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED:
        if (
            selection.consumed_tested_matrix
            or selection.selection_quantitative_matrix_fingerprint is not None
            or _policy_bool(
                selection.threshold_top_k_policy, "data_adaptive_membership"
            )
        ):
            return _unavailable_decision(
                (
                    KSEA_MEMBERSHIP_CONSUMED_TESTED_MATRIX_REASON
                    if selection.consumed_tested_matrix
                    else KSEA_MEMBERSHIP_QUANTITATIVE_SELECTION_REASON
                ),
                missing_evidence=missing_evidence,
            )
        return _unavailable_decision(
            KSEA_MEMBERSHIP_PREDICTION_SELECTED_REASON,
            missing_evidence=missing_evidence,
        )
    if source_category == ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE:
        return _independent_reference_decision(
            selection,
            expected_policy=KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_FIXED_EXTERNAL_REFERENCE,
            extra_missing_evidence=missing_evidence,
        )
    if source_category == ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF:
        extra_missing = list(missing_evidence)
        if selection.score_source != _KSEA_SEQUENCE_ONLY_MOTIF_SCORE_SOURCE:
            extra_missing.append("score_source.kinase_library_motif_scores")
        return _independent_reference_decision(
            selection,
            expected_policy=KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_SEQUENCE_ONLY_MOTIF,
            extra_missing_evidence=tuple(extra_missing),
        )
    return _unavailable_decision(
        KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON,
        missing_evidence=("source_category", *missing_evidence),
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


def _coerce_optional_table_fingerprint(
    value: object,
    *,
    field_name: str,
) -> TableFingerprint | None:
    if value is None:
        return None
    if not isinstance(value, TableFingerprint):
        raise WorkflowBoundaryError(f"{field_name} must be TableFingerprint or None")
    return value


def _optional_fingerprint_from_payload(
    payload: object,
    *,
    field_name: str,
) -> TableFingerprint | None:
    if payload is None:
        return None
    return table_fingerprint_from_payload(
        _require_mapping(payload, field_name=field_name)
    )


def _validate_requested_inferential_fields(
    *,
    requested_eligible: object,
    requested_reason: object,
    decision: KseaMembershipInferentialDecision,
) -> None:
    if requested_eligible is not None:
        requested = _require_bool(
            requested_eligible,
            field_name="activity_membership_selection.inferential_eligible",
        )
        if requested is not bool(decision.ordinary_p_q_available):
            raise WorkflowBoundaryError(
                "activity_membership_selection.inferential_eligible is "
                "science-derived and cannot contradict the KSEA membership "
                "inferential policy; "
                f"requested={requested!r}, derived="
                f"{decision.ordinary_p_q_available!r}, reason={decision.reason!r}"
            )
    if requested_reason is not None:
        requested = _require_non_empty_text(
            requested_reason,
            field_name="activity_membership_selection.inferential_eligibility_reason",
        )
        if requested != decision.reason:
            raise WorkflowBoundaryError(
                "activity_membership_selection.inferential_eligibility_reason "
                "is science-derived and cannot contradict the KSEA membership "
                "inferential policy; "
                f"requested={requested!r}, derived={decision.reason!r}"
            )


def _validate_serialized_inferential_fields(
    payload: Mapping[str, object],
    *,
    selection: ActivityMembershipSelection,
) -> None:
    serialized_status = payload.get("inferential_status")
    if serialized_status is not None:
        status = _require_non_empty_text(
            serialized_status,
            field_name="activity_membership_selection.inferential_status",
        )
        if status != selection.inferential_status:
            raise WorkflowBoundaryError(
                "activity_membership_selection.inferential_status is "
                "science-derived and cannot contradict the reconstructed "
                "membership facts; "
                f"serialized={status!r}, derived={selection.inferential_status!r}"
            )
    raw_decision = payload.get("inferential_decision")
    if raw_decision is None:
        return
    decision_payload = _require_mapping(
        raw_decision,
        field_name="activity_membership_selection.inferential_decision",
    )
    _compare_optional_decision_bool(
        decision_payload,
        key="ordinary_p_q_available",
        expected=selection.inferential_decision.ordinary_p_q_available,
    )
    _compare_optional_decision_text(
        decision_payload,
        key="status",
        expected=selection.inferential_decision.status,
    )
    _compare_optional_decision_text(
        decision_payload,
        key="reason",
        expected=selection.inferential_decision.reason,
    )
    if "missing_evidence" in decision_payload:
        serialized_missing = _string_tuple(
            decision_payload.get("missing_evidence", ()),
            field_name=(
                "activity_membership_selection.inferential_decision.missing_evidence"
            ),
        )
        if serialized_missing != selection.inferential_decision.missing_evidence:
            raise WorkflowBoundaryError(
                "activity_membership_selection.inferential_decision."
                "missing_evidence cannot contradict the reconstructed "
                "membership facts; "
                f"serialized={serialized_missing!r}, derived="
                f"{selection.inferential_decision.missing_evidence!r}"
            )


def _compare_optional_decision_bool(
    payload: Mapping[str, object],
    *,
    key: str,
    expected: bool,
) -> None:
    if key not in payload:
        return
    observed = _require_bool(
        payload.get(key),
        field_name=f"activity_membership_selection.inferential_decision.{key}",
    )
    if observed is expected:
        return
    raise WorkflowBoundaryError(
        "activity_membership_selection.inferential_decision cannot contradict "
        f"the reconstructed membership facts for {key}; "
        f"serialized={observed!r}, derived={expected!r}"
    )


def _compare_optional_decision_text(
    payload: Mapping[str, object],
    *,
    key: str,
    expected: str,
) -> None:
    if key not in payload:
        return
    observed = _require_non_empty_text(
        payload.get(key),
        field_name=f"activity_membership_selection.inferential_decision.{key}",
    )
    if observed == expected:
        return
    raise WorkflowBoundaryError(
        "activity_membership_selection.inferential_decision cannot contradict "
        f"the reconstructed membership facts for {key}; "
        f"serialized={observed!r}, derived={expected!r}"
    )


def _base_missing_evidence(
    selection: ActivityMembershipSelection,
) -> tuple[str, ...]:
    missing: list[str] = []
    if selection.consumed_tested_matrix and (
        selection.selection_quantitative_matrix_fingerprint is None
    ):
        missing.append("selection_quantitative_matrix_fingerprint")
    return tuple(missing)


def _independent_reference_decision(
    selection: ActivityMembershipSelection,
    *,
    expected_policy: str,
    extra_missing_evidence: tuple[str, ...],
) -> KseaMembershipInferentialDecision:
    missing = list(extra_missing_evidence)
    if selection.consumed_tested_matrix:
        missing.append("consumed_tested_matrix")
    if selection.selection_quantitative_matrix_fingerprint is not None:
        missing.append("selection_quantitative_matrix_fingerprint")
    if selection.tested_quantitative_matrix_fingerprint is None:
        missing.append("tested_quantitative_matrix_fingerprint")
    if not selection.source_reference_fingerprints:
        missing.append("source_reference_fingerprints")
    if not selection.selected_kinase_universe:
        missing.append("selected_kinase_universe")
    if not selection.selected_substrate_universe:
        missing.append("selected_substrate_universe")
    independence_policy = str(
        selection.threshold_top_k_policy.get("independent_membership_policy", "")
    )
    if independence_policy != expected_policy:
        missing.append("threshold_top_k_policy.independent_membership_policy")
    independence_policy_version = str(
        selection.threshold_top_k_policy.get(
            "independent_membership_policy_version",
            "",
        )
    )
    if independence_policy_version != KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION:
        missing.append("threshold_top_k_policy.independent_membership_policy_version")
    if missing:
        return _unavailable_decision(
            KSEA_MEMBERSHIP_INCOMPLETE_INDEPENDENCE_EVIDENCE_REASON,
            missing_evidence=tuple(dict.fromkeys(missing)),
        )
    return KseaMembershipInferentialDecision(
        ordinary_p_q_available=True,
        status=KSEA_MEMBERSHIP_INFERENCE_STATUS_ELIGIBLE,
        reason=KSEA_MEMBERSHIP_ELIGIBLE_REASON,
        missing_evidence=(),
        policy_version=KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION,
    )


def _unavailable_decision(
    reason: str,
    *,
    missing_evidence: tuple[str, ...],
) -> KseaMembershipInferentialDecision:
    return KseaMembershipInferentialDecision(
        ordinary_p_q_available=False,
        status=KSEA_MEMBERSHIP_INFERENCE_STATUS_UNAVAILABLE,
        reason=reason,
        missing_evidence=tuple(dict.fromkeys(missing_evidence)),
        policy_version=KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION,
    )


def _policy_bool(policy: Mapping[str, object], key: str) -> bool:
    value = policy.get(key)
    return value if isinstance(value, bool) else False


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


def _optional_bool(value: object, *, field_name: str) -> bool | None:
    if value is None:
        return None
    return _require_bool(value, field_name=field_name)


def _optional_text(value: object, *, field_name: str) -> str | None:
    if value is None:
        return None
    return _require_non_empty_text(value, field_name=field_name)


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
    "derive_ksea_membership_inferential_decision",
    "fingerprint_ksea_selection_quantitative_matrix",
    "fingerprint_ksea_tested_quantitative_matrix",
    "KSEA_MEMBERSHIP_FUSED_PROFILE_MOTIF_REASON",
    "KSEA_MEMBERSHIP_INCOMPLETE_INDEPENDENCE_EVIDENCE_REASON",
    "KSEA_MEMBERSHIP_CONSUMED_TESTED_MATRIX_REASON",
    "KSEA_MEMBERSHIP_ELIGIBLE_REASON",
    "KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_FIXED_EXTERNAL_REFERENCE",
    "KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_SEQUENCE_ONLY_MOTIF",
    "KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION",
    "KSEA_MEMBERSHIP_INFERENCE_STATUS_ELIGIBLE",
    "KSEA_MEMBERSHIP_INFERENCE_STATUS_UNAVAILABLE",
    "KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION",
    "KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON",
    "KSEA_MEMBERSHIP_PREDICTION_SELECTED_REASON",
    "KSEA_MEMBERSHIP_PROFILE_DERIVED_REASON",
    "KSEA_MEMBERSHIP_QUANTITATIVE_SELECTION_REASON",
    "KSEA_SELECTION_QUANTITATIVE_MATRIX_FINGERPRINT_NAME",
    "KSEA_TESTED_QUANTITATIVE_MATRIX_FINGERPRINT_NAME",
    "KseaMembershipInferentialDecision",
    "selected_substrate_universe_from_prediction_matrix",
]
