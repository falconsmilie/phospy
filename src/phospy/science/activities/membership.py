"""Typed substrate-membership provenance for activity inference."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Final, TypeVar

import numpy as np
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.policies import PolicyEnum, coerce_policy_enum
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
ACTIVITY_MEMBERSHIP_SOURCE_INCOMPLETE: Final = "incomplete"

ACTIVITY_MEMBERSHIP_SOURCE_CATEGORIES: Final = frozenset(
    {
        ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE,
        ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF,
        ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED,
        ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF,
        ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED,
        ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN,
        ACTIVITY_MEMBERSHIP_SOURCE_INCOMPLETE,
    }
)

ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION: Final = "4"
ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION: Final = "2"
KSEA_MEMBERSHIP_INFERENTIAL_POLICY_VERSION: Final = "4"
KSEA_TESTED_QUANTITATIVE_MATRIX_FINGERPRINT_NAME: Final = (
    "dataset.ksea_background_phospho_matrix"
)
KSEA_SELECTION_QUANTITATIVE_MATRIX_FINGERPRINT_NAME: Final = (
    "dataset.scoring_phospho_matrix"
)
KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION: Final = "2"
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
_KSEA_INDEPENDENCE_POLICY_FIELD: Final = "independent_membership_policy"
_KSEA_INDEPENDENCE_POLICY_VERSION_FIELD: Final = "independent_membership_policy_version"
_KSEA_DATA_ADAPTIVE_MEMBERSHIP_FIELD: Final = "data_adaptive_membership"
_PolicyEnumT = TypeVar("_PolicyEnumT", bound=PolicyEnum)


class ActivityMembershipSelectionProcessKind(PolicyEnum):
    """Closed scientific selection-process classification for KSEA membership."""

    FIXED_EXTERNAL_REFERENCE = ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE
    SEQUENCE_ONLY_MOTIF = ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF
    PROFILE_DERIVED = ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED
    FUSED_PROFILE_MOTIF = ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF
    PREDICTION_SELECTED = ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED
    UNKNOWN = ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN
    INCOMPLETE = ACTIVITY_MEMBERSHIP_SOURCE_INCOMPLETE


class ActivityMembershipScoreSourceKind(PolicyEnum):
    """Closed scientific score-source kind used by membership selection."""

    EXTERNAL_REFERENCE = "external_reference"
    KINASE_LIBRARY_MOTIF = "kinase_library_motif"
    PROFILE_DERIVED = "profile_derived"
    FUSED_PROFILE_MOTIF = "fused_profile_motif"
    PREDICTION_DERIVED = "prediction_derived"
    UNKNOWN = "unknown"


class ActivityMembershipIndependencePolicyKind(PolicyEnum):
    """Closed independence-policy evidence accepted by the activity domain."""

    FIXED_EXTERNAL_REFERENCE = (
        KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_FIXED_EXTERNAL_REFERENCE
    )
    SEQUENCE_ONLY_MOTIF = KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_SEQUENCE_ONLY_MOTIF


_PROCESS_KIND_SOURCE_CATEGORY: Final = {
    ActivityMembershipSelectionProcessKind.FIXED_EXTERNAL_REFERENCE: (
        ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE
    ),
    ActivityMembershipSelectionProcessKind.SEQUENCE_ONLY_MOTIF: (
        ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF
    ),
    ActivityMembershipSelectionProcessKind.PROFILE_DERIVED: (
        ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED
    ),
    ActivityMembershipSelectionProcessKind.FUSED_PROFILE_MOTIF: (
        ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF
    ),
    ActivityMembershipSelectionProcessKind.PREDICTION_SELECTED: (
        ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED
    ),
    ActivityMembershipSelectionProcessKind.UNKNOWN: ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN,
    ActivityMembershipSelectionProcessKind.INCOMPLETE: (
        ACTIVITY_MEMBERSHIP_SOURCE_INCOMPLETE
    ),
}
_KSEA_DECISION_POLICY_MAPPING_FIELDS: Final = frozenset(
    {
        _KSEA_DATA_ADAPTIVE_MEMBERSHIP_FIELD,
        _KSEA_INDEPENDENCE_POLICY_FIELD,
        _KSEA_INDEPENDENCE_POLICY_VERSION_FIELD,
    }
)


@dataclass(frozen=True, slots=True)
class ActivityMembershipIndependenceEvidence:
    """Typed source-specific evidence that membership is independent of KSEA data."""

    policy_kind: ActivityMembershipIndependencePolicyKind
    policy_version: str

    def __post_init__(self) -> None:
        policy_kind = _coerce_policy_enum(
            ActivityMembershipIndependencePolicyKind,
            self.policy_kind,
            field_name="activity_membership_selection.selection_evidence."
            "independence_evidence.policy_kind",
        )
        object.__setattr__(self, "policy_kind", policy_kind)
        object.__setattr__(
            self,
            "policy_version",
            _require_non_empty_text(
                self.policy_version,
                field_name="activity_membership_selection.selection_evidence."
                "independence_evidence.policy_version",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible independence-evidence payload."""

        return {
            "policy_kind": self.policy_kind.value,
            "policy_version": self.policy_version,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMembershipIndependenceEvidence:
        """Reconstruct typed independence evidence from a payload."""

        return cls(
            policy_kind=_coerce_policy_enum(
                ActivityMembershipIndependencePolicyKind,
                payload.get("policy_kind", ""),
                field_name="activity_membership_selection.selection_evidence."
                "independence_evidence.policy_kind",
            ),
            policy_version=str(payload.get("policy_version", "")),
        )


@dataclass(frozen=True, slots=True)
class ActivityMembershipSelectionEvidence:
    """Closed decision-bearing evidence for activity membership provenance."""

    selection_process_kind: ActivityMembershipSelectionProcessKind
    selection_contract_version: str
    score_source_kind: ActivityMembershipScoreSourceKind
    data_adaptive_membership: bool | None
    consumed_tested_matrix: bool
    independence_evidence: ActivityMembershipIndependenceEvidence | None = None

    def __post_init__(self) -> None:
        process_kind = _coerce_policy_enum(
            ActivityMembershipSelectionProcessKind,
            self.selection_process_kind,
            field_name="activity_membership_selection.selection_evidence."
            "selection_process_kind",
        )
        score_source_kind = _coerce_policy_enum(
            ActivityMembershipScoreSourceKind,
            self.score_source_kind,
            field_name="activity_membership_selection.selection_evidence."
            "score_source_kind",
        )
        data_adaptive_membership = (
            None
            if self.data_adaptive_membership is None
            else _require_bool(
                self.data_adaptive_membership,
                field_name="activity_membership_selection.selection_evidence."
                "data_adaptive_membership",
            )
        )
        consumed_tested_matrix = _require_bool(
            self.consumed_tested_matrix,
            field_name="activity_membership_selection.selection_evidence."
            "consumed_tested_matrix",
        )
        independence_evidence = _coerce_independence_evidence(
            self.independence_evidence
        )
        object.__setattr__(self, "selection_process_kind", process_kind)
        object.__setattr__(
            self,
            "selection_contract_version",
            _require_non_empty_text(
                self.selection_contract_version,
                field_name="activity_membership_selection.selection_evidence."
                "selection_contract_version",
            ),
        )
        object.__setattr__(self, "score_source_kind", score_source_kind)
        object.__setattr__(
            self,
            "data_adaptive_membership",
            data_adaptive_membership,
        )
        object.__setattr__(self, "consumed_tested_matrix", consumed_tested_matrix)
        object.__setattr__(self, "independence_evidence", independence_evidence)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible selection-evidence payload."""

        return {
            "selection_process_kind": self.selection_process_kind.value,
            "selection_contract_version": self.selection_contract_version,
            "score_source_kind": self.score_source_kind.value,
            "data_adaptive_membership": self.data_adaptive_membership,
            "consumed_tested_matrix": bool(self.consumed_tested_matrix),
            "independence_evidence": (
                None
                if self.independence_evidence is None
                else self.independence_evidence.to_payload()
            ),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityMembershipSelectionEvidence:
        """Reconstruct typed selection evidence from a payload."""

        return cls(
            selection_process_kind=_coerce_policy_enum(
                ActivityMembershipSelectionProcessKind,
                payload.get("selection_process_kind", ""),
                field_name="activity_membership_selection.selection_evidence."
                "selection_process_kind",
            ),
            selection_contract_version=str(
                payload.get("selection_contract_version", "")
            ),
            score_source_kind=_coerce_policy_enum(
                ActivityMembershipScoreSourceKind,
                payload.get("score_source_kind", ""),
                field_name="activity_membership_selection.selection_evidence."
                "score_source_kind",
            ),
            data_adaptive_membership=(
                None
                if payload.get("data_adaptive_membership") is None
                else _require_bool(
                    payload.get("data_adaptive_membership"),
                    field_name="activity_membership_selection.selection_evidence."
                    "data_adaptive_membership",
                )
            ),
            consumed_tested_matrix=_require_bool(
                payload.get("consumed_tested_matrix"),
                field_name="activity_membership_selection.selection_evidence."
                "consumed_tested_matrix",
            ),
            independence_evidence=(
                None
                if payload.get("independence_evidence") is None
                else ActivityMembershipIndependenceEvidence.from_payload(
                    _require_mapping(
                        payload.get("independence_evidence"),
                        field_name="activity_membership_selection.selection_evidence."
                        "independence_evidence",
                    )
                )
            ),
        )


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
    membership_selection_schema_version: str = (
        ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
    )
    threshold_top_k_policy: Mapping[str, object] = field(default_factory=dict)
    selection_evidence: ActivityMembershipSelectionEvidence | None = None
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
        requested_source_category = str(self.source_category).strip()
        if requested_source_category not in ACTIVITY_MEMBERSHIP_SOURCE_CATEGORIES:
            allowed = ", ".join(sorted(ACTIVITY_MEMBERSHIP_SOURCE_CATEGORIES))
            raise WorkflowBoundaryError(
                "activity membership source_category must be one of: "
                f"{allowed}; got {requested_source_category!r}"
            )
        schema_version = _require_supported_membership_payload_schema_version(
            self.membership_selection_schema_version
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
        threshold_top_k_policy = _freeze_policy_mapping(self.threshold_top_k_policy)
        selection_evidence = _coerce_selection_evidence(
            self.selection_evidence,
            consumed_tested_matrix=consumed_tested_matrix,
        )
        derived_source_category = _source_category_from_selection_evidence(
            selection_evidence
        )
        if self.selection_evidence is not None and (
            requested_source_category != derived_source_category
        ):
            raise WorkflowBoundaryError(
                "activity_membership_selection.source_category is derived from "
                "selection_evidence and cannot contradict it; "
                f"serialized={requested_source_category!r}, derived="
                f"{derived_source_category!r}"
            )
        if selection_evidence.consumed_tested_matrix is not consumed_tested_matrix:
            raise WorkflowBoundaryError(
                "activity_membership_selection.consumed_tested_matrix must match "
                "selection_evidence.consumed_tested_matrix"
            )
        object.__setattr__(self, "source_category", derived_source_category)
        object.__setattr__(
            self,
            "membership_selection_schema_version",
            schema_version,
        )
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
            threshold_top_k_policy,
        )
        object.__setattr__(
            self,
            "selection_evidence",
            selection_evidence,
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
        _validate_membership_source_fact_coherence(self)
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

    @property
    def selection_process_kind(self) -> ActivityMembershipSelectionProcessKind:
        """Return the closed scientific membership-selection process kind."""

        return _selection_evidence_model(self.selection_evidence).selection_process_kind

    @property
    def score_source_kind(self) -> ActivityMembershipScoreSourceKind:
        """Return the closed scientific score-source kind."""

        return _selection_evidence_model(self.selection_evidence).score_source_kind

    @property
    def data_adaptive_membership(self) -> bool | None:
        """Return whether the membership selection was explicitly data-adaptive."""

        return _selection_evidence_model(
            self.selection_evidence
        ).data_adaptive_membership

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible membership-selection payload."""

        source_reference_fingerprints = [
            table_fingerprint_to_payload(item)
            for item in self.source_reference_fingerprints
        ]
        return {
            "membership_selection_schema_version": (
                self.membership_selection_schema_version
            ),
            "source_category": self.source_category,
            "selection_method": self.selection_method,
            "selection_method_version": self.selection_method_version,
            "score_source": self.score_source,
            "selection_evidence": _selection_evidence_model(
                self.selection_evidence
            ).to_payload(),
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
        raw_schema_version = payload.get("membership_selection_schema_version")
        if raw_schema_version is None:
            if "selection_evidence" in payload:
                raise WorkflowBoundaryError(
                    "activity_membership_selection.membership_selection_schema_version "
                    "is required when selection_evidence is present"
                )
            schema_version = ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
            selection_evidence = None
            consumed_tested_matrix = _require_bool(
                payload.get("consumed_tested_matrix", False),
                field_name=("activity_membership_selection.consumed_tested_matrix"),
            )
        else:
            schema_version = _require_supported_membership_payload_schema_version(
                raw_schema_version
            )
            raw_selection_evidence = payload.get("selection_evidence")
            selection_evidence = ActivityMembershipSelectionEvidence.from_payload(
                _require_mapping(
                    raw_selection_evidence,
                    field_name="activity_membership_selection.selection_evidence",
                )
            )
            consumed_tested_matrix = _require_bool(
                payload.get("consumed_tested_matrix"),
                field_name=("activity_membership_selection.consumed_tested_matrix"),
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
            membership_selection_schema_version=schema_version,
            threshold_top_k_policy=_require_mapping(
                payload.get("threshold_top_k_policy", {}),
                field_name="activity_membership_selection.threshold_top_k_policy",
            ),
            selection_evidence=selection_evidence,
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=(
                selection_quantitative_fingerprint
            ),
            tested_quantitative_matrix_fingerprint=tested_quantitative_fingerprint,
            quantitative_dataset_fingerprint=legacy_quantitative_fingerprint,
            consumed_tested_matrix=consumed_tested_matrix,
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
            membership_selection_schema_version=(
                ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
            ),
            threshold_top_k_policy={},
            selection_evidence=ActivityMembershipSelectionEvidence(
                selection_process_kind=ActivityMembershipSelectionProcessKind.UNKNOWN,
                selection_contract_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
                score_source_kind=ActivityMembershipScoreSourceKind.UNKNOWN,
                data_adaptive_membership=None,
                consumed_tested_matrix=False,
                independence_evidence=None,
            ),
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

    @classmethod
    def fixed_external_reference(
        cls,
        *,
        provider_method_identifier: str,
        provider_method_version: str,
        provider_score_source_identifier: str,
        threshold_top_k_policy: Mapping[str, object] | None = None,
        source_reference_fingerprints: tuple[TableFingerprint, ...] = (),
        tested_quantitative_matrix_fingerprint: TableFingerprint | None = None,
        selected_kinase_universe: Iterable[object] = (),
        selected_substrate_universe: Iterable[object] = (),
    ) -> ActivityMembershipSelection:
        """Build current fixed-external membership evidence.

        Provider labels are retained as descriptive provenance. The scientific
        classification comes only from the closed evidence object.
        """

        return cls(
            source_category=ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE,
            selection_method=provider_method_identifier,
            selection_method_version=provider_method_version,
            score_source=provider_score_source_identifier,
            membership_selection_schema_version=(
                ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
            ),
            threshold_top_k_policy=(
                {} if threshold_top_k_policy is None else threshold_top_k_policy
            ),
            selection_evidence=ActivityMembershipSelectionEvidence(
                selection_process_kind=(
                    ActivityMembershipSelectionProcessKind.FIXED_EXTERNAL_REFERENCE
                ),
                selection_contract_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
                score_source_kind=ActivityMembershipScoreSourceKind.EXTERNAL_REFERENCE,
                data_adaptive_membership=False,
                consumed_tested_matrix=False,
                independence_evidence=ActivityMembershipIndependenceEvidence(
                    policy_kind=(
                        ActivityMembershipIndependencePolicyKind.FIXED_EXTERNAL_REFERENCE
                    ),
                    policy_version=KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION,
                ),
            ),
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=None,
            tested_quantitative_matrix_fingerprint=(
                tested_quantitative_matrix_fingerprint
            ),
            quantitative_dataset_fingerprint=None,
            consumed_tested_matrix=False,
            selected_kinase_universe=_string_tuple(
                selected_kinase_universe,
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
            selected_substrate_universe=_string_tuple(
                selected_substrate_universe,
                field_name="activity_membership_selection.selected_substrate_universe",
            ),
        )

    @classmethod
    def sequence_only_motif(
        cls,
        *,
        provider_method_identifier: str,
        provider_method_version: str,
        threshold_top_k_policy: Mapping[str, object] | None = None,
        source_reference_fingerprints: tuple[TableFingerprint, ...] = (),
        tested_quantitative_matrix_fingerprint: TableFingerprint | None = None,
        selected_kinase_universe: Iterable[object] = (),
        selected_substrate_universe: Iterable[object] = (),
    ) -> ActivityMembershipSelection:
        """Build current sequence-only motif membership evidence."""

        return cls(
            source_category=ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF,
            selection_method=provider_method_identifier,
            selection_method_version=provider_method_version,
            score_source=_KSEA_SEQUENCE_ONLY_MOTIF_SCORE_SOURCE,
            membership_selection_schema_version=(
                ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
            ),
            threshold_top_k_policy=(
                {} if threshold_top_k_policy is None else threshold_top_k_policy
            ),
            selection_evidence=ActivityMembershipSelectionEvidence(
                selection_process_kind=(
                    ActivityMembershipSelectionProcessKind.SEQUENCE_ONLY_MOTIF
                ),
                selection_contract_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
                score_source_kind=(
                    ActivityMembershipScoreSourceKind.KINASE_LIBRARY_MOTIF
                ),
                data_adaptive_membership=False,
                consumed_tested_matrix=False,
                independence_evidence=ActivityMembershipIndependenceEvidence(
                    policy_kind=(
                        ActivityMembershipIndependencePolicyKind.SEQUENCE_ONLY_MOTIF
                    ),
                    policy_version=KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION,
                ),
            ),
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=None,
            tested_quantitative_matrix_fingerprint=(
                tested_quantitative_matrix_fingerprint
            ),
            quantitative_dataset_fingerprint=None,
            consumed_tested_matrix=False,
            selected_kinase_universe=_string_tuple(
                selected_kinase_universe,
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
            selected_substrate_universe=_string_tuple(
                selected_substrate_universe,
                field_name="activity_membership_selection.selected_substrate_universe",
            ),
        )

    @classmethod
    def profile_derived(
        cls,
        *,
        selection_method: str,
        score_source: str,
        threshold_top_k_policy: Mapping[str, object] | None = None,
        source_reference_fingerprints: tuple[TableFingerprint, ...] = (),
        selection_quantitative_matrix_fingerprint: TableFingerprint | None = None,
        tested_quantitative_matrix_fingerprint: TableFingerprint | None = None,
        consumed_tested_matrix: bool,
        selected_kinase_universe: Iterable[object] = (),
        selected_substrate_universe: Iterable[object] = (),
    ) -> ActivityMembershipSelection:
        """Build current profile-derived membership evidence."""

        return cls._adaptive(
            source_category=ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED,
            process_kind=ActivityMembershipSelectionProcessKind.PROFILE_DERIVED,
            score_source_kind=ActivityMembershipScoreSourceKind.PROFILE_DERIVED,
            selection_method=selection_method,
            score_source=score_source,
            threshold_top_k_policy=threshold_top_k_policy,
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=(
                selection_quantitative_matrix_fingerprint
            ),
            tested_quantitative_matrix_fingerprint=(
                tested_quantitative_matrix_fingerprint
            ),
            consumed_tested_matrix=consumed_tested_matrix,
            selected_kinase_universe=_string_tuple(
                selected_kinase_universe,
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
            selected_substrate_universe=_string_tuple(
                selected_substrate_universe,
                field_name="activity_membership_selection.selected_substrate_universe",
            ),
        )

    @classmethod
    def fused_profile_motif(
        cls,
        *,
        selection_method: str,
        score_source: str,
        threshold_top_k_policy: Mapping[str, object] | None = None,
        source_reference_fingerprints: tuple[TableFingerprint, ...] = (),
        selection_quantitative_matrix_fingerprint: TableFingerprint | None = None,
        tested_quantitative_matrix_fingerprint: TableFingerprint | None = None,
        consumed_tested_matrix: bool,
        selected_kinase_universe: Iterable[object] = (),
        selected_substrate_universe: Iterable[object] = (),
    ) -> ActivityMembershipSelection:
        """Build current fused profile/motif membership evidence."""

        return cls._adaptive(
            source_category=ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF,
            process_kind=ActivityMembershipSelectionProcessKind.FUSED_PROFILE_MOTIF,
            score_source_kind=ActivityMembershipScoreSourceKind.FUSED_PROFILE_MOTIF,
            selection_method=selection_method,
            score_source=score_source,
            threshold_top_k_policy=threshold_top_k_policy,
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=(
                selection_quantitative_matrix_fingerprint
            ),
            tested_quantitative_matrix_fingerprint=(
                tested_quantitative_matrix_fingerprint
            ),
            consumed_tested_matrix=consumed_tested_matrix,
            selected_kinase_universe=_string_tuple(
                selected_kinase_universe,
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
            selected_substrate_universe=_string_tuple(
                selected_substrate_universe,
                field_name="activity_membership_selection.selected_substrate_universe",
            ),
        )

    @classmethod
    def prediction_selected(
        cls,
        *,
        selection_method: str,
        score_source: str,
        score_source_kind: ActivityMembershipScoreSourceKind | str = (
            ActivityMembershipScoreSourceKind.PREDICTION_DERIVED
        ),
        threshold_top_k_policy: Mapping[str, object] | None = None,
        source_reference_fingerprints: tuple[TableFingerprint, ...] = (),
        selection_quantitative_matrix_fingerprint: TableFingerprint | None = None,
        tested_quantitative_matrix_fingerprint: TableFingerprint | None = None,
        consumed_tested_matrix: bool,
        data_adaptive_membership: bool,
        selected_kinase_universe: Iterable[object] = (),
        selected_substrate_universe: Iterable[object] = (),
    ) -> ActivityMembershipSelection:
        """Build current prediction-selected membership evidence."""

        return cls(
            source_category=ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED,
            selection_method=selection_method,
            selection_method_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
            score_source=score_source,
            membership_selection_schema_version=(
                ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
            ),
            threshold_top_k_policy=(
                {} if threshold_top_k_policy is None else threshold_top_k_policy
            ),
            selection_evidence=ActivityMembershipSelectionEvidence(
                selection_process_kind=(
                    ActivityMembershipSelectionProcessKind.PREDICTION_SELECTED
                ),
                selection_contract_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
                score_source_kind=_coerce_policy_enum(
                    ActivityMembershipScoreSourceKind,
                    score_source_kind,
                    field_name="activity_membership_selection.selection_evidence."
                    "score_source_kind",
                ),
                data_adaptive_membership=data_adaptive_membership,
                consumed_tested_matrix=consumed_tested_matrix,
                independence_evidence=None,
            ),
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=(
                selection_quantitative_matrix_fingerprint
            ),
            tested_quantitative_matrix_fingerprint=(
                tested_quantitative_matrix_fingerprint
            ),
            quantitative_dataset_fingerprint=None,
            consumed_tested_matrix=consumed_tested_matrix,
            selected_kinase_universe=_string_tuple(
                selected_kinase_universe,
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
            selected_substrate_universe=_string_tuple(
                selected_substrate_universe,
                field_name="activity_membership_selection.selected_substrate_universe",
            ),
        )

    @classmethod
    def _adaptive(
        cls,
        *,
        source_category: str,
        process_kind: ActivityMembershipSelectionProcessKind,
        score_source_kind: ActivityMembershipScoreSourceKind,
        selection_method: str,
        score_source: str,
        threshold_top_k_policy: Mapping[str, object] | None,
        source_reference_fingerprints: tuple[TableFingerprint, ...],
        selection_quantitative_matrix_fingerprint: TableFingerprint | None,
        tested_quantitative_matrix_fingerprint: TableFingerprint | None,
        consumed_tested_matrix: bool,
        selected_kinase_universe: Iterable[object],
        selected_substrate_universe: Iterable[object],
    ) -> ActivityMembershipSelection:
        return cls(
            source_category=source_category,
            selection_method=selection_method,
            selection_method_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
            score_source=score_source,
            membership_selection_schema_version=(
                ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION
            ),
            threshold_top_k_policy=(
                {} if threshold_top_k_policy is None else threshold_top_k_policy
            ),
            selection_evidence=ActivityMembershipSelectionEvidence(
                selection_process_kind=process_kind,
                selection_contract_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
                score_source_kind=_coerce_policy_enum(
                    ActivityMembershipScoreSourceKind,
                    score_source_kind,
                    field_name="activity_membership_selection.selection_evidence."
                    "score_source_kind",
                ),
                data_adaptive_membership=bool(consumed_tested_matrix),
                consumed_tested_matrix=consumed_tested_matrix,
                independence_evidence=None,
            ),
            source_reference_fingerprints=source_reference_fingerprints,
            selection_quantitative_matrix_fingerprint=(
                selection_quantitative_matrix_fingerprint
            ),
            tested_quantitative_matrix_fingerprint=(
                tested_quantitative_matrix_fingerprint
            ),
            quantitative_dataset_fingerprint=None,
            consumed_tested_matrix=consumed_tested_matrix,
            selected_kinase_universe=_string_tuple(
                selected_kinase_universe,
                field_name="activity_membership_selection.selected_kinase_universe",
            ),
            selected_substrate_universe=_string_tuple(
                selected_substrate_universe,
                field_name="activity_membership_selection.selected_substrate_universe",
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

    selection_evidence = _selection_evidence_model(selection.selection_evidence)
    process_kind = selection_evidence.selection_process_kind
    missing_evidence = _base_missing_evidence(selection)
    if process_kind == ActivityMembershipSelectionProcessKind.UNKNOWN:
        return _unavailable_decision(
            KSEA_MEMBERSHIP_MISSING_PROVENANCE_REASON,
            missing_evidence=("selection_evidence", *missing_evidence),
        )
    if process_kind == ActivityMembershipSelectionProcessKind.INCOMPLETE:
        return _unavailable_decision(
            KSEA_MEMBERSHIP_INCOMPLETE_INDEPENDENCE_EVIDENCE_REASON,
            missing_evidence=("selection_evidence", *missing_evidence),
        )
    if process_kind == ActivityMembershipSelectionProcessKind.PROFILE_DERIVED:
        return _unavailable_decision(
            (
                KSEA_MEMBERSHIP_CONSUMED_TESTED_MATRIX_REASON
                if selection.consumed_tested_matrix
                else KSEA_MEMBERSHIP_PROFILE_DERIVED_REASON
            ),
            missing_evidence=missing_evidence,
        )
    if process_kind == ActivityMembershipSelectionProcessKind.FUSED_PROFILE_MOTIF:
        return _unavailable_decision(
            (
                KSEA_MEMBERSHIP_CONSUMED_TESTED_MATRIX_REASON
                if selection.consumed_tested_matrix
                else KSEA_MEMBERSHIP_FUSED_PROFILE_MOTIF_REASON
            ),
            missing_evidence=missing_evidence,
        )
    if process_kind == ActivityMembershipSelectionProcessKind.PREDICTION_SELECTED:
        if (
            selection.consumed_tested_matrix
            or selection.selection_quantitative_matrix_fingerprint is not None
            or selection_evidence.data_adaptive_membership is True
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
    if process_kind == ActivityMembershipSelectionProcessKind.FIXED_EXTERNAL_REFERENCE:
        return _independent_reference_decision(
            selection,
            expected_policy=KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_FIXED_EXTERNAL_REFERENCE,
            extra_missing_evidence=missing_evidence,
        )
    if process_kind == ActivityMembershipSelectionProcessKind.SEQUENCE_ONLY_MOTIF:
        extra_missing = list(missing_evidence)
        if (
            selection_evidence.score_source_kind
            != ActivityMembershipScoreSourceKind.KINASE_LIBRARY_MOTIF
        ):
            extra_missing.append("selection_evidence.score_source_kind")
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


def _validate_membership_source_fact_coherence(
    selection: ActivityMembershipSelection,
) -> None:
    """Reject source-category states with contradictory membership facts."""

    evidence = _selection_evidence_model(selection.selection_evidence)
    process_kind = evidence.selection_process_kind
    if (
        evidence.selection_contract_version
        == ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION
    ):
        _validate_current_threshold_policy_has_no_decision_fields(selection)

    if process_kind == ActivityMembershipSelectionProcessKind.FIXED_EXTERNAL_REFERENCE:
        _validate_fixed_external_source_facts(
            selection,
        )
        return

    if process_kind == ActivityMembershipSelectionProcessKind.SEQUENCE_ONLY_MOTIF:
        _validate_sequence_only_source_facts(
            selection,
        )
        return

    if process_kind in {
        ActivityMembershipSelectionProcessKind.PROFILE_DERIVED,
        ActivityMembershipSelectionProcessKind.FUSED_PROFILE_MOTIF,
        ActivityMembershipSelectionProcessKind.PREDICTION_SELECTED,
    }:
        _validate_adaptive_source_facts(
            selection,
        )
        return

    if process_kind in {
        ActivityMembershipSelectionProcessKind.UNKNOWN,
        ActivityMembershipSelectionProcessKind.INCOMPLETE,
    }:
        _validate_unknown_source_facts(
            selection,
        )


def _validate_fixed_external_source_facts(
    selection: ActivityMembershipSelection,
) -> None:
    evidence = _selection_evidence_model(selection.selection_evidence)
    _require_current_selection_contract_version(
        evidence,
        source_category=selection.source_category,
    )
    if evidence.data_adaptive_membership is None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.data_adaptive_membership",
            "fixed_external_reference membership requires explicit non-adaptive evidence",
        )
    if evidence.data_adaptive_membership is True:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.data_adaptive_membership",
            "fixed_external_reference membership cannot be data-adaptive",
        )
    if selection.consumed_tested_matrix:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.consumed_tested_matrix",
            "fixed_external_reference membership cannot consume the tested matrix",
        )
    if selection.selection_quantitative_matrix_fingerprint is not None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_quantitative_matrix_fingerprint",
            (
                "fixed_external_reference membership cannot carry a selection "
                "quantitative-matrix fingerprint"
            ),
        )
    if (
        evidence.score_source_kind
        != ActivityMembershipScoreSourceKind.EXTERNAL_REFERENCE
    ):
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.score_source_kind",
            (
                "fixed_external_reference membership requires an external "
                "reference-derived score-source kind"
            ),
        )
    independence_evidence = evidence.independence_evidence
    if independence_evidence is None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.independence_evidence",
            "fixed_external_reference membership requires typed independence evidence",
        )
    _validate_independence_policy_for_category(
        source_category=selection.source_category,
        expected_policy=(
            ActivityMembershipIndependencePolicyKind.FIXED_EXTERNAL_REFERENCE
        ),
        independence_evidence=independence_evidence,
    )


def _validate_sequence_only_source_facts(
    selection: ActivityMembershipSelection,
) -> None:
    evidence = _selection_evidence_model(selection.selection_evidence)
    _require_current_selection_contract_version(
        evidence,
        source_category=selection.source_category,
    )
    if (
        evidence.score_source_kind
        != ActivityMembershipScoreSourceKind.KINASE_LIBRARY_MOTIF
    ):
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.score_source_kind",
            (
                "sequence_only_motif membership requires kinase-library motif "
                "score-source kind"
            ),
        )
    if evidence.data_adaptive_membership is None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.data_adaptive_membership",
            "sequence_only_motif membership requires explicit non-adaptive evidence",
        )
    if evidence.data_adaptive_membership is True:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.data_adaptive_membership",
            "sequence_only_motif membership cannot be data-adaptive",
        )
    if selection.consumed_tested_matrix:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.consumed_tested_matrix",
            "sequence_only_motif membership cannot consume the tested matrix",
        )
    if selection.selection_quantitative_matrix_fingerprint is not None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_quantitative_matrix_fingerprint",
            (
                "sequence_only_motif membership cannot carry a selection "
                "quantitative-matrix fingerprint"
            ),
        )
    independence_evidence = evidence.independence_evidence
    if independence_evidence is None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.independence_evidence",
            "sequence_only_motif membership requires typed independence evidence",
        )
    _validate_independence_policy_for_category(
        source_category=selection.source_category,
        expected_policy=ActivityMembershipIndependencePolicyKind.SEQUENCE_ONLY_MOTIF,
        independence_evidence=independence_evidence,
    )


def _validate_adaptive_source_facts(
    selection: ActivityMembershipSelection,
) -> None:
    evidence = _selection_evidence_model(selection.selection_evidence)
    _require_current_selection_contract_version(
        evidence,
        source_category=selection.source_category,
    )
    if evidence.data_adaptive_membership is None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.data_adaptive_membership",
            f"{selection.source_category} membership requires explicit adaptive-state evidence",
        )
    if evidence.independence_evidence is not None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.independence_evidence",
            (
                f"{selection.source_category} membership cannot carry "
                "fixed-external or sequence-only independence-policy evidence"
            ),
        )
    _validate_adaptive_score_source_kind(selection)


def _validate_unknown_source_facts(
    selection: ActivityMembershipSelection,
) -> None:
    evidence = _selection_evidence_model(selection.selection_evidence)
    if evidence.independence_evidence is not None:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.independence_evidence",
            (
                "unknown membership provenance cannot carry fixed-external or "
                "sequence-only independence-policy evidence"
            ),
        )


def _validate_independence_policy_for_category(
    *,
    source_category: str,
    expected_policy: ActivityMembershipIndependencePolicyKind,
    independence_evidence: ActivityMembershipIndependenceEvidence | None,
) -> None:
    if independence_evidence is None:
        return
    if independence_evidence.policy_kind != expected_policy:
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence."
            "independence_evidence.policy_kind",
            (
                f"{source_category} membership cannot carry independence-policy "
                f"kind {independence_evidence.policy_kind.value!r}; "
                f"expected {expected_policy.value!r}"
            ),
        )
    if (
        independence_evidence.policy_version
        != KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION
    ):
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence."
            "independence_evidence.policy_version",
            (
                f"{source_category} membership carries unsupported "
                "independence-policy version "
                f"{independence_evidence.policy_version!r}; "
                f"expected {KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION!r}"
            ),
        )


def _require_current_selection_contract_version(
    evidence: ActivityMembershipSelectionEvidence,
    *,
    source_category: str,
) -> None:
    if (
        evidence.selection_contract_version
        == ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION
    ):
        return
    _raise_membership_fact_contradiction(
        "activity_membership_selection.selection_evidence.selection_contract_version",
        (
            f"{source_category} membership carries unsupported selection-contract "
            f"version {evidence.selection_contract_version!r}; expected "
            f"{ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION!r}"
        ),
    )


def _validate_current_threshold_policy_has_no_decision_fields(
    selection: ActivityMembershipSelection,
) -> None:
    present = sorted(
        key
        for key in _KSEA_DECISION_POLICY_MAPPING_FIELDS
        if key in selection.threshold_top_k_policy
    )
    if not present:
        return
    joined = ", ".join(present)
    raise WorkflowBoundaryError(
        "activity_membership_selection.threshold_top_k_policy must not carry "
        "decision-bearing membership evidence in the current schema; move these "
        f"field(s) to selection_evidence: {joined}"
    )


def _validate_adaptive_score_source_kind(
    selection: ActivityMembershipSelection,
) -> None:
    evidence = _selection_evidence_model(selection.selection_evidence)
    if evidence.selection_process_kind == (
        ActivityMembershipSelectionProcessKind.PROFILE_DERIVED
    ):
        expected = ActivityMembershipScoreSourceKind.PROFILE_DERIVED
        if evidence.score_source_kind == expected:
            return
    elif evidence.selection_process_kind == (
        ActivityMembershipSelectionProcessKind.FUSED_PROFILE_MOTIF
    ):
        expected = ActivityMembershipScoreSourceKind.FUSED_PROFILE_MOTIF
        if evidence.score_source_kind == expected:
            return
    elif evidence.selection_process_kind == (
        ActivityMembershipSelectionProcessKind.PREDICTION_SELECTED
    ):
        allowed = {
            ActivityMembershipScoreSourceKind.PREDICTION_DERIVED,
            ActivityMembershipScoreSourceKind.PROFILE_DERIVED,
            ActivityMembershipScoreSourceKind.FUSED_PROFILE_MOTIF,
        }
        if evidence.score_source_kind in allowed:
            return
        allowed_text = ", ".join(sorted(item.value for item in allowed))
        _raise_membership_fact_contradiction(
            "activity_membership_selection.selection_evidence.score_source_kind",
            (
                "prediction_selected membership requires one of these "
                f"score-source kinds: {allowed_text}"
            ),
        )
        return
    else:
        return
    _raise_membership_fact_contradiction(
        "activity_membership_selection.selection_evidence.score_source_kind",
        (
            f"{selection.source_category} membership carries score-source kind "
            f"{evidence.score_source_kind.value!r}; expected {expected.value!r}"
        ),
    )


def _coerce_policy_enum(
    enum_type: type[_PolicyEnumT],
    value: object,
    *,
    field_name: str,
) -> _PolicyEnumT:
    try:
        return coerce_policy_enum(
            enum_type,
            value,
            field_name=field_name,
            error_type=WorkflowBoundaryError,
        )
    except TypeError as exc:
        raise WorkflowBoundaryError(f"{field_name} has unsupported enum type") from exc


def _coerce_selection_evidence(
    value: object,
    *,
    consumed_tested_matrix: bool,
) -> ActivityMembershipSelectionEvidence:
    if value is None:
        return ActivityMembershipSelectionEvidence(
            selection_process_kind=ActivityMembershipSelectionProcessKind.INCOMPLETE,
            selection_contract_version=ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION,
            score_source_kind=ActivityMembershipScoreSourceKind.UNKNOWN,
            data_adaptive_membership=None,
            consumed_tested_matrix=consumed_tested_matrix,
            independence_evidence=None,
        )
    if isinstance(value, ActivityMembershipSelectionEvidence):
        return value
    if isinstance(value, Mapping):
        return ActivityMembershipSelectionEvidence.from_payload(
            _require_mapping(
                value,
                field_name="activity_membership_selection.selection_evidence",
            )
        )
    raise WorkflowBoundaryError(
        "activity_membership_selection.selection_evidence must be "
        "ActivityMembershipSelectionEvidence, a mapping, or None"
    )


def _selection_evidence_model(
    value: ActivityMembershipSelectionEvidence | None,
) -> ActivityMembershipSelectionEvidence:
    if isinstance(value, ActivityMembershipSelectionEvidence):
        return value
    raise WorkflowBoundaryError(
        "activity_membership_selection.selection_evidence must be resolved before "
        "deriving KSEA membership eligibility"
    )


def _coerce_independence_evidence(
    value: object,
) -> ActivityMembershipIndependenceEvidence | None:
    if value is None:
        return None
    if isinstance(value, ActivityMembershipIndependenceEvidence):
        return value
    if isinstance(value, Mapping):
        return ActivityMembershipIndependenceEvidence.from_payload(
            _require_mapping(
                value,
                field_name="activity_membership_selection.selection_evidence."
                "independence_evidence",
            )
        )
    raise WorkflowBoundaryError(
        "activity_membership_selection.selection_evidence.independence_evidence "
        "must be ActivityMembershipIndependenceEvidence, a mapping, or None"
    )


def _source_category_from_selection_evidence(
    evidence: ActivityMembershipSelectionEvidence,
) -> str:
    return _PROCESS_KIND_SOURCE_CATEGORY[evidence.selection_process_kind]


def _require_supported_membership_payload_schema_version(value: object) -> str:
    schema_version = _require_non_empty_text(
        value,
        field_name=(
            "activity_membership_selection.membership_selection_schema_version"
        ),
    )
    if schema_version == ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION:
        return schema_version
    raise WorkflowBoundaryError(
        "activity_membership_selection.membership_selection_schema_version "
        "is unsupported; "
        f"expected {ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION!r}, "
        f"got {schema_version!r}"
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
    _compare_optional_decision_text(
        decision_payload,
        key="policy_version",
        expected=selection.inferential_decision.policy_version,
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
    evidence = _selection_evidence_model(selection.selection_evidence)
    if (
        evidence.selection_contract_version
        != ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION
    ):
        missing.append("selection_evidence.selection_contract_version")
    if evidence.data_adaptive_membership is not False:
        missing.append("selection_evidence.data_adaptive_membership")
    if evidence.consumed_tested_matrix:
        missing.append("selection_evidence.consumed_tested_matrix")
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
    independence_evidence = evidence.independence_evidence
    if independence_evidence is None:
        missing.append("selection_evidence.independence_evidence")
    else:
        if independence_evidence.policy_kind.value != expected_policy:
            missing.append("selection_evidence.independence_evidence.policy_kind")
        if (
            independence_evidence.policy_version
            != KSEA_MEMBERSHIP_INDEPENDENCE_POLICY_VERSION
        ):
            missing.append("selection_evidence.independence_evidence.policy_version")
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


def _raise_membership_fact_contradiction(
    field_name: str,
    detail: str,
) -> None:
    raise WorkflowBoundaryError(
        f"{field_name} contradicts activity_membership_selection.source_category: "
        f"{detail}"
    )


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
    "ACTIVITY_MEMBERSHIP_SELECTION_PAYLOAD_SCHEMA_VERSION",
    "ACTIVITY_MEMBERSHIP_SELECTION_POLICY_VERSION",
    "ACTIVITY_MEMBERSHIP_SOURCE_CATEGORIES",
    "ACTIVITY_MEMBERSHIP_SOURCE_FIXED_EXTERNAL_REFERENCE",
    "ACTIVITY_MEMBERSHIP_SOURCE_FUSED_PROFILE_MOTIF",
    "ACTIVITY_MEMBERSHIP_SOURCE_INCOMPLETE",
    "ACTIVITY_MEMBERSHIP_SOURCE_PREDICTION_SELECTED",
    "ACTIVITY_MEMBERSHIP_SOURCE_PROFILE_DERIVED",
    "ACTIVITY_MEMBERSHIP_SOURCE_SEQUENCE_ONLY_MOTIF",
    "ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN",
    "ActivityMembershipIndependenceEvidence",
    "ActivityMembershipIndependencePolicyKind",
    "ActivityMembershipSelection",
    "ActivityMembershipSelectionEvidence",
    "ActivityMembershipSelectionProcessKind",
    "ActivityMembershipScoreSourceKind",
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
