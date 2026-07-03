"""Structured kinase workflow site attrition metrics."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

import pandas as pd

from phospy.contracts.configs import KinaseAttritionPolicy
from phospy.errors.workflows import PhosPyWorkflowError

KINASE_ATTRITION_POLICY_CAVEAT_CODE = "kinase_attrition_policy_violation"
KINASE_ATTRITION_POLICY_OUTCOME_PASSED = "passed"
KINASE_ATTRITION_POLICY_OUTCOME_WARNED = "warned"
KINASE_ATTRITION_POLICY_OUTCOME_FAILED = "failed"
KinaseAttritionPolicyOutcome = Literal["passed", "warned", "failed"]


@dataclass(frozen=True, slots=True)
class KinaseAttritionMetrics:
    """Pre-scoring kinase workflow site attrition counts and fractions."""

    total_dataset_sites: int
    reference_overlap_sites: int
    sequence_supported_sites: int
    scored_sites: int
    reference_overlap_fraction: float
    sequence_supported_fraction: float
    scored_fraction: float

    @classmethod
    def from_counts(
        cls,
        *,
        total_dataset_sites: int,
        reference_overlap_sites: int,
        sequence_supported_sites: int,
        scored_sites: int,
    ) -> KinaseAttritionMetrics:
        total = int(total_dataset_sites)
        if total <= 0:
            raise PhosPyWorkflowError(
                "kinase workflow internal attrition validation failed at "
                "seam=kinase.attrition_metrics.total_dataset_sites; "
                "total_dataset_sites must be greater than zero before "
                "fractional attrition metrics can be computed"
            )
        reference_overlap = _validate_count(
            value=reference_overlap_sites,
            field_name="reference_overlap_sites",
            total=total,
        )
        sequence_supported = _validate_count(
            value=sequence_supported_sites,
            field_name="sequence_supported_sites",
            total=total,
        )
        scored = _validate_count(
            value=scored_sites,
            field_name="scored_sites",
            total=total,
        )
        return cls(
            total_dataset_sites=total,
            reference_overlap_sites=reference_overlap,
            sequence_supported_sites=sequence_supported,
            scored_sites=scored,
            reference_overlap_fraction=reference_overlap / total,
            sequence_supported_fraction=sequence_supported / total,
            scored_fraction=scored / total,
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "total_dataset_sites": int(self.total_dataset_sites),
            "reference_overlap_sites": int(self.reference_overlap_sites),
            "sequence_supported_sites": int(self.sequence_supported_sites),
            "scored_sites": int(self.scored_sites),
            "reference_overlap_fraction": float(self.reference_overlap_fraction),
            "sequence_supported_fraction": float(self.sequence_supported_fraction),
            "scored_fraction": float(self.scored_fraction),
        }


@dataclass(frozen=True, slots=True)
class KinaseAttritionPolicyViolation:
    """A configured attrition-policy threshold violation."""

    threshold_name: str
    configured_threshold: float
    observed_value: float
    total_dataset_sites: int
    reference_overlap_sites: int
    sequence_supported_sites: int
    scored_sites: int
    message: str

    def to_payload(self) -> dict[str, object]:
        return {
            "code": KINASE_ATTRITION_POLICY_CAVEAT_CODE,
            "threshold_name": self.threshold_name,
            "configured_threshold": float(self.configured_threshold),
            "observed_value": float(self.observed_value),
            "total_dataset_sites": int(self.total_dataset_sites),
            "reference_overlap_sites": int(self.reference_overlap_sites),
            "sequence_supported_sites": int(self.sequence_supported_sites),
            "scored_sites": int(self.scored_sites),
            "message": self.message,
        }


def evaluate_kinase_attrition_policy(
    *,
    metrics: KinaseAttritionMetrics,
    policy: KinaseAttritionPolicy,
) -> tuple[KinaseAttritionPolicyViolation, ...]:
    """Evaluate policy thresholds against one shared metrics object."""

    checks = (
        (
            "minimum_reference_overlap_fraction",
            float(policy.minimum_reference_overlap_fraction),
            float(metrics.reference_overlap_fraction),
            "Kinase reference overlap",
            "reference_overlap_sites",
            int(metrics.reference_overlap_sites),
        ),
        (
            "minimum_sequence_supported_fraction",
            float(policy.minimum_sequence_supported_fraction),
            float(metrics.sequence_supported_fraction),
            "Kinase site-sequence support",
            "sequence_supported_sites",
            int(metrics.sequence_supported_sites),
        ),
        (
            "minimum_scored_fraction",
            float(policy.minimum_scored_fraction),
            float(metrics.scored_fraction),
            "Kinase scoring",
            "scored_sites",
            int(metrics.scored_sites),
        ),
    )
    violations: list[KinaseAttritionPolicyViolation] = []
    for (
        threshold_name,
        configured_threshold,
        observed_value,
        stage_label,
        retained_count_name,
        retained_count,
    ) in checks:
        if observed_value >= configured_threshold:
            continue
        message = (
            f"{stage_label} retained {_format_percent(observed_value)} of "
            "dataset sites, below the configured "
            f"{threshold_name}={_format_percent(configured_threshold)} "
            f"({retained_count_name}={retained_count}, "
            f"total_dataset_sites={int(metrics.total_dataset_sites)})."
        )
        violations.append(
            KinaseAttritionPolicyViolation(
                threshold_name=threshold_name,
                configured_threshold=configured_threshold,
                observed_value=observed_value,
                total_dataset_sites=int(metrics.total_dataset_sites),
                reference_overlap_sites=int(metrics.reference_overlap_sites),
                sequence_supported_sites=int(metrics.sequence_supported_sites),
                scored_sites=int(metrics.scored_sites),
                message=message,
            )
        )
    return tuple(violations)


def kinase_attrition_policy_to_payload(
    policy: KinaseAttritionPolicy,
) -> dict[str, object]:
    """Return the configured kinase attrition policy as JSON-safe provenance."""

    return {
        "minimum_reference_overlap_fraction": float(
            policy.minimum_reference_overlap_fraction
        ),
        "minimum_sequence_supported_fraction": float(
            policy.minimum_sequence_supported_fraction
        ),
        "minimum_scored_fraction": float(policy.minimum_scored_fraction),
        "on_violation": policy.on_violation,
    }


def resolve_kinase_attrition_policy_outcome(
    *,
    policy: KinaseAttritionPolicy,
    violations: tuple[KinaseAttritionPolicyViolation, ...],
) -> KinaseAttritionPolicyOutcome:
    """Resolve pass/warn/fail status from the configured policy and violations."""

    if not violations:
        return KINASE_ATTRITION_POLICY_OUTCOME_PASSED
    if policy.on_violation == "error":
        return KINASE_ATTRITION_POLICY_OUTCOME_FAILED
    return KINASE_ATTRITION_POLICY_OUTCOME_WARNED


def build_kinase_attrition_provenance_payload(
    *,
    metrics: KinaseAttritionMetrics,
    policy: KinaseAttritionPolicy,
    violations: tuple[KinaseAttritionPolicyViolation, ...],
) -> dict[str, object]:
    """Build the structured attrition payload shared by results and provenance."""

    warning_messages = tuple(violation.message for violation in violations)
    return {
        "metrics": metrics.to_payload(),
        "policy": kinase_attrition_policy_to_payload(policy),
        "policy_outcome": resolve_kinase_attrition_policy_outcome(
            policy=policy,
            violations=violations,
        ),
        "policy_violations": [violation.to_payload() for violation in violations],
        "warning_messages": list(warning_messages),
    }


def build_kinase_attrition_metrics(
    *,
    dataset_site_index: pd.Index,
    reference_site_index: Iterable[object],
    sequence_supported_site_index: pd.Index,
) -> KinaseAttritionMetrics:
    """Build attrition metrics from resolved site identity sets."""

    dataset_sites = set(_string_values(dataset_site_index))
    reference_sites = set(_string_values(reference_site_index))
    sequence_supported_sites = set(_string_values(sequence_supported_site_index))
    reference_overlap_sites = dataset_sites.intersection(reference_sites)
    scored_sites = reference_overlap_sites.intersection(sequence_supported_sites)
    return KinaseAttritionMetrics.from_counts(
        total_dataset_sites=len(dataset_sites),
        reference_overlap_sites=len(reference_overlap_sites),
        sequence_supported_sites=len(
            dataset_sites.intersection(sequence_supported_sites)
        ),
        scored_sites=len(scored_sites),
    )


def build_kinase_attrition_metrics_from_overlap(
    *,
    total_dataset_sites: int,
    reference_overlap_site_ids: Iterable[object],
    sequence_supported_site_index: pd.Index,
) -> KinaseAttritionMetrics:
    """Build attrition metrics from pre-resolved reference-overlap identities."""

    reference_overlap_sites = set(_string_values(reference_overlap_site_ids))
    sequence_supported_sites = set(_string_values(sequence_supported_site_index))
    scored_sites = reference_overlap_sites.intersection(sequence_supported_sites)
    return KinaseAttritionMetrics.from_counts(
        total_dataset_sites=int(total_dataset_sites),
        reference_overlap_sites=len(reference_overlap_sites),
        sequence_supported_sites=len(sequence_supported_sites),
        scored_sites=len(scored_sites),
    )


def _validate_count(*, value: int, field_name: str, total: int) -> int:
    count = int(value)
    if count < 0:
        raise PhosPyWorkflowError(
            "kinase workflow internal attrition validation failed at "
            f"seam=kinase.attrition_metrics.{field_name}; "
            f"{field_name} must be greater than or equal to zero"
        )
    if count > total:
        raise PhosPyWorkflowError(
            "kinase workflow internal attrition validation failed at "
            f"seam=kinase.attrition_metrics.{field_name}; "
            f"{field_name} must not exceed total_dataset_sites; "
            f"{field_name}={count}; total_dataset_sites={total}"
        )
    return count


def _string_values(values: Iterable[object]) -> list[str]:
    return [str(value) for value in values]


def _format_percent(value: float) -> str:
    return f"{float(value) * 100:.1f}%"


__all__ = [
    "KINASE_ATTRITION_POLICY_CAVEAT_CODE",
    "KINASE_ATTRITION_POLICY_OUTCOME_FAILED",
    "KINASE_ATTRITION_POLICY_OUTCOME_PASSED",
    "KINASE_ATTRITION_POLICY_OUTCOME_WARNED",
    "KinaseAttritionMetrics",
    "KinaseAttritionPolicyOutcome",
    "KinaseAttritionPolicyViolation",
    "build_kinase_attrition_provenance_payload",
    "build_kinase_attrition_metrics",
    "build_kinase_attrition_metrics_from_overlap",
    "evaluate_kinase_attrition_policy",
    "kinase_attrition_policy_to_payload",
    "resolve_kinase_attrition_policy_outcome",
]
