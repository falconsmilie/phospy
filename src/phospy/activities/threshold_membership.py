"""Shared threshold-membership policy for activity substrate selection.

Policy ownership:
- This policy is defined by the activities domain.
- Activity methods consume this policy.
- Diagnostics report this policy and do not define threshold behavior.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.scoring.policy_models import ThresholdMode

THRESHOLD_MEMBERSHIP_MODE = ThresholdMode.GREATER_THAN_OR_EQUAL


@dataclass(frozen=True, slots=True)
class ActivityThresholdMembershipPolicy:
    """Standard activity threshold-membership policy metadata."""

    mode: ThresholdMode
    operator: str
    description: str

    @property
    def rule(self) -> str:
        return self.mode.value

    def to_payload(self) -> dict[str, str]:
        return {
            "operator": self.operator,
            "rule": self.rule,
            "description": self.description,
        }


@dataclass(frozen=True, slots=True)
class ActivityThresholdMembershipDiagnostics:
    """Method-level threshold membership diagnostics metadata."""

    threshold_parameter: str
    threshold_value: float
    operator: str
    rule: str
    description: str

    def to_payload(self) -> dict[str, object]:
        return {
            "threshold_parameter": self.threshold_parameter,
            "threshold_value": float(self.threshold_value),
            "operator": self.operator,
            "rule": self.rule,
            "description": self.description,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityThresholdMembershipDiagnostics:
        threshold_parameter = str(payload.get("threshold_parameter", "")).strip()
        operator = str(payload.get("operator", "")).strip()
        rule = str(payload.get("rule", "")).strip()
        description = str(payload.get("description", "")).strip()
        threshold_value = _resolve_threshold_value(payload)
        if not threshold_parameter:
            raise ValueError("threshold_parameter must be a non-empty string")
        if not operator:
            raise ValueError("operator must be a non-empty string")
        if not rule:
            raise ValueError("rule must be a non-empty string")
        if not description:
            raise ValueError("description must be a non-empty string")
        return cls(
            threshold_parameter=threshold_parameter,
            threshold_value=threshold_value,
            operator=operator,
            rule=rule,
            description=description,
        )


def resolve_activity_threshold_membership_policy(
    threshold_mode: ThresholdMode | str = THRESHOLD_MEMBERSHIP_MODE,
) -> ActivityThresholdMembershipPolicy:
    mode = ThresholdMode.parse(
        threshold_mode,
        field_name="activity threshold membership mode",
    )
    if mode is ThresholdMode.GREATER_THAN_OR_EQUAL:
        return ActivityThresholdMembershipPolicy(
            mode=mode,
            operator=">=",
            description=("scores greater than or equal to the threshold are included"),
        )
    if mode is ThresholdMode.GREATER_THAN:
        return ActivityThresholdMembershipPolicy(
            mode=mode,
            operator=">",
            description="scores greater than the threshold are included",
        )
    return ActivityThresholdMembershipPolicy(
        mode=mode,
        operator=">=",
        description="scores greater than or equal to the threshold are included",
    )


def build_activity_threshold_membership_diagnostics(
    *,
    threshold_parameter: str,
    threshold_value: float,
    threshold_mode: ThresholdMode | str = THRESHOLD_MEMBERSHIP_MODE,
) -> ActivityThresholdMembershipDiagnostics:
    policy = resolve_activity_threshold_membership_policy(threshold_mode)
    return ActivityThresholdMembershipDiagnostics(
        threshold_parameter=str(threshold_parameter),
        threshold_value=float(threshold_value),
        operator=policy.operator,
        rule=policy.rule,
        description=policy.description,
    )


_THRESHOLD_MEMBERSHIP_POLICY = resolve_activity_threshold_membership_policy()
THRESHOLD_MEMBERSHIP_RULE = _THRESHOLD_MEMBERSHIP_POLICY.rule
THRESHOLD_MEMBERSHIP_OPERATOR = _THRESHOLD_MEMBERSHIP_POLICY.operator
THRESHOLD_MEMBERSHIP_DESCRIPTION = _THRESHOLD_MEMBERSHIP_POLICY.description


def threshold_membership_mask_array(
    scores: np.ndarray,
    *,
    threshold: float,
    threshold_mode: ThresholdMode | str = THRESHOLD_MEMBERSHIP_MODE,
) -> np.ndarray:
    """Return membership mask using the standard activity-threshold rule."""

    threshold_value = float(threshold)
    return _threshold_comparison_mask(
        scores=scores,
        threshold_value=threshold_value,
        threshold_mode=threshold_mode,
    )


def threshold_membership_mask_frame(
    scores: pd.DataFrame,
    *,
    threshold: float,
    threshold_mode: ThresholdMode | str = THRESHOLD_MEMBERSHIP_MODE,
) -> pd.DataFrame:
    """Return a DataFrame mask using the standard activity-threshold rule."""

    mask = threshold_membership_mask_array(
        scores.to_numpy(dtype=float, copy=False),
        threshold=threshold,
        threshold_mode=threshold_mode,
    )
    return pd.DataFrame(
        mask,
        index=pd.Index(scores.index),
        columns=pd.Index(scores.columns),
    )


def threshold_membership_filtered_frame(
    scores: pd.DataFrame,
    *,
    threshold: float,
    threshold_mode: ThresholdMode | str = THRESHOLD_MEMBERSHIP_MODE,
) -> pd.DataFrame:
    """Return score values for members and NaN for non-members."""

    score_values = scores.to_numpy(dtype=float, copy=False)
    mask_values = threshold_membership_mask_array(
        score_values,
        threshold=threshold,
        threshold_mode=threshold_mode,
    )
    filtered_values = np.where(mask_values, score_values, np.nan)
    return pd.DataFrame(
        filtered_values,
        index=pd.Index(scores.index),
        columns=pd.Index(scores.columns),
    )


def _resolve_threshold_value(payload: Mapping[str, object]) -> float:
    raw_value = payload.get("threshold_value")
    if raw_value is None:
        return float("nan")
    if isinstance(raw_value, (int, float)):
        return float(raw_value)
    if isinstance(raw_value, str):
        try:
            return float(raw_value)
        except ValueError as exc:
            raise ValueError("threshold_value must be numeric") from exc
    raise ValueError("threshold_value must be numeric")


def _threshold_comparison_mask(
    *,
    scores: np.ndarray,
    threshold_value: float,
    threshold_mode: ThresholdMode | str,
) -> np.ndarray:
    mode = ThresholdMode.parse(
        threshold_mode,
        field_name="activity threshold membership mode",
    )
    finite = np.isfinite(scores)
    if mode is ThresholdMode.GREATER_THAN_OR_EQUAL:
        return finite & (scores >= threshold_value)
    if mode is ThresholdMode.GREATER_THAN:
        return finite & (scores > threshold_value)
    return finite & (scores >= threshold_value)


__all__ = [
    "ActivityThresholdMembershipDiagnostics",
    "ActivityThresholdMembershipPolicy",
    "THRESHOLD_MEMBERSHIP_DESCRIPTION",
    "THRESHOLD_MEMBERSHIP_MODE",
    "THRESHOLD_MEMBERSHIP_OPERATOR",
    "THRESHOLD_MEMBERSHIP_RULE",
    "build_activity_threshold_membership_diagnostics",
    "resolve_activity_threshold_membership_policy",
    "threshold_membership_filtered_frame",
    "threshold_membership_mask_array",
    "threshold_membership_mask_frame",
]
