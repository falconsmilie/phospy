"""Public localisation-confidence policy models for workflow validation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.contracts.configs.common import _require_real_between
from phospy.errors.validation import WorkflowValidationError

LOCALISATION_POLICY_ALLOW_UNKNOWN = "allow_unknown"
LOCALISATION_POLICY_REQUIRE_PRESENT = "require_present"
LOCALISATION_POLICY_REQUIRE_THRESHOLD = "require_threshold"
LocalisationPolicy = Literal[
    "allow_unknown",
    "require_present",
    "require_threshold",
]
LOCALISATION_POLICIES = frozenset(
    {
        LOCALISATION_POLICY_ALLOW_UNKNOWN,
        LOCALISATION_POLICY_REQUIRE_PRESENT,
        LOCALISATION_POLICY_REQUIRE_THRESHOLD,
    }
)


@dataclass(frozen=True, slots=True)
class LocalisationRequirement:
    """Workflow-facing localisation-confidence requirement."""

    require_present: bool = False
    minimum_probability: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.require_present, bool):
            raise WorkflowValidationError(
                "localisation_requirement.require_present must be a bool"
            )
        if self.minimum_probability is not None:
            _require_real_between(
                self.minimum_probability,
                field_name="localisation_requirement.minimum_probability",
                minimum=0.0,
                maximum=1.0,
                error_type=WorkflowValidationError,
            )

    @property
    def policy(self) -> LocalisationPolicy:
        if self.minimum_probability is not None:
            return LOCALISATION_POLICY_REQUIRE_THRESHOLD
        if self.require_present:
            return LOCALISATION_POLICY_REQUIRE_PRESENT
        return LOCALISATION_POLICY_ALLOW_UNKNOWN

    @property
    def requires_probability_column(self) -> bool:
        return self.require_present or self.minimum_probability is not None


__all__ = [
    "LOCALISATION_POLICIES",
    "LOCALISATION_POLICY_ALLOW_UNKNOWN",
    "LOCALISATION_POLICY_REQUIRE_PRESENT",
    "LOCALISATION_POLICY_REQUIRE_THRESHOLD",
    "LocalisationPolicy",
    "LocalisationRequirement",
]
