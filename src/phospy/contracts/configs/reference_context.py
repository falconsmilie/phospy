"""Reference-context compatibility policy configuration."""

from __future__ import annotations

from phospy.policies import PolicyEnum


class ReferenceContextCompatibilityPolicy(PolicyEnum):
    """Workflow policy for unknown biological reference-context compatibility."""

    REQUIRE_KNOWN_MATCH = "require_known_match"
    ALLOW_UNKNOWN_WITH_CAVEAT = "allow_unknown_with_caveat"


REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH = (
    ReferenceContextCompatibilityPolicy.REQUIRE_KNOWN_MATCH
)
REFERENCE_CONTEXT_COMPATIBILITY_POLICY_ALLOW_UNKNOWN_WITH_CAVEAT = (
    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
)
REFERENCE_CONTEXT_COMPATIBILITY_POLICIES = frozenset(
    {
        REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH,
        REFERENCE_CONTEXT_COMPATIBILITY_POLICY_ALLOW_UNKNOWN_WITH_CAVEAT,
    }
)


__all__ = [
    "REFERENCE_CONTEXT_COMPATIBILITY_POLICIES",
    "REFERENCE_CONTEXT_COMPATIBILITY_POLICY_ALLOW_UNKNOWN_WITH_CAVEAT",
    "REFERENCE_CONTEXT_COMPATIBILITY_POLICY_REQUIRE_KNOWN_MATCH",
    "ReferenceContextCompatibilityPolicy",
]
