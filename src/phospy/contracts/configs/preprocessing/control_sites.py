"""Control-site constants for native SPS/RUV-style preprocessing correction."""

from __future__ import annotations

from typing import Literal

CONTROL_SITE_STATUS_CONTROL = "control"
CONTROL_SITE_STATUS_NON_CONTROL = "non_control"
CONTROL_SITE_STATUS_EXCLUDED = "excluded"
CONTROL_SITE_STATUS_UNKNOWN = "unknown"
CONTROL_SITE_STATUS_INVALID = "invalid"
ControlSiteStatusValue = Literal[
    "control",
    "non_control",
    "excluded",
    "unknown",
    "invalid",
]
CONTROL_SITE_STATUSES = frozenset(
    {
        CONTROL_SITE_STATUS_CONTROL,
        CONTROL_SITE_STATUS_NON_CONTROL,
        CONTROL_SITE_STATUS_EXCLUDED,
        CONTROL_SITE_STATUS_UNKNOWN,
        CONTROL_SITE_STATUS_INVALID,
    }
)

CONTROL_SITE_SOURCE_CALLER_SUPPLIED = "caller_supplied"
CONTROL_SITE_SELECTION_METHOD_CALLER_SUPPLIED = "caller_supplied"

__all__ = [
    "CONTROL_SITE_SELECTION_METHOD_CALLER_SUPPLIED",
    "CONTROL_SITE_SOURCE_CALLER_SUPPLIED",
    "CONTROL_SITE_STATUS_CONTROL",
    "CONTROL_SITE_STATUS_EXCLUDED",
    "CONTROL_SITE_STATUS_INVALID",
    "CONTROL_SITE_STATUS_NON_CONTROL",
    "CONTROL_SITE_STATUS_UNKNOWN",
    "CONTROL_SITE_STATUSES",
    "ControlSiteStatusValue",
]
