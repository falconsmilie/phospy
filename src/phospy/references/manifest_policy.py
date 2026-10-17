"""Reference manifest policy constants and enums."""

from __future__ import annotations

from enum import Enum

REFERENCE_MANIFEST_SCHEMA_VERSION = "1.1"


class RedistributionStatus(str, Enum):
    """Machine-readable redistribution status for reference manifests."""

    APPROVED = "approved"
    EXTERNAL_ONLY = "external_only"
    UNRESOLVED = "unresolved"


class RedistributionEvidenceType(str, Enum):
    """Machine-readable evidence category for exact-file redistribution approval."""

    UPSTREAM_PACKAGE_LICENSE = "upstream_package_license"


__all__ = [
    "REFERENCE_MANIFEST_SCHEMA_VERSION",
    "RedistributionEvidenceType",
    "RedistributionStatus",
]
