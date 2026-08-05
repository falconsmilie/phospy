"""Compatibility import route for reference manifest models.

Owned implementations are split by responsibility under this package. This
module preserves the historical ``phospy.science.references.manifest`` import
path by re-exporting the same class, enum, and constant objects.
"""

from phospy.science.references.manifest_files import (
    ReferenceFileManifest,
    SequenceWindowDefinition,
)
from phospy.science.references.manifest_model import ReferenceManifest
from phospy.science.references.manifest_policy import (
    REFERENCE_MANIFEST_SCHEMA_VERSION,
    RedistributionEvidenceType,
    RedistributionStatus,
)
from phospy.science.references.redistribution import (
    RedistributionAttribution,
    RedistributionEvidence,
    RedistributionScope,
    UpstreamPackageLicenseEvidence,
)

__all__ = [
    "REFERENCE_MANIFEST_SCHEMA_VERSION",
    "ReferenceFileManifest",
    "ReferenceManifest",
    "RedistributionAttribution",
    "RedistributionEvidence",
    "RedistributionEvidenceType",
    "RedistributionScope",
    "RedistributionStatus",
    "SequenceWindowDefinition",
    "UpstreamPackageLicenseEvidence",
]
