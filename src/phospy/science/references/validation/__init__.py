"""Compatibility facade for reference-manifest validation imports."""

from __future__ import annotations

from phospy.science.references.validation.bundle_semantics import (
    validate_bundled_reference_manifests as validate_bundled_reference_manifests,
)
from phospy.science.references.validation.bundle_semantics import (
    validate_reference_manifest as validate_reference_manifest,
)
from phospy.science.references.validation.manifest_schema import (
    _REQUIRED_FILE_FIELDS as _REQUIRED_FILE_FIELDS,  # noqa: F401
)
from phospy.science.references.validation.manifest_schema import (
    _REQUIRED_MANIFEST_FIELDS as _REQUIRED_MANIFEST_FIELDS,  # noqa: F401
)
from phospy.science.references.validation.manifest_schema import (
    load_reference_manifest as load_reference_manifest,
)
from phospy.science.references.validation.manifest_schema import (
    parse_reference_manifest_payload as parse_reference_manifest_payload,
)

__all__ = [
    "load_reference_manifest",
    "parse_reference_manifest_payload",
    "validate_bundled_reference_manifests",
    "validate_reference_manifest",
]
