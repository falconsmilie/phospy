"""Reference-manifest specific exceptions."""

from __future__ import annotations

from phospy.core.errors.references import ReferenceResolutionError


class ReferenceManifestError(ReferenceResolutionError):
    """A reference manifest is missing, malformed, or unverifiable."""


__all__ = ["ReferenceManifestError"]
