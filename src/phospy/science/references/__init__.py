"""Reference domain package."""

from __future__ import annotations

from typing import TYPE_CHECKING

__all__ = [
    "BundledReferenceLane",
    "Organism",
    "ReferenceBundle",
    "ReferenceBundleBuildRequest",
    "ReferenceBundleBuilder",
    "ReferenceManifest",
    "ReferencePreset",
    "SequenceWindowDefinition",
]

if TYPE_CHECKING:
    from phospy.science.references.builder import ReferenceBundleBuilder
    from phospy.science.references.models import (
        BundledReferenceLane,
        Organism,
        ReferenceBundle,
        ReferenceBundleBuildRequest,
        ReferenceManifest,
        ReferencePreset,
        SequenceWindowDefinition,
    )


def __getattr__(name: str) -> object:
    if name == "ReferenceBundleBuilder":
        from phospy.science.references.builder import ReferenceBundleBuilder

        return ReferenceBundleBuilder
    if name in __all__:
        from phospy.science.references import models as _models

        return getattr(_models, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
