"""Reference domain package."""

from __future__ import annotations

from typing import TYPE_CHECKING

__all__ = [
    "BundledReferenceLane",
    "KinaseLibraryMatrix",
    "KinaseLibraryResidueClass",
    "KinaseLibraryResource",
    "KinaseLibraryResourceLoadRequest",
    "KinaseLibraryResourceLoader",
    "Organism",
    "ReferenceBundle",
    "ReferenceBundleBuildRequest",
    "ReferenceBundleBuilder",
    "ReferenceBundleMissingValueCount",
    "ReferenceBundleSourceFileValidationReport",
    "ReferenceBundleTableValidationReport",
    "ReferenceBundleValidationReport",
    "ReferenceFileManifest",
    "ReferenceManifest",
    "ReferencePreset",
    "SequenceWindowDefinition",
    "load_kinase_library_resource",
]

if TYPE_CHECKING:
    from phospy.science.references.builder import ReferenceBundleBuilder
    from phospy.science.references.kinase_library import (
        KinaseLibraryMatrix,
        KinaseLibraryResidueClass,
        KinaseLibraryResource,
        KinaseLibraryResourceLoader,
        KinaseLibraryResourceLoadRequest,
        load_kinase_library_resource,
    )
    from phospy.science.references.models import (
        BundledReferenceLane,
        Organism,
        ReferenceBundle,
        ReferenceBundleBuildRequest,
        ReferenceBundleMissingValueCount,
        ReferenceBundleSourceFileValidationReport,
        ReferenceBundleTableValidationReport,
        ReferenceBundleValidationReport,
        ReferenceFileManifest,
        ReferenceManifest,
        ReferencePreset,
        SequenceWindowDefinition,
    )


def __getattr__(name: str) -> object:
    if name == "ReferenceBundleBuilder":
        from phospy.science.references.builder import ReferenceBundleBuilder

        return ReferenceBundleBuilder
    if name in {
        "KinaseLibraryMatrix",
        "KinaseLibraryResidueClass",
        "KinaseLibraryResource",
        "KinaseLibraryResourceLoadRequest",
        "KinaseLibraryResourceLoader",
        "load_kinase_library_resource",
    }:
        from phospy.science.references import kinase_library as _kinase_library

        return getattr(_kinase_library, name)
    if name in __all__:
        from phospy.science.references import models as _models

        return getattr(_models, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
