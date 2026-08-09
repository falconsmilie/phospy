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
    "KinaseLibraryResourceValidator",
    "Organism",
    "ReferenceBundle",
    "ReferenceBundleBuildRequest",
    "ReferenceBundleBuilder",
    "ReferenceBundleMissingValueCount",
    "ReferenceBundleSourceFileValidationReport",
    "ReferenceBundleTableValidationReport",
    "ReferenceBundleValidationReport",
    "ReferenceContext",
    "ReferenceFileManifest",
    "ReferenceManifest",
    "ReferencePreset",
    "RedistributionStatus",
    "SequenceWindowDefinition",
    "load_kinase_library_resource",
    "reference_context_from_provenance",
]

if TYPE_CHECKING:
    from phospy.science.references.builder import ReferenceBundleBuilder
    from phospy.science.references.kinase_library_loading import (
        KinaseLibraryResourceLoader,
        load_kinase_library_resource,
    )
    from phospy.science.references.kinase_library_models import (
        KinaseLibraryMatrix,
        KinaseLibraryResidueClass,
        KinaseLibraryResource,
        KinaseLibraryResourceLoadRequest,
    )
    from phospy.science.references.models import (
        BundledReferenceLane,
        Organism,
        RedistributionStatus,
        ReferenceBundle,
        ReferenceBundleBuildRequest,
        ReferenceBundleMissingValueCount,
        ReferenceBundleSourceFileValidationReport,
        ReferenceBundleTableValidationReport,
        ReferenceBundleValidationReport,
        ReferenceContext,
        ReferenceFileManifest,
        ReferenceManifest,
        ReferencePreset,
        SequenceWindowDefinition,
        reference_context_from_provenance,
    )
    from phospy.science.references.validation.kinase_library import (
        KinaseLibraryResourceValidator,
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
    }:
        from phospy.science.references import kinase_library_models as _kinase_library

        return getattr(_kinase_library, name)
    if name in {
        "KinaseLibraryResourceLoader",
        "load_kinase_library_resource",
    }:
        from phospy.science.references import kinase_library_loading as _kinase_library

        return getattr(_kinase_library, name)
    if name == "KinaseLibraryResourceValidator":
        from phospy.science.references.validation.kinase_library import (
            KinaseLibraryResourceValidator,
        )

        return KinaseLibraryResourceValidator
    if name in __all__:
        from phospy.science.references import models as _models

        return getattr(_models, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
