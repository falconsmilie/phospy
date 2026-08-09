"""Compatibility route for Kinase Library resource imports."""

from __future__ import annotations

from phospy.science.references.kinase_library_loading import (
    KinaseLibraryResourceLoader,
    ReferenceSourceTableReader,
    load_kinase_library_resource,
)
from phospy.science.references.kinase_library_models import (
    KinaseLibraryMatrix,
    KinaseLibraryPath,
    KinaseLibraryResidueClass,
    KinaseLibraryResource,
    KinaseLibraryResourceLoadRequest,
)
from phospy.science.references.validation.kinase_library import (
    KinaseLibraryResourceValidator,
)

__all__ = [
    "KinaseLibraryMatrix",
    "KinaseLibraryPath",
    "KinaseLibraryResidueClass",
    "KinaseLibraryResource",
    "KinaseLibraryResourceLoadRequest",
    "KinaseLibraryResourceLoader",
    "KinaseLibraryResourceValidator",
    "ReferenceSourceTableReader",
    "load_kinase_library_resource",
]
