"""IO adapters for Kinase Library-style local reference resources."""

from __future__ import annotations

from phospy.io.bundles.reference_sources import ReferenceSourceTableReader
from phospy.science.references.kinase_library import (
    KinaseLibraryPath,
    KinaseLibraryResource,
    KinaseLibraryResourceLoadRequest,
)
from phospy.science.references.kinase_library import (
    KinaseLibraryResourceLoader as _KinaseLibraryResourceLoader,
)


class KinaseLibraryResourceLoader(_KinaseLibraryResourceLoader):
    """Public local-file loader wired with the concrete source table reader."""

    def __init__(self) -> None:
        self._init_components(source_reader=None)

    @classmethod
    def _with_components(
        cls,
        *,
        source_reader: ReferenceSourceTableReader | None = None,
    ) -> KinaseLibraryResourceLoader:
        loader = cls.__new__(cls)
        loader._init_components(source_reader=source_reader)
        return loader

    def _init_components(
        self,
        *,
        source_reader: ReferenceSourceTableReader | None = None,
    ) -> None:
        super().__init__(source_reader=source_reader or ReferenceSourceTableReader())


def load_kinase_library_resource(
    request: KinaseLibraryResourceLoadRequest | KinaseLibraryPath,
) -> KinaseLibraryResource:
    """Load a local Kinase Library-style resource with default IO wiring."""

    return KinaseLibraryResourceLoader().run(request)


__all__ = [
    "KinaseLibraryResourceLoader",
    "load_kinase_library_resource",
]
