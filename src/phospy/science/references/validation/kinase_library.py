"""Kinase Library resource validation service."""

from __future__ import annotations

from phospy.science.references.kinase_library_models import (
    KinaseLibraryResource,
    _validate_kinase_library_resource_contract,
)


class KinaseLibraryResourceValidator:
    """Validate the stable Kinase Library-style resource contract."""

    def run(self, resource: KinaseLibraryResource) -> None:
        _validate_kinase_library_resource_contract(resource)


__all__ = ["KinaseLibraryResourceValidator"]
