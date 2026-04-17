"""Biological reference handling domain.

This package owns bundled reference assets, species and reference resolution,
substrate maps, motif resources, and site-sequence resources used across the
application.
"""

from .models import (
    ReferenceBundle,
    ReferenceBundleProvenance,
    ReferenceBundleSourceMetadata,
    ReferenceProvider,
)
from .resolution import (
    BundledReferenceProvider,
    normalize_bundled_reference,
    normalize_bundled_species,
)
from .resources import (
    build_reference_motif_sequences,
    bundled_reference_resource_path,
    load_bundled_site_sequences,
    load_bundled_substrate_map,
    load_grouped_mapping_file,
    load_string_mapping_file,
)

__all__ = [
    "BundledReferenceProvider",
    "ReferenceBundle",
    "ReferenceBundleProvenance",
    "ReferenceBundleSourceMetadata",
    "ReferenceProvider",
    "build_reference_motif_sequences",
    "bundled_reference_resource_path",
    "load_bundled_site_sequences",
    "load_bundled_substrate_map",
    "load_grouped_mapping_file",
    "load_string_mapping_file",
    "normalize_bundled_reference",
    "normalize_bundled_species",
]
