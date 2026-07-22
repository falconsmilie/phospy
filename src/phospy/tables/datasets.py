"""Compatibility re-exports for dataset table schemas."""

from phospy.science.tables.datasets import (
    PhosphoIntensityMatrix,
    SampleMetadataTable,
    SiteMetadataTable,
    TotalProteinMatrix,
)
from phospy.science.tables.datasets import (
    _build_identity_coherence_frame as _build_identity_coherence_frame,
)
from phospy.science.tables.datasets import (
    _display_ids_are_canonical as _display_ids_are_canonical,
)
from phospy.science.tables.datasets import (
    _drop_signalome_grouping_metadata as _drop_signalome_grouping_metadata,
)

__all__ = [
    "PhosphoIntensityMatrix",
    "SampleMetadataTable",
    "SiteMetadataTable",
    "TotalProteinMatrix",
]
