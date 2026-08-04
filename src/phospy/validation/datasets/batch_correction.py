"""Compatibility facade for private batch-correction dataset validation.

Focused implementation lives in sibling modules for design metadata,
applied provenance, and control-site provenance checks. This module keeps the
legacy private import path available for internal callers.
"""

from __future__ import annotations

from phospy.validation.datasets.batch_correction_controls import (
    normalize_applied_selected_site_key_rows,
)
from phospy.validation.datasets.batch_correction_design import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
    BatchStructureValidator,
    ConditionStructureValidator,
    DesignRankValidator,
    ReplicateStructureDiagnosticHelper,
    ReplicateStructureDiagnostics,
    ReplicateStructureValidator,
    ResolvedBatchDesignMetadata,
    SampleMetadataAlignmentValidator,
)
from phospy.validation.datasets.batch_correction_provenance import (
    validate_applied_native_sps_ruv_correction_provenance,
)

__all__ = [
    "BatchCorrectionAdequacyValidator",
    "BatchDesignMetadataValidator",
    "BatchStructureValidator",
    "ConditionStructureValidator",
    "DesignRankValidator",
    "ReplicateStructureDiagnosticHelper",
    "ReplicateStructureDiagnostics",
    "ReplicateStructureValidator",
    "ResolvedBatchDesignMetadata",
    "SampleMetadataAlignmentValidator",
    "normalize_applied_selected_site_key_rows",
    "validate_applied_native_sps_ruv_correction_provenance",
]
