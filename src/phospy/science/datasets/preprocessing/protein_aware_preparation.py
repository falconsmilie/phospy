"""Compatibility route for protein-aware preparation imports."""

from __future__ import annotations

from phospy.science.datasets.preprocessing.protein_aware_models import (
    PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS,
    PROTEIN_AWARE_MATCHED_PAIR_COLUMNS,
    PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS,
    PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION,
    PROTEIN_AWARE_SITE_ELIGIBILITY_COLUMNS,
    ProteinAwareMappingDiagnostics,
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)
from phospy.science.datasets.preprocessing.stages.protein_aware_preparation import (
    ProteinAwarePreparationStage,
)

__all__ = [
    "PROTEIN_AWARE_AMBIGUOUS_MAPPING_DIAGNOSTIC_COLUMNS",
    "PROTEIN_AWARE_MATCHED_PAIR_COLUMNS",
    "PROTEIN_AWARE_MISSING_PROTEIN_ABUNDANCE_DIAGNOSTIC_COLUMNS",
    "PROTEIN_AWARE_PREPARATION_SCHEMA_VERSION",
    "PROTEIN_AWARE_SITE_ELIGIBILITY_COLUMNS",
    "ProteinAwareMappingDiagnostics",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
    "ProteinAwarePreparationStage",
    "ProteinAwareSiteEligibility",
]
