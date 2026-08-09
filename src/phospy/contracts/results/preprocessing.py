"""Public preprocessing result contract aliases."""

from __future__ import annotations

from phospy.science.datasets.preprocessing.batch_correction_models import (
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.protein_aware_models import (
    ProteinAwareMappingDiagnostics,
    ProteinAwarePreparationReport,
    ProteinAwarePreparationResult,
    ProteinAwareSiteEligibility,
)

__all__ = [
    "BatchCorrectionDiagnostics",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "ProteinAwareMappingDiagnostics",
    "ProteinAwarePreparationReport",
    "ProteinAwarePreparationResult",
    "ProteinAwareSiteEligibility",
]
