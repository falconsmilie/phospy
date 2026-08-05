"""Compatibility import route for batch-correction preprocessing models."""

from phospy.science.datasets.preprocessing.batch_correction_engine import (
    BatchCorrectionEngine,
    LinearResidualizeBatchCorrectionEngine,
)
from phospy.science.datasets.preprocessing.batch_correction_models import (
    BATCH_CORRECTION_CONFOUNDING_CONFOUNDED,
    BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
    BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED,
    BATCH_CORRECTION_CONFOUNDING_PASSED,
    BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS,
    BATCH_CORRECTION_STATUS_APPLIED,
    BATCH_CORRECTION_STATUS_DISABLED,
    BATCH_CORRECTION_STATUS_REJECTED,
    BatchCorrectionConfoundingCheckStatus,
    BatchCorrectionDiagnostics,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
    BatchCorrectionResult,
    BatchCorrectionStatus,
    MatrixShape,
)

__all__ = [
    "BATCH_CORRECTION_CONFOUNDING_CONFOUNDED",
    "BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE",
    "BATCH_CORRECTION_CONFOUNDING_NOT_CHECKED",
    "BATCH_CORRECTION_CONFOUNDING_PASSED",
    "BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS",
    "BATCH_CORRECTION_STATUS_APPLIED",
    "BATCH_CORRECTION_STATUS_DISABLED",
    "BATCH_CORRECTION_STATUS_REJECTED",
    "BatchCorrectionConfoundingCheckStatus",
    "BatchCorrectionDiagnostics",
    "BatchCorrectionEngine",
    "BatchCorrectionPolicy",
    "BatchCorrectionReport",
    "BatchCorrectionResult",
    "BatchCorrectionStatus",
    "LinearResidualizeBatchCorrectionEngine",
    "MatrixShape",
]
