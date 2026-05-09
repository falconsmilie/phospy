"""Differential-analysis domain exports."""

from phospy.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    DifferentialAnalysisRequest,
    DifferentialAnalysisResult,
    EmpiricalBayesConfig,
)
from phospy.differential.public import DifferentialAnalysis

__all__ = [
    "ContrastMatrix",
    "DesignMatrix",
    "DifferentialAnalysis",
    "DifferentialAnalysisRequest",
    "DifferentialAnalysisResult",
    "EmpiricalBayesConfig",
]
