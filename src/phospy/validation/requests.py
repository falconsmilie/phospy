from __future__ import annotations

from .analysis import KinaseActivityRequest, ValidatedAnalysisRequest
from .pipeline import CorePipelineRequest, ValidatedPipelineRequest
from .prediction import PredictionRequest
from .workflow import KinaseWorkflowRequest, ValidatedWorkflowRequest

"""Compatibility re-exports for legacy validation request imports.

New validation rules should live in the categorized boundary modules:
- phospy.validation.dataset
- phospy.validation.pipeline
- phospy.validation.workflow
- phospy.validation.analysis
- phospy.validation.prediction
"""

__all__ = [
    "CorePipelineRequest",
    "KinaseActivityRequest",
    "KinaseWorkflowRequest",
    "PredictionRequest",
    "ValidatedAnalysisRequest",
    "ValidatedPipelineRequest",
    "ValidatedWorkflowRequest",
]
