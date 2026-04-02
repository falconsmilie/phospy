from __future__ import annotations

from .analysis import KinaseActivityRequest
from .pipeline import CorePipelineRequest
from .prediction import PredictionRequest
from .workflow import KinaseWorkflowRequest

"""Compatibility re-exports for legacy validation request imports.

New validation rules should live in the categorized boundary modules:
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
]
