from __future__ import annotations

from ..activities.analysis import KinaseActivityAnalyzer
from ..activities.results import KinaseActivityResult
from ..preprocessing.core import CoreProcessingResult
from ..validation.requests.pipeline import (
    PipelineInputs,
    validate_pipeline_runtime_compatibility,
)

__all__ = ["PipelineRuntimeService"]


class PipelineRuntimeService:
    """Execution service for the internal core+activity pipeline runtime path."""

    def run(
        self,
        *,
        request: PipelineInputs,
        kinase_activity_analyzer: KinaseActivityAnalyzer,
    ) -> tuple[CoreProcessingResult, KinaseActivityResult | None]:
        core = request.dataset.preprocessing.run(config=request.preprocessing_config)

        kinase_activity = None
        kinase_activity_request = validate_pipeline_runtime_compatibility(
            request=request,
            site_matrix=core.site_matrix.matrix,
        )
        if kinase_activity_request is not None:
            kinase_activity = kinase_activity_analyzer.run_validated(
                kinase_activity_request
            )
        return core, kinase_activity
