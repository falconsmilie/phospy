"""Internal validator for simple kinase workflow requests."""

from __future__ import annotations

from phospy.api.configs import (
    KinaseActivityConfig,
    KinasePredictionConfig,
    KinaseScoringConfig,
)
from phospy.api.requests import SimpleKinaseWorkflowRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.validation import WorkflowValidationError
from phospy.references.models import ReferenceBundle, ReferencePreset


class SimpleKinaseWorkflowValidator:
    """Validate `SimpleKinaseWorkflowRequest` before interpretation."""

    def run(self, request: SimpleKinaseWorkflowRequest) -> SimpleKinaseWorkflowRequest:
        if not isinstance(request, SimpleKinaseWorkflowRequest):
            raise WorkflowValidationError(
                "simple kinase workflow input must be a SimpleKinaseWorkflowRequest"
            )
        if not isinstance(request.dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "simple kinase workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        if not isinstance(request.references, (ReferencePreset, ReferenceBundle)):
            raise WorkflowValidationError(
                "simple kinase workflow request references must be ReferencePreset or ReferenceBundle"
            )
        if not isinstance(request.scoring_config, KinaseScoringConfig):
            raise WorkflowValidationError(
                "simple kinase workflow request scoring_config must be KinaseScoringConfig"
            )
        if not isinstance(request.prediction_config, KinasePredictionConfig):
            raise WorkflowValidationError(
                "simple kinase workflow request prediction_config must be KinasePredictionConfig"
            )
        if request.activity_config is not None and not isinstance(
            request.activity_config, KinaseActivityConfig
        ):
            raise WorkflowValidationError(
                "simple kinase workflow request activity_config must be KinaseActivityConfig or None"
            )
        return request
