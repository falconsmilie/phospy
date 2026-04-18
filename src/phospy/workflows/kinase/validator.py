"""Internal validator for kinase workflow requests."""

from __future__ import annotations

from phospy.api.requests import KinaseWorkflowRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.validation import WorkflowValidationError
from phospy.references.models import ReferenceBundle, ReferencePreset
from phospy.validation.references.compatibility import ReferenceCompatibilityValidator
from phospy.validation.workflows.configs import WorkflowConfigValidator


class KinaseWorkflowValidator:
    """Validate `KinaseWorkflowRequest` before interpretation."""

    def __init__(
        self,
        *,
        config_validator: WorkflowConfigValidator | None = None,
        reference_compatibility: ReferenceCompatibilityValidator | None = None,
    ) -> None:
        self._config_validator = config_validator or WorkflowConfigValidator()
        self._reference_compatibility = (
            reference_compatibility or ReferenceCompatibilityValidator()
        )

    def run(self, request: KinaseWorkflowRequest) -> KinaseWorkflowRequest:
        if not isinstance(request, KinaseWorkflowRequest):
            raise WorkflowValidationError(
                "kinase workflow input must be a KinaseWorkflowRequest"
            )
        if not isinstance(request.dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "kinase workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        if not isinstance(request.references, (ReferencePreset, ReferenceBundle)):
            raise WorkflowValidationError(
                "kinase workflow request references must be ReferencePreset or ReferenceBundle"
            )
        self._reference_compatibility.run(
            request.references,
            dataset_organism=request.dataset.organism,
        )
        self._config_validator.run_kinase_scoring(request.scoring_config)
        self._config_validator.run_kinase_prediction(request.prediction_config)
        self._config_validator.run_kinase_activity(request.activity_config)
        return request
