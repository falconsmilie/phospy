"""Internal validator for kinase workflow requests."""

from __future__ import annotations

from phospy.api.requests import KinaseWorkflowRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.validation import WorkflowValidationError
from phospy.references.models import ReferenceBundle, ReferencePreset
from phospy.validation.workflows.configs import (
    KinaseWorkflowConfigValidator,
    reject_mixed_total_protein_quantitative_meaning,
)


class KinaseWorkflowValidator:
    """Validate `KinaseWorkflowRequest` before interpretation."""

    def __init__(
        self,
        *,
        config_validator: KinaseWorkflowConfigValidator | None = None,
    ) -> None:
        self._config_validator = config_validator or KinaseWorkflowConfigValidator()

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
        scoring_config, _, _ = self._config_validator.run(
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )
        reject_mixed_total_protein_quantitative_meaning(
            dataset=request.dataset,
            allow_mixed=scoring_config.allow_mixed_total_protein_quantitative_meaning,
            context="kinase workflow request dataset",
        )
        return request
