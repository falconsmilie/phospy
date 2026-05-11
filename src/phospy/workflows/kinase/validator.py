"""Internal validator for kinase workflow requests."""

from __future__ import annotations

from typing import cast

from phospy.api.requests import KinaseWorkflowRequest
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.errors.validation import WorkflowValidationError
from phospy.references.models import ReferenceBundle, ReferencePreset
from phospy.validation.common.dataframes import require_dataframe
from phospy.validation.datasets.site_metadata import enforce_localisation_requirement
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

    def run(self, request: object) -> KinaseWorkflowRequest:
        if not isinstance(request, KinaseWorkflowRequest):
            raise WorkflowValidationError(
                "kinase workflow input must be a KinaseWorkflowRequest"
            )
        dataset = cast(object, request.dataset)
        if not isinstance(dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "kinase workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        references = cast(object, request.references)
        if not isinstance(references, (ReferencePreset, ReferenceBundle)):
            raise WorkflowValidationError(
                "kinase workflow request references must be ReferencePreset or ReferenceBundle"
            )
        scoring_config, _, _ = self._config_validator.run(
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )
        reject_mixed_total_protein_quantitative_meaning(
            dataset=dataset,
            allow_mixed=scoring_config.allow_mixed_total_protein_quantitative_meaning,
            context="kinase workflow request dataset",
        )
        site_metadata = require_dataframe(
            dataset._borrow_site_metadata_frame(),
            field_name="kinase workflow request dataset.site_metadata",
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        enforce_localisation_requirement(
            site_metadata=site_metadata,
            field_name="kinase workflow request dataset.site_metadata",
            workflow_name="kinase workflow request",
            requirement=scoring_config.localisation_requirement,
            error_type=WorkflowValidationError,
        )
        return request
