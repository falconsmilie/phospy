"""Internal validator for kinase workflow requests."""

from __future__ import annotations

from typing import cast

from phospy.contracts.configs import (
    KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES,
    KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY,
)
from phospy.contracts.requests import KinaseWorkflowRequest
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.kinase_library import KinaseLibraryResource
from phospy.science.references.models import ReferenceBundle, ReferencePreset
from phospy.validation.common.dataframes import require_dataframe
from phospy.validation.datasets.site_metadata import (
    enforce_localisation_requirement,
)
from phospy.validation.workflows.configs import (
    KinaseWorkflowConfigValidator,
    reject_mixed_total_protein_quantitative_meaning,
)
from phospy.validation.workflows.identity import (
    KINASE_IDENTITY_CONTRACT,
    enforce_workflow_site_identity_contract,
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
        if (
            request.reference_display_ambiguity_policy
            not in KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES
        ):
            supported = ", ".join(sorted(KINASE_REFERENCE_DISPLAY_AMBIGUITY_POLICIES))
            raise WorkflowValidationError(
                "kinase workflow request reference_display_ambiguity_policy "
                f"must be one of: {supported}"
            )
        scoring_config, _, _ = self._config_validator.run(
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )
        if scoring_config.scoring_mode in KINASE_SCORING_MODES_REQUIRING_KINASE_LIBRARY:
            if request.kinase_library_resource is None:
                raise WorkflowValidationError(
                    "kinase workflow request kinase_library_resource is required "
                    f"when scoring_config.scoring_mode={scoring_config.scoring_mode!r}"
                )
            if not isinstance(request.kinase_library_resource, KinaseLibraryResource):
                raise WorkflowValidationError(
                    "kinase workflow request kinase_library_resource must be "
                    "KinaseLibraryResource when Kinase Library scoring is selected"
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
        enforce_workflow_site_identity_contract(
            site_metadata=site_metadata,
            expected_index=dataset._borrow_phospho_frame().index,
            expected_index_field_name="kinase workflow request dataset.phospho.index",
            field_name="kinase workflow request dataset.site_metadata",
            contract=KINASE_IDENTITY_CONTRACT,
            error_type=WorkflowValidationError,
            allow_opaque_site_values=dataset.opaque_site_values_allowed,
        )
        enforce_localisation_requirement(
            site_metadata=site_metadata,
            field_name="kinase workflow request dataset.site_metadata",
            workflow_name="kinase workflow request",
            requirement=scoring_config.localisation_requirement,
            error_type=WorkflowValidationError,
        )
        return request
