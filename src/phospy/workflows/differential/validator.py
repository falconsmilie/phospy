"""Internal validator for differential workflow requests."""

from __future__ import annotations

from phospy.api.requests import (
    DifferentialAnalysisRequest,
    TechnicalReplicatePolicy,
)
from phospy.datasets.models import AnalysisReadyPhosphoDataset
from phospy.differential.models import ContrastMatrix, DesignMatrix
from phospy.errors.validation import WorkflowValidationError
from phospy.validation.workflows.differential import (
    ExperimentalDesignContractValidator,
)
from phospy.workflows.differential.models import (
    ValidatedDifferentialAnalysisRequest,
)
from phospy.workflows.differential.replicates import TechnicalReplicateResolver


class DifferentialAnalysisValidator:
    """Validate `DifferentialAnalysisRequest` before interpretation."""

    def __init__(
        self,
        *,
        design_validator: ExperimentalDesignContractValidator | None = None,
        technical_replicate_resolver: TechnicalReplicateResolver | None = None,
    ) -> None:
        self._design_validator = (
            design_validator or ExperimentalDesignContractValidator()
        )
        self._technical_replicate_resolver = (
            technical_replicate_resolver or TechnicalReplicateResolver()
        )

    def run(self, request: object) -> ValidatedDifferentialAnalysisRequest:
        if not isinstance(request, DifferentialAnalysisRequest):
            raise WorkflowValidationError(
                "differential workflow input must be a DifferentialAnalysisRequest"
            )
        if not isinstance(request.dataset, AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "differential workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        technical_replicate_policy = request.technical_replicate_policy
        if not isinstance(technical_replicate_policy, TechnicalReplicatePolicy):
            raise WorkflowValidationError(
                "differential workflow request technical_replicate_policy must be "
                "TechnicalReplicatePolicy"
            )
        technical_replicate_resolution = self._technical_replicate_resolver.run(
            dataset=request.dataset,
            design=request.design,
            technical_replicate_policy=technical_replicate_policy,
        )
        validated_design_contract = self._design_validator.run(
            dataset=technical_replicate_resolution.dataset,
            design=technical_replicate_resolution.design,
            contrasts=request.contrasts,
            allow_design_subset=request.allow_design_subset,
            minimum_condition_replicates=request.minimum_condition_replicates,
        )

        design_matrix = DesignMatrix(validated_design_contract.design_frame)
        contrast_matrix = ContrastMatrix(validated_design_contract.contrast_frame)
        return ValidatedDifferentialAnalysisRequest(
            dataset=technical_replicate_resolution.dataset,
            design=validated_design_contract.design,
            contrasts=validated_design_contract.contrasts,
            analysis_sample_ids=validated_design_contract.analysis_sample_ids,
            design_matrix=design_matrix,
            contrast_matrix=contrast_matrix,
            empirical_bayes=request.empirical_bayes,
            multiple_testing=request.multiple_testing,
            workflow_provenance=technical_replicate_resolution.workflow_provenance,
        )


__all__ = ["DifferentialAnalysisValidator"]
