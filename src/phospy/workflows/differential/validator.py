"""Internal validator for differential workflow requests."""

from __future__ import annotations

from typing import cast

from phospy.api.configs import DifferentialAnalysisConfig, MultipleTestingConfig
from phospy.api.requests import DifferentialAnalysisRequest
from phospy.errors.validation import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.differential.models import (
    ContrastMatrix,
    DesignMatrix,
    EmpiricalBayesConfig,
)
from phospy.science.differential.policy_models import TechnicalReplicatePolicy
from phospy.validation.common.dataframes import require_dataframe
from phospy.validation.workflows.differential import (
    DifferentialDatasetEligibilityValidator,
    ExperimentalDesignContractValidator,
)
from phospy.validation.workflows.identity import (
    DIFFERENTIAL_IDENTITY_CONTRACT,
    enforce_workflow_site_identity_contract,
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
        dataset_eligibility_validator: (
            DifferentialDatasetEligibilityValidator | None
        ) = None,
        technical_replicate_resolver: TechnicalReplicateResolver | None = None,
    ) -> None:
        self._design_validator = (
            design_validator or ExperimentalDesignContractValidator()
        )
        self._dataset_eligibility_validator = (
            dataset_eligibility_validator or DifferentialDatasetEligibilityValidator()
        )
        self._technical_replicate_resolver = (
            technical_replicate_resolver or TechnicalReplicateResolver()
        )

    def run(self, request: object) -> ValidatedDifferentialAnalysisRequest:
        if not isinstance(request, DifferentialAnalysisRequest):
            raise WorkflowValidationError(
                "differential workflow input must be a DifferentialAnalysisRequest"
            )
        if not isinstance(cast(object, request.dataset), AnalysisReadyPhosphoDataset):
            raise WorkflowValidationError(
                "differential workflow request dataset must be AnalysisReadyPhosphoDataset"
            )
        site_metadata = require_dataframe(
            request.dataset._borrow_site_metadata_frame(),  # pyright: ignore[reportPrivateUsage] - workflow boundary reads trusted internal dataset snapshots
            field_name="differential workflow request dataset.site_metadata",
            allow_empty=False,
            error_type=WorkflowValidationError,
        )
        enforce_workflow_site_identity_contract(
            site_metadata=site_metadata,
            field_name="differential workflow request dataset.site_metadata",
            contract=DIFFERENTIAL_IDENTITY_CONTRACT,
            error_type=WorkflowValidationError,
            allow_opaque_site_values=request.dataset.opaque_site_values_allowed,
        )
        self._dataset_eligibility_validator.run(dataset=request.dataset)
        config = request.config
        if not isinstance(cast(object, config), DifferentialAnalysisConfig):
            raise WorkflowValidationError(
                "differential workflow request config must be DifferentialAnalysisConfig"
            )
        technical_replicate_policy = config.technical_replicate_policy
        if not isinstance(
            cast(object, technical_replicate_policy), TechnicalReplicatePolicy
        ):
            raise WorkflowValidationError(
                "differential workflow request technical_replicate_policy must be "
                "TechnicalReplicatePolicy"
            )
        if not isinstance(cast(object, config.allow_design_subset), bool):
            raise WorkflowValidationError(
                "differential workflow request allow_design_subset must be a bool"
            )
        if not isinstance(
            cast(object, config.minimum_condition_replicates), int
        ) or isinstance(cast(object, config.minimum_condition_replicates), bool):
            raise WorkflowValidationError(
                "differential workflow request minimum_condition_replicates must be an int"
            )
        if config.minimum_condition_replicates < 1:
            raise WorkflowValidationError(
                "differential workflow request minimum_condition_replicates must be >= 1"
            )
        if not isinstance(cast(object, config.empirical_bayes), EmpiricalBayesConfig):
            raise WorkflowValidationError(
                "differential workflow request empirical_bayes must be EmpiricalBayesConfig"
            )
        if not isinstance(cast(object, config.multiple_testing), MultipleTestingConfig):
            raise WorkflowValidationError(
                "differential workflow request multiple_testing must be MultipleTestingConfig"
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
            allow_design_subset=config.allow_design_subset,
            minimum_condition_replicates=config.minimum_condition_replicates,
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
            config=config,
            workflow_provenance=technical_replicate_resolution.workflow_provenance,
        )


__all__ = ["DifferentialAnalysisValidator"]
