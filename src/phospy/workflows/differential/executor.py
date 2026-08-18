"""Internal executor coordinator for differential workflow requests."""

from __future__ import annotations

from typing import cast

from phospy.contracts.configs.differential import (
    PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION,
)
from phospy.errors.input import PhosPyInputError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.executor import (
    DuplicateCorrelationDifferentialAnalysisExecutor,
)
from phospy.science.differential.models import DifferentialAnalysisResult
from phospy.workflows.differential.eligibility import (
    DifferentialComputationEligibilityResolver,
    DifferentialPostFitEligibilityResolver,
)
from phospy.workflows.differential.fitting import DifferentialModelFitter
from phospy.workflows.differential.models import InterpretedDifferentialAnalysisRequest
from phospy.workflows.differential.provenance import (
    DifferentialWorkflowProvenanceAssembler,
    build_duplicate_correlation_workflow_provenance,
)
from phospy.workflows.differential.result_assembly import DifferentialResultAssembler


class DifferentialAnalysisExecutor:
    """Coordinate eligible differential fitting and public result assembly."""

    def __init__(
        self,
        *,
        computation_executor: DifferentialComputationExecutor | None = None,
        post_fit_eligibility_resolver: DifferentialPostFitEligibilityResolver
        | None = None,
        eligibility_resolver: DifferentialComputationEligibilityResolver | None = None,
        model_fitter: DifferentialModelFitter | None = None,
        duplicate_correlation_executor: (
            DuplicateCorrelationDifferentialAnalysisExecutor | None
        ) = None,
        result_assembler: DifferentialResultAssembler | None = None,
        provenance_assembler: DifferentialWorkflowProvenanceAssembler | None = None,
    ) -> None:
        self._eligibility_resolver = (
            eligibility_resolver
            or DifferentialComputationEligibilityResolver(
                post_fit_eligibility_resolver=post_fit_eligibility_resolver,
            )
        )
        self._model_fitter = model_fitter or DifferentialModelFitter(
            computation_executor=computation_executor,
        )
        self._duplicate_correlation_executor = (
            duplicate_correlation_executor
            or DuplicateCorrelationDifferentialAnalysisExecutor()
        )
        self._result_assembler = result_assembler or DifferentialResultAssembler()
        self._provenance_assembler = (
            provenance_assembler or DifferentialWorkflowProvenanceAssembler()
        )

    def run(
        self, request: InterpretedDifferentialAnalysisRequest
    ) -> DifferentialAnalysisResult:
        if not isinstance(
            cast(object, request), InterpretedDifferentialAnalysisRequest
        ):
            raise WorkflowBoundaryError(
                seam="differential.executor.interpreted_request_type",
                next_action=(
                    "pass interpreter output into DifferentialAnalysisExecutor.run"
                ),
                message_prefix="differential workflow boundary validation failed",
            )

        eligibility = self._eligibility_resolver.run(request)
        duplicate_correlation_provenance = None
        if (
            request.execution_config.paired_design_policy
            == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
        ):
            execution_design = request.execution_design
            if execution_design is None or execution_design.block_ids is None:
                raise WorkflowBoundaryError(
                    seam="differential.executor.duplicate_correlation_blocks",
                    next_action=(
                        "interpret duplicate_correlation with a validated block_id "
                        "vector before execution"
                    ),
                    message_prefix="differential workflow boundary validation failed",
                )
            try:
                duplicate_fit = self._duplicate_correlation_executor.run(
                    eligibility.computation_request,
                    block_ids=execution_design.block_ids,
                )
            except PhosPyInputError as error:
                raise WorkflowBoundaryError(
                    seam="differential.executor.duplicate_correlation_fit",
                    next_action=(
                        "provide a full-rank non-block fixed-effects design, "
                        "explicit repeated block IDs, and enough eligible features "
                        "for REML duplicate-correlation estimation"
                    ),
                    details={"error": str(error)},
                    message_prefix=(
                        "differential duplicate-correlation fitting failed"
                    ),
                ) from error
            computation_result = duplicate_fit.computation_result
            duplicate_correlation_provenance = (
                build_duplicate_correlation_workflow_provenance(
                    request=request,
                    computation_request=eligibility.computation_request,
                    consensus_result=duplicate_fit.consensus,
                    gls_fit=duplicate_fit.gls_fit,
                    imputation_policy_inputs=request.imputation_policy_inputs,
                    feature_eligibility_inputs=eligibility.feature_eligibility_inputs,
                )
            )
        else:
            computation_result = self._model_fitter.run(eligibility.computation_request)
        workflow_provenance = self._provenance_assembler.run(
            workflow_provenance=request.workflow_provenance,
            input_feature_ids=eligibility.input_feature_ids,
            model_fit_feature_ids=eligibility.model_fit_feature_ids,
            failed_model_fit_feature_ids=eligibility.failed_model_fit_feature_ids,
            multiple_testing_feature_ids=eligibility.multiple_testing_feature_ids,
            imputation_policy_inputs=request.imputation_policy_inputs,
            feature_eligibility_inputs=eligibility.feature_eligibility_inputs,
            duplicate_correlation=duplicate_correlation_provenance,
        )
        return self._result_assembler.run(
            request=request,
            computation_result=computation_result,
            eligibility=eligibility,
            workflow_provenance=workflow_provenance,
            duplicate_correlation=duplicate_correlation_provenance,
        )


__all__ = ["DifferentialAnalysisExecutor"]
