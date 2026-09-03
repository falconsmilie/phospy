"""Internal executor coordinator for differential workflow requests."""

from __future__ import annotations

from collections.abc import Mapping
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
from phospy.science.differential.protein_covariate_adjusted import (
    ProteinCovariateAdjustedDifferentialKernel,
)
from phospy.workflows.differential.eligibility import (
    DifferentialComputationEligibilityResolver,
    DifferentialPostFitEligibilityResolver,
)
from phospy.workflows.differential.fitting import DifferentialModelFitter
from phospy.workflows.differential.models import InterpretedDifferentialAnalysisRequest
from phospy.workflows.differential.protein_aware_eligibility import (
    protein_aware_feature_eligibility_after_computation,
)
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
        protein_aware_kernel: ProteinCovariateAdjustedDifferentialKernel | None = None,
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
        self._protein_aware_kernel = (
            protein_aware_kernel or ProteinCovariateAdjustedDifferentialKernel()
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

        if request.execution_config.protein_aware_method is not None:
            return self._run_protein_aware(request)

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

    def _run_protein_aware(
        self,
        request: InterpretedDifferentialAnalysisRequest,
    ) -> DifferentialAnalysisResult:
        resolved_inputs = request.protein_aware_inputs
        if resolved_inputs is None:
            raise WorkflowBoundaryError(
                seam="differential.executor.protein_aware_inputs",
                next_action=(
                    "pass interpreter output with resolved protein-aware inputs into "
                    "DifferentialAnalysisExecutor.run"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        if (
            request.execution_config.paired_design_policy
            == PAIRED_DESIGN_POLICY_DUPLICATE_CORRELATION
        ):
            raise WorkflowBoundaryError(
                seam="differential.executor.protein_aware_duplicate_correlation",
                next_action=(
                    "use fixed_block for paired protein-aware differential analysis; "
                    "duplicate-correlation modelling is not supported for "
                    "protein-aware adjustment"
                ),
                message_prefix="differential workflow boundary validation failed",
            )
        try:
            computation_result = self._protein_aware_kernel.run(
                resolved_inputs.computation_request
            )
        except PhosPyInputError as error:
            raise WorkflowBoundaryError(
                seam="differential.executor.protein_aware_fit",
                next_action=(
                    "provide at least one protein-aware eligible phosphosite with "
                    "usable matched total-protein covariates and an admissible "
                    "augmented design"
                ),
                details=_protein_aware_fit_error_details(error),
                message_prefix="differential protein-aware fitting failed",
            ) from error

        failed_site_ids = tuple(
            str(value)
            for value in computation_result.site_failure_diagnostics.index.tolist()
        )
        feature_eligibility_inputs = (
            protein_aware_feature_eligibility_after_computation(
                resolved_inputs=resolved_inputs,
                computation_result=computation_result,
            )
        )
        workflow_provenance = self._provenance_assembler.run(
            workflow_provenance=request.workflow_provenance,
            input_feature_ids=resolved_inputs.full_site_ids,
            model_fit_feature_ids=resolved_inputs.tested_site_ids,
            failed_model_fit_feature_ids=failed_site_ids,
            multiple_testing_feature_ids=tuple(
                str(value) for value in computation_result.tested_site_ids
            ),
            imputation_policy_inputs=request.imputation_policy_inputs,
            feature_eligibility_inputs=feature_eligibility_inputs,
        )
        return self._result_assembler.run_protein_aware(
            request=request,
            resolved_inputs=resolved_inputs,
            computation_result=computation_result,
            workflow_provenance=workflow_provenance,
            feature_eligibility_inputs=feature_eligibility_inputs,
        )


def _protein_aware_fit_error_details(error: PhosPyInputError) -> dict[str, object]:
    details: dict[str, object] = {"error": str(error)}
    diagnostics = error.diagnostics
    if diagnostics is None:
        return details
    details["diagnostics_type"] = type(diagnostics).__name__
    if not isinstance(diagnostics, Mapping):
        return details

    diagnostic_mapping = cast(Mapping[object, object], diagnostics)
    details["diagnostic_keys"] = tuple(str(key) for key in diagnostic_mapping)
    row_counts: dict[str, int] = {}
    for key, value in diagnostic_mapping.items():
        row_count = _diagnostic_row_count(value)
        if row_count is not None:
            row_counts[str(key)] = row_count
        if str(key) in {"status_counts", "reason_counts"} and isinstance(
            value,
            Mapping,
        ):
            details[str(key)] = {
                str(name): int(count)
                for name, count in cast(Mapping[object, int], value).items()
            }
    if row_counts:
        details["diagnostic_row_counts"] = row_counts
    return details


def _diagnostic_row_count(value: object) -> int | None:
    shape = getattr(value, "shape", None)
    if not isinstance(shape, tuple) or not shape:
        return None
    shape_tuple = cast(tuple[object, ...], shape)
    row_count = shape_tuple[0]
    if not isinstance(row_count, int):
        return None
    return int(row_count)


__all__ = ["DifferentialAnalysisExecutor"]
