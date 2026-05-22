"""Internal executor for differential workflow requests."""

from __future__ import annotations

from typing import cast

from phospy.contracts.configs import MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.differential.executor import (
    DifferentialAnalysisExecutor as DifferentialComputationExecutor,
)
from phospy.science.differential.models import DifferentialAnalysisResult
from phospy.workflows.differential.models import (
    InterpretedDifferentialAnalysisRequest,
)


class DifferentialAnalysisExecutor:
    """Run differential computation for interpreted execution inputs."""

    def __init__(
        self,
        *,
        computation_executor: DifferentialComputationExecutor | None = None,
    ) -> None:
        self._computation_executor = (
            computation_executor or DifferentialComputationExecutor()
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
        if (
            request.config.multiple_testing.method
            != MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG
        ):
            raise WorkflowBoundaryError(
                seam="differential.executor.multiple_testing_method",
                next_action=(
                    "use multiple_testing.method='benjamini_hochberg' for this release"
                ),
                details={"method": request.config.multiple_testing.method},
                message_prefix="differential workflow boundary validation failed",
            )
        result = self._computation_executor.run(request.computation_request)
        if (
            request.workflow_provenance is None
            and request.policy_provenance is None
            and request.dataset_preprocessing_report is None
        ):
            return result
        return DifferentialAnalysisResult._from_owned(  # pyright: ignore[reportPrivateUsage] - trusted internal ownership-preserving constructor
            residual_variance=result.residual_variance,
            posterior_residual_variance=result.posterior_residual_variance,
            prior_residual_variance=result.prior_residual_variance,
            prior_degrees_of_freedom_series_value=(
                result.prior_degrees_of_freedom_series_value
            ),
            prior_variance=result.prior_variance,
            prior_degrees_of_freedom=result.prior_degrees_of_freedom,
            residual_degrees_of_freedom=result.residual_degrees_of_freedom,
            empirical_bayes_method=result.empirical_bayes_method,
            empirical_bayes_robust=result.empirical_bayes_robust,
            empirical_bayes_trend=result.empirical_bayes_trend,
            prior_diagnostics=result.prior_diagnostics,
            mean_variance_trend_diagnostics=result.mean_variance_trend_diagnostics,
            policy_provenance=request.policy_provenance,
            contrast_tables=result._contrast_tables,  # pyright: ignore[reportPrivateUsage] - owned contrast tables are forwarded without copying
            workflow_provenance=request.workflow_provenance,
            input_dataset_preprocessing_report=request.dataset_preprocessing_report,
        )


__all__ = ["DifferentialAnalysisExecutor"]
