"""Batch-correction workflow orchestrator."""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

from phospy.provenance.models import JsonValue
from phospy.science.batch_correction import SpsRuvStyleExecutor
from phospy.validation.workflows.batch_correction.control_site_workflow import (
    BatchCorrectionWorkflowControlSiteValidator,
)
from phospy.validation.workflows.batch_correction.design import (
    BatchCorrectionWorkflowDesignValidator,
)
from phospy.validation.workflows.batch_correction.missingness import (
    BatchCorrectionWorkflowMissingnessValidator,
)
from phospy.validation.workflows.batch_correction.request import (
    BatchCorrectionWorkflowRequestValidator,
)
from phospy.validation.workflows.batch_correction.stage_order import (
    BatchCorrectionWorkflowStageOrderValidator,
)
from phospy.workflows.batch_correction.contracts import (
    BatchCorrectionControlSiteValidatorContract,
    BatchCorrectionDesignValidatorContract,
    BatchCorrectionExecutorContract,
    BatchCorrectionInterpreterContract,
    BatchCorrectionMissingnessValidatorContract,
    BatchCorrectionProvenanceRecorderContract,
    BatchCorrectionRequestValidatorContract,
    BatchCorrectionStageOrderValidatorContract,
    BatchCorrectionWorkflowResult,
)
from phospy.workflows.batch_correction.interpreter import BatchCorrectionPlanInterpreter
from phospy.workflows.batch_correction.provenance import (
    BatchCorrectionProvenanceRecorder,
)


class BatchCorrectionWorkflow:
    """Coordinate validation, interpretation, execution, and provenance."""

    def __init__(
        self,
        *,
        request_validator: BatchCorrectionRequestValidatorContract | None = None,
        design_validator: BatchCorrectionDesignValidatorContract | None = None,
        control_site_validator: (
            BatchCorrectionControlSiteValidatorContract | None
        ) = None,
        stage_order_validator: (
            BatchCorrectionStageOrderValidatorContract | None
        ) = None,
        missingness_validator: (
            BatchCorrectionMissingnessValidatorContract | None
        ) = None,
        interpreter: BatchCorrectionInterpreterContract | None = None,
        executor: BatchCorrectionExecutorContract | None = None,
        provenance_recorder: BatchCorrectionProvenanceRecorderContract | None = None,
    ) -> None:
        self._request_validator = (
            request_validator or BatchCorrectionWorkflowRequestValidator()
        )
        self._design_validator = (
            design_validator or BatchCorrectionWorkflowDesignValidator()
        )
        self._control_site_validator = (
            control_site_validator or BatchCorrectionWorkflowControlSiteValidator()
        )
        self._stage_order_validator = (
            stage_order_validator or BatchCorrectionWorkflowStageOrderValidator()
        )
        self._missingness_validator = (
            missingness_validator or BatchCorrectionWorkflowMissingnessValidator()
        )
        self._interpreter = interpreter or BatchCorrectionPlanInterpreter()
        self._executor = executor or cast(
            BatchCorrectionExecutorContract,
            SpsRuvStyleExecutor(),
        )
        self._provenance_recorder = (
            provenance_recorder or BatchCorrectionProvenanceRecorder()
        )

    def run(self, request: object) -> BatchCorrectionWorkflowResult:
        """Validate, interpret, execute, and record batch correction."""

        validated_request = self._request_validator.run(request)
        dataset_metadata = self._design_validator.run(request=validated_request)
        self._stage_order_validator.run(config=validated_request.config)
        control_site_mapping = self._control_site_validator.run(
            request=validated_request
        )
        missingness_policy = self._missingness_validator.run(request=validated_request)
        plan = self._interpreter.run(
            config=validated_request.config,
            dataset_metadata=dataset_metadata,
            control_site_mapping=control_site_mapping,
            missingness_policy=missingness_policy,
        )
        executor_result = self._executor.run(
            phospho=validated_request.phospho,
            plan=plan,
        )
        provenance = self._provenance_recorder.run(
            request=validated_request,
            dataset_metadata=dataset_metadata,
            control_site_mapping=control_site_mapping,
            missingness_policy=missingness_policy,
            plan=plan,
            executor_result=executor_result,
        )
        diagnostics = cast(
            Mapping[str, JsonValue],
            {
                "interpreter": plan.diagnostic_requirements.to_payload(),
                "executor": executor_result.diagnostics.to_payload(),
            },
        )
        return BatchCorrectionWorkflowResult(
            corrected_matrix=executor_result.corrected_matrix,
            corrected_preprocessing_output=getattr(
                executor_result,
                "corrected_preprocessing_output",
                None,
            ),
            diagnostics=diagnostics,
            warnings=tuple(str(warning) for warning in executor_result.warnings),
            provenance=provenance,
        )


__all__ = ["BatchCorrectionWorkflow"]
