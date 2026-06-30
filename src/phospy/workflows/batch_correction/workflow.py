"""Batch-correction workflow orchestrator."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance import fingerprint_matrix
from phospy.provenance.models import JsonValue
from phospy.science.batch_correction import SpsRuvStyleExecutor
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.validation.workflows.batch_correction.control_site_workflow import (
    BatchCorrectionWorkflowControlSiteValidator,
)
from phospy.validation.workflows.batch_correction.design import (
    BatchCorrectionWorkflowDesignValidator,
    BatchCorrectionWorkflowFactorFeasibilityValidator,
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
    BatchCorrectionExecutorDiagnosticsContract,
    BatchCorrectionExecutorResultContract,
    BatchCorrectionFactorFeasibilityValidatorContract,
    BatchCorrectionInterpreterContract,
    BatchCorrectionMissingnessValidatorContract,
    BatchCorrectionProvenanceRecorderContract,
    BatchCorrectionRequestValidatorContract,
    BatchCorrectionStageOrderValidatorContract,
    BatchCorrectionWorkflowRequest,
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
        factor_feasibility_validator: (
            BatchCorrectionFactorFeasibilityValidatorContract | None
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
        self._factor_feasibility_validator = (
            factor_feasibility_validator
            or BatchCorrectionWorkflowFactorFeasibilityValidator()
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
        self._factor_feasibility_validator.run(
            request=validated_request,
            dataset_metadata=dataset_metadata,
            control_site_mapping=control_site_mapping,
            missingness_policy=missingness_policy,
        )
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
        executor_result = _with_conservative_upstream_observation_mask(
            request=validated_request,
            executor_result=executor_result,
        )
        corrected_preprocessing_output = getattr(
            executor_result,
            "corrected_preprocessing_output",
            None,
        )
        if not isinstance(corrected_preprocessing_output, CorrectedPreprocessingOutput):
            raise PhosPyInputError(
                "batch-correction workflow did not produce complete "
                "analysis-ready output; corrected_preprocessing_output is "
                "unavailable because the corrected matrix is diagnostic-only "
                "or contains restored missing cells"
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
        corrected_preprocessing_output = corrected_preprocessing_output.with_provenance(
            provenance
        )
        return BatchCorrectionWorkflowResult(
            corrected_matrix=corrected_preprocessing_output.corrected_matrix,
            corrected_preprocessing_output=corrected_preprocessing_output,
            diagnostics=diagnostics,
            warnings=tuple(str(warning) for warning in executor_result.warnings),
            provenance=provenance,
        )


@dataclass(frozen=True, slots=True)
class _ExecutorResultWithCombinedMask:
    inner: BatchCorrectionExecutorResultContract
    upstream_observation_mask: pd.DataFrame
    output_observation_mask: pd.DataFrame
    corrected_cell_status: pd.DataFrame
    corrected_preprocessing_output: CorrectedPreprocessingOutput | None

    @property
    def corrected_matrix(self) -> pd.DataFrame:
        return self.inner.corrected_matrix

    @property
    def diagnostics(self) -> BatchCorrectionExecutorDiagnosticsContract:
        return self.inner.diagnostics

    @property
    def warnings(self) -> tuple[str, ...]:
        return tuple(str(warning) for warning in self.inner.warnings)

    @property
    def executor_output_observation_mask(self) -> pd.DataFrame:
        return self.inner.output_observation_mask

    @property
    def provenance_payload(self) -> Mapping[str, object]:
        payload = dict(self.inner.provenance_payload)
        payload["observation_mask_lineage"] = {
            "upstream_observation_mask_fingerprint": _fingerprint_payload(
                self.upstream_observation_mask,
                name="batch_correction.workflow.upstream_observation_mask",
            ),
            "executor_output_observation_mask_fingerprint": _fingerprint_payload(
                self.inner.output_observation_mask,
                name="batch_correction.workflow.executor_output_observation_mask",
            ),
            "final_combined_observation_mask_fingerprint": _fingerprint_payload(
                self.output_observation_mask,
                name="batch_correction.workflow.final_combined_observation_mask",
            ),
            "combination_rule": (
                "final_combined_observation_mask = "
                "upstream_observation_mask & executor_output_observation_mask"
            ),
            "final_observation_mask_source": "combined_upstream_and_executor_masks",
        }
        return payload

    @property
    def rejected_rows(self) -> tuple[str, ...]:
        return tuple(str(row) for row in self.inner.rejected_rows)

    @property
    def rejected_cells(self) -> tuple[tuple[str, str], ...]:
        return tuple(
            (str(row), str(column)) for row, column in self.inner.rejected_cells
        )

    @property
    def withheld_rows(self) -> tuple[str, ...]:
        return tuple(str(row) for row in self.inner.withheld_rows)

    @property
    def withheld_cells(self) -> tuple[tuple[str, str], ...]:
        return tuple(
            (str(row), str(column)) for row, column in self.inner.withheld_cells
        )


def _with_conservative_upstream_observation_mask(
    *,
    request: object,
    executor_result: BatchCorrectionExecutorResultContract,
) -> BatchCorrectionExecutorResultContract:
    if not isinstance(request, BatchCorrectionWorkflowRequest):
        return executor_result
    upstream_mask = request.upstream_observation_mask
    if upstream_mask is None:
        return executor_result

    combined_mask = _combine_observation_masks(
        upstream_mask=upstream_mask,
        correction_mask=executor_result.output_observation_mask,
        corrected_matrix=executor_result.corrected_matrix,
    )
    corrected_status = _corrected_cell_status(
        observed_mask=combined_mask,
        existing=getattr(executor_result, "corrected_cell_status", None),
    )
    corrected_output = getattr(executor_result, "corrected_preprocessing_output", None)
    if isinstance(corrected_output, CorrectedPreprocessingOutput):
        corrected_output = replace(
            corrected_output,
            output_observation_mask=combined_mask,
            corrected_cell_status=corrected_status,
        )
    else:
        corrected_output = None
    return _ExecutorResultWithCombinedMask(
        inner=executor_result,
        upstream_observation_mask=upstream_mask,
        output_observation_mask=combined_mask,
        corrected_cell_status=corrected_status,
        corrected_preprocessing_output=corrected_output,
    )


def _combine_observation_masks(
    *,
    upstream_mask: pd.DataFrame,
    correction_mask: pd.DataFrame,
    corrected_matrix: pd.DataFrame,
) -> pd.DataFrame:
    _require_aligned_boolean_mask(
        mask=upstream_mask,
        expected=corrected_matrix,
        field_name="upstream observation mask",
    )
    _require_aligned_boolean_mask(
        mask=correction_mask,
        expected=corrected_matrix,
        field_name="correction output observation mask",
    )
    # Validation above guarantees boolean-aligned masks; combine the raw arrays
    # once to avoid repeated full-frame defensive copies on large imputation masks.
    combined = np.logical_and(
        upstream_mask.to_numpy(dtype=bool, copy=False),
        correction_mask.to_numpy(dtype=bool, copy=False),
    )
    return pd.DataFrame(
        combined,
        index=corrected_matrix.index.copy(),
        columns=corrected_matrix.columns.copy(),
    )


def _require_aligned_boolean_mask(
    *,
    mask: pd.DataFrame,
    expected: pd.DataFrame,
    field_name: str,
) -> None:
    if not mask.index.equals(expected.index):
        raise PhosPyInputError(
            f"batch-correction workflow {field_name} alignment failed: "
            "mask index must match corrected matrix index"
        )
    if not mask.columns.equals(expected.columns):
        raise PhosPyInputError(
            f"batch-correction workflow {field_name} alignment failed: "
            "mask columns must match corrected matrix columns"
        )
    values = mask.to_numpy(dtype="object", copy=False)
    missing_values: npt.NDArray[np.bool_] = np.asarray(
        pd.isna(values),
        dtype=bool,
    )
    if bool(missing_values.any()):
        raise PhosPyInputError(
            f"batch-correction workflow {field_name} must contain only "
            "boolean values; missing values are not allowed"
        )
    if all(pd.api.types.is_bool_dtype(dtype) for dtype in mask.dtypes):
        return

    boolean_result = np.frompyfunc(
        lambda value: isinstance(value, (bool, np.bool_)),
        1,
        1,
    )(values)
    boolean_cells: npt.NDArray[np.bool_] = np.asarray(boolean_result, dtype=bool)
    invalid_locations = np.argwhere(~boolean_cells)
    if invalid_locations.size == 0:
        return
    row_position, column_position = invalid_locations[0]
    row_id = mask.index[int(row_position)]
    column_id = mask.columns[int(column_position)]
    raise PhosPyInputError(
        f"batch-correction workflow {field_name} must contain only "
        f"boolean values; invalid value at "
        f"({str(row_id)!r}, {str(column_id)!r})"
    )


def _corrected_cell_status(
    *,
    observed_mask: pd.DataFrame,
    existing: object,
) -> pd.DataFrame:
    if (
        isinstance(existing, pd.DataFrame)
        and existing.index.equals(observed_mask.index)
        and existing.columns.equals(observed_mask.columns)
    ):
        status = existing.copy(deep=True)
    else:
        status = pd.DataFrame(
            "corrected_observed",
            index=observed_mask.index.copy(),
            columns=observed_mask.columns.copy(),
        )
    return status.mask(~observed_mask.astype(bool), "restored_missing")


def _fingerprint_payload(mask: pd.DataFrame, *, name: str) -> dict[str, object]:
    fingerprint = fingerprint_matrix(mask.astype("int8"), name=name)
    return {
        "name": fingerprint.name,
        "rows": int(fingerprint.rows),
        "columns": int(fingerprint.columns),
        "index_name": fingerprint.index_name,
        "column_names": list(fingerprint.column_names),
        "dtypes": list(fingerprint.dtypes),
        "exact_hash_algorithm": fingerprint.exact_hash_algorithm,
        "exact_hash_value": fingerprint.exact_hash_value,
        "tolerance_hash_algorithm": fingerprint.tolerance_hash_algorithm,
        "tolerance_hash_value": fingerprint.tolerance_hash_value,
        "index_structure": fingerprint.index_structure,
        "column_index_structure": fingerprint.column_index_structure,
    }


__all__ = ["BatchCorrectionWorkflow"]
