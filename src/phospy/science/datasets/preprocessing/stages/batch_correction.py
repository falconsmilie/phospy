"""Batch-correction stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace

from phospy.contracts.configs import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    DatasetBatchCorrectionConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.batch_correction import (
    BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
    BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS,
    BATCH_CORRECTION_STATUS_DISABLED,
    BatchCorrectionDiagnostics,
    BatchCorrectionEngine,
    BatchCorrectionPolicy,
    BatchCorrectionReport,
)
from phospy.science.datasets.preprocessing.batch_correction_metadata import (
    BatchCorrectionMetadataResolver,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.science.datasets.preprocessing.stage_contract import (
    PreprocessingStageContract,
)
from phospy.validation.datasets.batch_correction import (
    BatchCorrectionAdequacyValidator,
)


class BatchCorrectionStage:
    """Resolve metadata, validate design adequacy, and apply batch correction."""

    stage_key = DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION

    def __init__(
        self,
        *,
        metadata_resolver: BatchCorrectionMetadataResolver | None = None,
        adequacy_validator: BatchCorrectionAdequacyValidator | None = None,
        engine: BatchCorrectionEngine | None = None,
    ) -> None:
        self._metadata_resolver = metadata_resolver or BatchCorrectionMetadataResolver()
        self._adequacy_validator = (
            adequacy_validator or BatchCorrectionAdequacyValidator()
        )
        self._engine = engine or BatchCorrectionEngine()

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        method = _resolve_method(state.plan)
        if method == DATASET_BATCH_CORRECTION_METHOD_NONE:
            report = _build_disabled_report(state)
            return PreprocessingStageResult(
                state=replace(state, batch_correction_report=report),
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {
                        "method": method,
                        "status": BATCH_CORRECTION_STATUS_DISABLED,
                        "matrix_shape_before": list(report.matrix_shape_before or ()),
                        "matrix_shape_after": list(report.matrix_shape_after or ()),
                    },
                },
            )
        if method != DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH:
            raise PhosPyInputError(
                "dataset preprocessing plan contains unsupported "
                f"batch_correction_method={method!r}"
            )

        metadata = self._metadata_resolver.run(
            phospho=state.phospho,
            sample_metadata=state.sample_metadata,
            batch_column=state.plan.batch_correction_batch_column,
            condition_column=state.plan.batch_correction_condition_column,
        )
        self._adequacy_validator.run(
            batch_by_sample=metadata.batch_by_sample,
            condition_by_sample=metadata.condition_by_sample,
            sample_order=metadata.sample_order,
            preserve_condition_effects=(
                state.plan.batch_correction_preserve_condition_effects
            ),
        )
        result = self._engine.run(
            phospho=state.phospho,
            batch_labels=metadata.batch_labels,
            condition_labels=metadata.condition_labels,
            config=DatasetBatchCorrectionConfig(
                method=method,
                batch_column=state.plan.batch_correction_batch_column,
                condition_column=state.plan.batch_correction_condition_column,
                preserve_condition_effects=True,
            ),
        )
        return PreprocessingStageResult(
            state=replace(
                state,
                phospho=result.corrected_matrix,
                batch_correction_metadata=metadata,
                batch_correction_report=result.report,
            ),
            diagnostics={
                "dropped_row_ids": (),
                "dropped_row_count": 0,
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "stage executed",
                "diagnostics": dict(result.diagnostics),
            },
        )


def _resolve_method(plan: PreprocessingPlan) -> str:
    method = str(plan.batch_correction_method).strip()
    if not method:
        return DATASET_BATCH_CORRECTION_METHOD_NONE
    return method


def _build_disabled_report(state: PreprocessingState) -> BatchCorrectionReport:
    shape = (int(state.phospho.shape[0]), int(state.phospho.shape[1]))
    return BatchCorrectionReport(
        status=BATCH_CORRECTION_STATUS_DISABLED,
        policy=BatchCorrectionPolicy(
            method=DATASET_BATCH_CORRECTION_METHOD_NONE,
            batch_column=state.plan.batch_correction_batch_column,
            condition_column=state.plan.batch_correction_condition_column,
            design_preservation_policy=(
                BATCH_CORRECTION_DESIGN_PRESERVATION_PRESERVE_CONDITION_EFFECTS
            ),
            preserve_condition_effects=bool(
                state.plan.batch_correction_preserve_condition_effects
            ),
        ),
        diagnostics=BatchCorrectionDiagnostics(
            confounding_check_status=BATCH_CORRECTION_CONFOUNDING_NOT_APPLICABLE,
            matrix_shape_before=shape,
            matrix_shape_after=shape,
            limitations=("batch correction disabled by preprocessing configuration",),
        ),
    )


def _include_when(plan: PreprocessingPlan) -> bool:
    return _resolve_method(plan) != DATASET_BATCH_CORRECTION_METHOD_NONE


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return _resolve_method(plan)


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    return {
        "method": _resolve_method(plan),
        "batch_column": plan.batch_correction_batch_column,
        "condition_column": plan.batch_correction_condition_column,
        "preserve_condition_effects": bool(
            plan.batch_correction_preserve_condition_effects
        ),
    }


BATCH_CORRECTION_STAGE_CONTRACT = PreprocessingStageContract(
    stage_key=DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    display_label=DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    provenance_stage=DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    operation_name=_resolve_operation,
    serialize_parameters=_resolve_parameters,
    consumed_input_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SAMPLE_METADATA,
    ),
    produced_output_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
    stage_factory=BatchCorrectionStage,
    backend="numpy",
    include_when=_include_when,
    diagnostics_metadata={
        "known_diagnostics_fields": (
            "method",
            "status",
            "number_of_sites",
            "number_of_samples",
            "number_of_batches",
            "batch_levels",
            "condition_levels",
            "condition_design_columns",
            "batch_design_columns",
            "full_design_rank",
            "residual_degrees_of_freedom",
            "matrix_shape_before",
            "matrix_shape_after",
            "max_abs_estimated_batch_contribution",
            "mean_abs_estimated_batch_contribution",
        )
    },
)


__all__ = ["BATCH_CORRECTION_STAGE_CONTRACT", "BatchCorrectionStage"]
