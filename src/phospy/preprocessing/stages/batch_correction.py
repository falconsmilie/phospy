"""Thin core adapter and contract for batch-correction preprocessing."""

from __future__ import annotations

from dataclasses import replace

from phospy.errors.input import PhosPyInputError
from phospy.preprocessing.batch_correction.contracts.configuration import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    SPS_RUV_BATCH_CORRECTION_METHODS,
)
from phospy.preprocessing.batch_correction.contracts.request import (
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
    SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER,
)
from phospy.preprocessing.batch_correction.contracts.stage import (
    BatchCorrectionAdequacyValidatorProtocol,
    BatchDesignMetadataValidatorProtocol,
    SpsRuvStyleBatchCorrectionRunner,
)
from phospy.preprocessing.batch_correction.execution.linear import (
    LinearResidualizeBatchCorrectionEngine,
)
from phospy.preprocessing.batch_correction.execution.stage import (
    BatchCorrectionStageExecutionRequest,
    BatchCorrectionStageExecutor,
    resolve_batch_correction_method,
)
from phospy.preprocessing.batch_correction.metadata import (
    BatchCorrectionMetadataResolver,
)
from phospy.preprocessing.contracts import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.preprocessing.contracts.stage import (
    DeterminismKind,
    PreprocessingStageContract,
    PreprocessingStageFactoryContext,
)
from phospy.quantitative.models import IntensityScaleKind, QuantitativeMeaning
from phospy.quantitative.quantitative_contracts import (
    NegativeDomainPolicy,
    QuantitativeEvidenceRequirement,
    QuantitativeInformationLossKind,
    QuantitativeOperationContract,
    QuantitativeReversibilityKind,
    preserve_meaning_transition,
    preserve_quantitative_contract,
    preserve_scale_transition,
)


class BatchCorrectionStage:
    """Adapt core pipeline state to capability-owned correction execution."""

    stage_key = DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION

    def __init__(
        self,
        *,
        metadata_resolver: BatchCorrectionMetadataResolver | None = None,
        metadata_validator: BatchDesignMetadataValidatorProtocol | None = None,
        adequacy_validator: BatchCorrectionAdequacyValidatorProtocol | None = None,
        engine: LinearResidualizeBatchCorrectionEngine | None = None,
        sps_ruv_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
    ) -> None:
        self._executor = BatchCorrectionStageExecutor(
            metadata_resolver=metadata_resolver,
            metadata_validator=metadata_validator,
            adequacy_validator=adequacy_validator,
            engine=engine,
            sps_ruv_runner=sps_ruv_runner,
        )
        # Preserve the established collaborator inspection seams on the thin
        # core adapter while execution ownership resides in the capability.
        self._metadata_validator = self._executor._metadata_validator
        self._adequacy_validator = self._executor._adequacy_validator

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        result = self._executor.run(_execution_request(state))
        return PreprocessingStageResult(
            state=replace(
                state,
                phospho=result.corrected_matrix,
                imputation_observation_mask=result.output_observation_mask,
                batch_correction_metadata=(
                    result.metadata
                    if result.metadata is not None
                    else state.batch_correction_metadata
                ),
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
            batch_correction_provenance=result.provenance,
        )

    def validate_before_quantitative_contract(self, state: PreprocessingState) -> None:
        self._executor.validate_linear(_execution_request(state))


def _execution_request(
    state: PreprocessingState,
) -> BatchCorrectionStageExecutionRequest:
    plan = state.plan
    return BatchCorrectionStageExecutionRequest(
        phospho=state.phospho,
        sample_metadata=state.sample_metadata,
        site_metadata=state.site_metadata,
        batch_correction_method=plan.batch_correction_method,
        batch_correction_batch_column=plan.batch_correction_batch_column,
        batch_correction_condition_column=plan.batch_correction_condition_column,
        batch_correction_condition_columns=plan.batch_correction_condition_columns,
        batch_correction_preserve_condition_effects=(
            plan.batch_correction_preserve_condition_effects
        ),
        batch_correction_internal_request=plan.batch_correction_internal_request,
        batch_correction_control_site_set=plan.batch_correction_control_site_set,
        batch_correction_missingness_policy=plan.batch_correction_missingness_policy,
        imputation_observation_mask=state.imputation_observation_mask,
        stage_order=plan.stage_order,
    )


def _include_when(plan: PreprocessingPlan) -> bool:
    return (
        resolve_batch_correction_method(plan.batch_correction_method)
        != DATASET_BATCH_CORRECTION_METHOD_NONE
    )


def _resolve_operation(plan: PreprocessingPlan) -> str:
    return resolve_batch_correction_method(plan.batch_correction_method)


def _resolve_parameters(plan: PreprocessingPlan) -> dict[str, object]:
    parameters: dict[str, object] = {
        "method": resolve_batch_correction_method(plan.batch_correction_method),
        "batch_column": plan.batch_correction_batch_column,
        "condition_column": plan.batch_correction_condition_column,
        "condition_columns": list(plan.batch_correction_condition_columns),
        "preserve_condition_effects": bool(
            plan.batch_correction_preserve_condition_effects
        ),
    }
    if (
        resolve_batch_correction_method(plan.batch_correction_method)
        in SPS_RUV_BATCH_CORRECTION_METHODS
    ):
        request = plan.batch_correction_internal_request
        executed_stage_order = (
            None
            if request is None
            or request.stage_order
            is not SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER
            else list(SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER)
        )
        parameters.update(
            {
                "replicate_column": plan.batch_correction_replicate_column,
                "n_unwanted_factors": (
                    None if request is None else request.n_unwanted_factors
                ),
                "diagnostics_enabled": (
                    None if request is None else request.diagnostics_enabled
                ),
                "requested_stage_order": (
                    None if request is None else request.stage_order.value
                ),
                "executed_stage_order": executed_stage_order,
                "stage_order_policy": (
                    None if request is None else request.stage_order.value
                ),
                "control_site_source": (
                    None if request is None else request.control_site_source.value
                ),
                "missing_value_policy": (
                    None if request is None else request.missing_value_policy.value
                ),
                "imputation_policy": (
                    None if request is None else request.imputation_policy.value
                ),
            }
        )
    return parameters


def _resolve_quantitative_contract(
    plan: PreprocessingPlan,
) -> QuantitativeOperationContract:
    method = resolve_batch_correction_method(plan.batch_correction_method)
    if method == DATASET_BATCH_CORRECTION_METHOD_NONE:
        return preserve_quantitative_contract(
            required_evidence=frozenset({QuantitativeEvidenceRequirement.NONE}),
            negative_domain_policy=NegativeDomainPolicy.PRESERVES_INPUT_DOMAIN,
            reversibility=QuantitativeReversibilityKind.REVERSIBLE,
            information_loss=QuantitativeInformationLossKind.NONE,
        )
    if (
        method == DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH
        or method in SPS_RUV_BATCH_CORRECTION_METHODS
    ):
        accepted_meanings = frozenset(
            {
                QuantitativeMeaning.PHOSPHOSITE_LOG_ABUNDANCE,
                QuantitativeMeaning.UNKNOWN,
            }
        )
        required_evidence = {
            QuantitativeEvidenceRequirement.SAMPLE_METADATA_DESIGN,
        }
        if method in SPS_RUV_BATCH_CORRECTION_METHODS:
            required_evidence.add(QuantitativeEvidenceRequirement.CONTROL_SITE_SET)
            required_evidence.add(QuantitativeEvidenceRequirement.MISSINGNESS_MASK)
        return QuantitativeOperationContract(
            accepted_input_scale_kinds=frozenset({IntensityScaleKind.LOG2}),
            accepted_quantitative_meanings=accepted_meanings,
            output_scale_transition=preserve_scale_transition(
                frozenset({IntensityScaleKind.LOG2}),
                output_scale_label="log2",
            ),
            output_meaning_transition=preserve_meaning_transition(accepted_meanings),
            preserves_abundance=False,
            negative_domain_policy=NegativeDomainPolicy.MAY_INTRODUCE_NEGATIVE_VALUES,
            required_evidence=frozenset(required_evidence),
            reversibility=QuantitativeReversibilityKind.IRREVERSIBLE,
            information_loss=QuantitativeInformationLossKind.ADDITIVE_RESIDUALIZATION,
        )
    raise PhosPyInputError(
        "dataset preprocessing plan contains unsupported "
        f"batch_correction_method={method!r}"
    )


def _build_batch_correction_stage(
    context: PreprocessingStageFactoryContext,
) -> BatchCorrectionStage:
    return BatchCorrectionStage(
        sps_ruv_runner=context.batch_correction_runner,
        metadata_validator=context.batch_correction_metadata_validator,
        adequacy_validator=context.batch_correction_adequacy_validator,
    )


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
    produced_output_tables=(
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_IMPUTATION_OBSERVATION_MASK,
    ),
    quantitative_contract=_resolve_quantitative_contract,
    stage_factory=_build_batch_correction_stage,
    backend="numpy",
    determinism_kind=DeterminismKind.DETERMINISTIC,
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


__all__ = [
    "BATCH_CORRECTION_STAGE_CONTRACT",
    "BatchCorrectionStage",
]
