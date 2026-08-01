"""Batch-correction stage for dataset preprocessing."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace
from typing import Protocol

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.provenance.models import BatchCorrectionProvenance, JsonValue
from phospy.science.configs import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
    SPS_RUV_BATCH_CORRECTION_METHODS,
    SUPPORTED_INTERNAL_BATCH_CORRECTION_EXECUTED_STAGE_ORDER,
    SUPPORTED_INTERNAL_BATCH_CORRECTION_STAGE_ORDER,
    DatasetBatchCorrectionConfig,
)
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
    ResolvedBatchCorrectionMetadata,
)
from phospy.science.datasets.preprocessing.batch_correction_provenance import (
    build_native_batch_correction_provenance,
)
from phospy.science.datasets.preprocessing.control_sites import ControlSiteSet
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
)
from phospy.science.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION,
    PreprocessingPlan,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.science.datasets.preprocessing.stage_contract import (
    DeterminismKind,
    PreprocessingStageContract,
    PreprocessingStageFactoryContext,
)
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
)
from phospy.science.transformations.quantitative_contracts import (
    NegativeDomainPolicy,
    QuantitativeEvidenceRequirement,
    QuantitativeInformationLossKind,
    QuantitativeOperationContract,
    QuantitativeReversibilityKind,
    preserve_meaning_transition,
    preserve_quantitative_contract,
    preserve_scale_transition,
)


class BatchDesignMetadataValidatorProtocol(Protocol):
    """Validator required before resolving batch-correction design metadata."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        batch_column: str,
        condition_columns: tuple[str, ...],
        context: str,
    ) -> object: ...


class BatchCorrectionAdequacyValidatorProtocol(Protocol):
    """Validator required before executing linear residual batch correction."""

    def run(
        self,
        *,
        batch_by_sample: Mapping[str, str],
        condition_by_sample: Mapping[str, str],
        sample_order: tuple[str, ...],
        preserve_condition_effects: bool,
    ) -> object: ...


class _MissingBatchDesignMetadataValidator:
    def run(self, **_: object) -> object:
        raise PhosPyInputError(
            "BatchCorrectionStage requires an injected BatchDesignMetadataValidator"
        )


class _MissingBatchCorrectionAdequacyValidator:
    def run(self, **_: object) -> object:
        raise PhosPyInputError(
            "BatchCorrectionStage requires an injected BatchCorrectionAdequacyValidator"
        )


class SpsRuvStyleBatchCorrectionResult(Protocol):
    """Result shape required by the preprocessing batch-correction stage."""

    @property
    def corrected_preprocessing_output(self) -> CorrectedPreprocessingOutput: ...

    @property
    def diagnostics(self) -> Mapping[str, JsonValue]: ...

    @property
    def provenance(self) -> BatchCorrectionProvenance: ...


class SpsRuvStyleBatchCorrectionRunner(Protocol):
    """Runner protocol for SPS/RUV-style batch correction."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        config: object,
        sample_metadata: pd.DataFrame | None,
        control_site_set: ControlSiteSet,
        missingness_policy: object,
        upstream_observation_mask: pd.DataFrame | None,
        site_metadata: pd.DataFrame,
    ) -> SpsRuvStyleBatchCorrectionResult: ...


class BatchCorrectionStage:
    """Resolve metadata, validate design adequacy, and apply batch correction."""

    stage_key = DATASET_PREPROCESSING_STAGE_BATCH_CORRECTION

    def __init__(
        self,
        *,
        metadata_resolver: BatchCorrectionMetadataResolver | None = None,
        metadata_validator: BatchDesignMetadataValidatorProtocol | None = None,
        adequacy_validator: BatchCorrectionAdequacyValidatorProtocol | None = None,
        engine: BatchCorrectionEngine | None = None,
        sps_ruv_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
    ) -> None:
        self._metadata_resolver = metadata_resolver or BatchCorrectionMetadataResolver()
        self._metadata_validator = (
            metadata_validator or _MissingBatchDesignMetadataValidator()
        )
        self._adequacy_validator = (
            adequacy_validator or _MissingBatchCorrectionAdequacyValidator()
        )
        self._engine = engine or BatchCorrectionEngine()
        self._sps_ruv_runner = sps_ruv_runner

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
        if method in SPS_RUV_BATCH_CORRECTION_METHODS:
            return _run_sps_ruv_style_correction(
                state,
                runner=self._sps_ruv_runner,
            )
        if method != DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH:
            raise PhosPyInputError(
                "dataset preprocessing plan contains unsupported "
                f"batch_correction_method={method!r}"
            )

        try:
            metadata = self._resolve_validated_linear_residualize_metadata(state)
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
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                _validation_failure_message(state=state, method=method, reason=str(exc))
            ) from exc
        provenance = build_native_batch_correction_provenance(
            input_matrix=state.phospho,
            output_matrix=result.corrected_matrix,
            plan=state.plan,
            report=result.report,
            metadata=metadata,
            diagnostics=result.diagnostics,
            warnings=result.report.warnings,
            observation_mask=state.imputation_observation_mask,
            source="native_preprocessing_stage",
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
            batch_correction_provenance=provenance,
        )

    def validate_before_quantitative_contract(self, state: PreprocessingState) -> None:
        method = _resolve_method(state.plan)
        if method != DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH:
            return
        try:
            self._resolve_validated_linear_residualize_metadata(state)
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                _validation_failure_message(state=state, method=method, reason=str(exc))
            ) from exc

    def _resolve_validated_linear_residualize_metadata(
        self,
        state: PreprocessingState,
    ) -> ResolvedBatchCorrectionMetadata:
        self._metadata_validator.run(
            phospho=state.phospho,
            sample_metadata=state.sample_metadata,
            batch_column=state.plan.batch_correction_batch_column,
            condition_columns=(state.plan.batch_correction_condition_column,),
            context="dataset build request preprocessing_config.batch_correction",
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
        return metadata


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
            condition_columns=state.plan.batch_correction_condition_columns,
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
    parameters: dict[str, object] = {
        "method": _resolve_method(plan),
        "batch_column": plan.batch_correction_batch_column,
        "condition_column": plan.batch_correction_condition_column,
        "condition_columns": list(plan.batch_correction_condition_columns),
        "preserve_condition_effects": bool(
            plan.batch_correction_preserve_condition_effects
        ),
    }
    if _resolve_method(plan) in SPS_RUV_BATCH_CORRECTION_METHODS:
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
    method = _resolve_method(plan)
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


def _run_sps_ruv_style_correction(
    state: PreprocessingState,
    *,
    runner: SpsRuvStyleBatchCorrectionRunner | None,
) -> PreprocessingStageResult:
    request = state.plan.batch_correction_internal_request
    control_site_set = state.plan.batch_correction_control_site_set
    missingness_policy = state.plan.batch_correction_missingness_policy
    if request is None or control_site_set is None or missingness_policy is None:
        raise PhosPyInputError(
            "dataset preprocessing plan contains incomplete SPS/RUV-style "
            "batch-correction configuration"
        )
    if not isinstance(control_site_set, ControlSiteSet):
        raise PhosPyInputError(
            "dataset preprocessing plan SPS/RUV-style batch correction requires "
            "a ControlSiteSet"
        )
    if runner is None:
        raise PhosPyInputError(
            "dataset preprocessing plan SPS/RUV-style batch correction requires "
            "an injected SpsRuvStyleBatchCorrectionRunner"
        )
    try:
        result = runner.run(
            phospho=state.phospho,
            config=request,
            sample_metadata=state.sample_metadata,
            control_site_set=control_site_set,
            missingness_policy=missingness_policy,
            upstream_observation_mask=state.imputation_observation_mask,
            site_metadata=state.site_metadata,
        )
    except PhosPyInputError as exc:
        raise PhosPyInputError(
            _validation_failure_message(
                state=state,
                method=_resolve_method(state.plan),
                reason=str(exc),
            )
        ) from exc

    corrected_output = result.corrected_preprocessing_output
    if not isinstance(corrected_output, CorrectedPreprocessingOutput):
        raise PhosPyInputError(
            "SPS/RUV-style batch correction did not produce complete "
            "analysis-ready preprocessing output; check missingness policy and "
            "observation-mask diagnostics"
        )
    diagnostics = dict(result.diagnostics)
    return PreprocessingStageResult(
        state=replace(
            state,
            phospho=corrected_output.corrected_matrix,
            imputation_observation_mask=corrected_output.output_observation_mask,
            batch_correction_report=corrected_output.batch_correction_report,
        ),
        diagnostics={
            "dropped_row_ids": (),
            "dropped_row_count": 0,
            "imputed_cell_count": 0,
            "imputed_row_ids": (),
            "notes": "stage executed",
            "diagnostics": diagnostics,
        },
        batch_correction_provenance=result.provenance,
    )


def _validation_failure_message(
    *,
    state: PreprocessingState,
    method: str,
    reason: str,
) -> str:
    return (
        "batch correction validation failed before correction execution; "
        f"method={method!r}; "
        f"input_shape={[int(state.phospho.shape[0]), int(state.phospho.shape[1])]}; "
        f"batch_column={state.plan.batch_correction_batch_column!r}; "
        f"condition_column={state.plan.batch_correction_condition_column!r}; "
        f"stage_order={list(state.plan.stage_order)!r}; "
        f"reason: {reason}"
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
    produced_output_tables=(PreprocessingStateTableKey.DATASET_PHOSPHO,),
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
    "SpsRuvStyleBatchCorrectionResult",
    "SpsRuvStyleBatchCorrectionRunner",
]
