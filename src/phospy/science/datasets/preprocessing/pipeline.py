"""Internal dataset preprocessing pipeline orchestration."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace

from phospy._deprecations import warn_deprecated
from phospy.errors.build import DatasetBuildError
from phospy.science.datasets.preprocessing.diagnostics_normalization import (
    _StageDiagnosticsDefaultsResolver,
    _StageDiagnosticsNormalizer,
)
from phospy.science.datasets.preprocessing.event_validation import (
    _TransformationEventSequenceValidator,
)
from phospy.science.datasets.preprocessing.fingerprints import (
    _resolve_state_table as _resolve_state_table,
)
from phospy.science.datasets.preprocessing.fingerprints import (
    _StageFingerprintService,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingReportRow,
    PreprocessingStage,
    PreprocessingStageExecution,
    PreprocessingStageResult,
    PreprocessingState,
    PreprocessingStateTableKey,
)
from phospy.science.datasets.preprocessing.quantitative_evidence import (
    QuantitativeOperationEvidenceContext,
    QuantitativeOperationEvidenceValidator,
)
from phospy.science.datasets.preprocessing.report_rows import (
    validate_preprocessing_report_row,
)
from phospy.science.datasets.preprocessing.stage_contract import (
    InterpretedPreprocessingStageContract,
    PreprocessingStageFactoryContext,
    validate_preprocessing_stage_instance,
)
from phospy.science.datasets.preprocessing.stage_registry import (
    PreprocessingStageMetadata,
    build_registered_preprocessing_stage_instances,
    resolve_registered_preprocessing_stages,
)
from phospy.science.datasets.preprocessing.stages import (
    BatchCorrectionAdequacyValidatorProtocol,
    BatchDesignMetadataValidatorProtocol,
    SpsRuvStyleBatchCorrectionRunner,
)
from phospy.science.datasets.preprocessing.trace_builder import _StageTraceBuilder
from phospy.science.transformations.models import (
    IntensityScaleKind,
    QuantitativeMeaning,
    QuantitativeMeaningEvidenceMode,
)
from phospy.science.transformations.quantitative_contracts import (
    QuantitativeContractState,
    initial_quantitative_contract_state,
)


class PreprocessingPipeline:
    """Apply ordered preprocessing stages for interpreted dataset input."""

    def __init__(
        self,
        *,
        stage_registry: tuple[PreprocessingStage, ...] | None = None,
        stage_contract_registry: tuple[PreprocessingStageMetadata, ...] | None = None,
        stage_metadata_registry: tuple[PreprocessingStageMetadata, ...] | None = None,
        batch_correction_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
        batch_correction_metadata_validator: (
            BatchDesignMetadataValidatorProtocol | None
        ) = None,
        batch_correction_adequacy_validator: (
            BatchCorrectionAdequacyValidatorProtocol | None
        ) = None,
    ) -> None:
        contract_overrides = _resolve_stage_contract_overrides(
            stage_contract_registry=stage_contract_registry,
            stage_metadata_registry=stage_metadata_registry,
        )
        resolved_metadata_registry = resolve_registered_preprocessing_stages(
            contract_overrides
        )
        factory_context = PreprocessingStageFactoryContext(
            batch_correction_runner=batch_correction_runner,
            batch_correction_metadata_validator=batch_correction_metadata_validator,
            batch_correction_adequacy_validator=batch_correction_adequacy_validator,
        )
        stages = stage_registry or build_registered_preprocessing_stage_instances(
            resolved_metadata_registry,
            context=factory_context,
        )
        validated_stages = tuple(
            validate_preprocessing_stage_instance(
                stage,
                context="dataset preprocessing stage registry",
            )
            for stage in stages
        )
        self._stage_contract_by_key = {
            metadata.stage_key: metadata for metadata in resolved_metadata_registry
        }
        self._stages_by_key = {stage.stage_key: stage for stage in validated_stages}
        if len(self._stages_by_key) != len(validated_stages):
            raise DatasetBuildError(
                "dataset preprocessing stage registry contains duplicate stage keys"
            )
        for stage in validated_stages:
            if stage.stage_key not in self._stage_contract_by_key:
                raise DatasetBuildError(
                    "dataset preprocessing stage registry contains stage "
                    f"{stage.stage_key!r} without registered stage contract"
                )

    def run(
        self,
        state: PreprocessingState,
        *,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> PreprocessingState:
        final_state, _ = self.run_with_trace(
            state,
            initial_quantitative_scale_kind=initial_quantitative_scale_kind,
            initial_quantitative_meaning=initial_quantitative_meaning,
        )
        return final_state

    def run_with_trace(
        self,
        state: PreprocessingState,
        *,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> tuple[PreprocessingState, tuple[PreprocessingStageExecution, ...]]:
        current = state
        trace: list[PreprocessingStageExecution] = []
        report_rows: list[PreprocessingReportRow] = list(current.report_rows)
        quantitative_state = initial_quantitative_contract_state(
            declared_input_scale_kind=initial_quantitative_scale_kind,
            explicit_quantitative_meaning=initial_quantitative_meaning,
        )
        event_validator = _TransformationEventSequenceValidator()
        diagnostics_defaults_resolver = _StageDiagnosticsDefaultsResolver()
        diagnostics_normalizer = _StageDiagnosticsNormalizer()
        fingerprint_service = _StageFingerprintService()
        trace_builder = _StageTraceBuilder()
        evidence_validator = QuantitativeOperationEvidenceValidator()
        quantitative_meaning_evidence_mode = (
            _initial_quantitative_meaning_evidence_mode(
                initial_quantitative_meaning=initial_quantitative_meaning
            )
        )

        interpreted_stage_contracts = self._interpret_stage_contracts_for_plan(
            current.plan
        )
        for stage_key, contract, interpreted_contract in interpreted_stage_contracts:
            stage = self._resolve_stage(stage_key)
            previous = current
            stage_input_quantitative_state = quantitative_state
            _run_stage_pre_quantitative_contract_validation(
                stage=stage,
                state=previous,
            )
            _validate_quantitative_contract_before_execution(
                quantitative_state=quantitative_state,
                stage_key=stage_key,
                interpreted_contract=interpreted_contract,
            )
            stage_result = stage.run(current)
            if not isinstance(stage_result, PreprocessingStageResult):
                raise DatasetBuildError(
                    "dataset preprocessing stage returned an invalid result payload: "
                    f"{stage_key}"
                )

            intensity_transformation_event = event_validator.run(
                stage_key=stage_key,
                event=stage_result.intensity_transformation_event,
            )
            current = stage_result.state
            report_rows.extend(_normalize_report_rows(stage_result.report_rows))
            diagnostics = diagnostics_normalizer.run(
                stage_key=stage_key,
                raw=stage_result.diagnostics,
                defaults=diagnostics_defaults_resolver.run(
                    previous=previous,
                    current=current,
                ),
            )
            fingerprints = fingerprint_service.run(
                stage_key=stage_key,
                previous=previous,
                current=current,
                consumed_input_tables=interpreted_contract.consumed_input_tables,
                produced_output_tables=interpreted_contract.produced_output_tables,
            )
            trace_record = trace_builder.run(
                stage_key=stage_key,
                contract=contract,
                interpreted_contract=interpreted_contract,
                previous=previous,
                current=current,
                stage_result=stage_result,
                diagnostics=diagnostics,
                fingerprints=fingerprints,
                intensity_transformation_event=intensity_transformation_event,
            )
            evidence_validator.validate(
                QuantitativeOperationEvidenceContext(
                    stage=trace_record.stage,
                    operation=trace_record.operation,
                    quantitative_contract=interpreted_contract.quantitative_contract,
                    trace_record=trace_record,
                    input_quantitative_state=stage_input_quantitative_state,
                    input_quantitative_meaning_evidence_mode=(
                        quantitative_meaning_evidence_mode
                    ),
                    interpreted_contract=interpreted_contract,
                    previous_preprocessing_state=previous,
                    current_preprocessing_state=current,
                )
            )
            quantitative_state = (
                interpreted_contract.quantitative_contract.validate_and_transition(
                    stage_input_quantitative_state,
                    stage=stage_key,
                    operation=interpreted_contract.operation,
                    evidence=trace_record.quantitative_transition_evidence,
                )
            )
            quantitative_meaning_evidence_mode = (
                _resolve_output_quantitative_meaning_evidence_mode(
                    input_state=stage_input_quantitative_state,
                    output_state=quantitative_state,
                    input_mode=quantitative_meaning_evidence_mode,
                )
            )
            trace.append(trace_record)
        if report_rows:
            current = replace(current, report_rows=tuple(report_rows))
        return current, tuple(trace)

    def validate_quantitative_contracts(
        self,
        *,
        plan: object,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> QuantitativeContractState:
        if not isinstance(plan, PreprocessingPlan):
            raise DatasetBuildError(
                "dataset preprocessing quantitative contract validation requires "
                "a PreprocessingPlan"
            )
        quantitative_state = initial_quantitative_contract_state(
            declared_input_scale_kind=initial_quantitative_scale_kind,
            explicit_quantitative_meaning=initial_quantitative_meaning,
        )
        for (
            stage_key,
            _contract,
            interpreted_contract,
        ) in self._interpret_stage_contracts_for_plan(plan):
            quantitative_state = _validate_quantitative_contract_before_execution(
                quantitative_state=quantitative_state,
                stage_key=stage_key,
                interpreted_contract=interpreted_contract,
            )
        return quantitative_state

    def _resolve_stage(self, stage_key: str) -> PreprocessingStage:
        stage = self._stages_by_key.get(stage_key)
        if stage is not None:
            return stage
        raise DatasetBuildError(
            f"dataset preprocessing plan references an unsupported stage: {stage_key}"
        )

    def _resolve_stage_contract(self, stage_key: str) -> PreprocessingStageMetadata:
        contract = self._stage_contract_by_key.get(stage_key)
        if contract is not None:
            return contract
        raise DatasetBuildError(
            "dataset preprocessing stage metadata is not registered for "
            f"stage {stage_key!r}"
        )

    def _interpret_stage_contracts_for_plan(
        self,
        plan: PreprocessingPlan,
    ) -> tuple[
        tuple[str, PreprocessingStageMetadata, InterpretedPreprocessingStageContract],
        ...,
    ]:
        available_tables = set(_INITIAL_PREPROCESSING_STATE_TABLES)
        interpreted_stage_contracts: list[
            tuple[
                str,
                PreprocessingStageMetadata,
                InterpretedPreprocessingStageContract,
            ]
        ] = []
        for stage_key in plan.stage_order:
            contract = self._resolve_stage_contract(stage_key)
            interpreted_contract = contract.interpret(plan)
            _validate_declared_table_dependencies(
                stage_key=stage_key,
                consumed_input_tables=interpreted_contract.consumed_input_tables,
                available_tables=available_tables,
            )
            available_tables.update(interpreted_contract.produced_output_tables)
            interpreted_stage_contracts.append(
                (stage_key, contract, interpreted_contract)
            )
        return tuple(interpreted_stage_contracts)


def _resolve_stage_contract_overrides(
    *,
    stage_contract_registry: tuple[PreprocessingStageMetadata, ...] | None,
    stage_metadata_registry: tuple[PreprocessingStageMetadata, ...] | None,
) -> tuple[PreprocessingStageMetadata, ...] | None:
    if stage_metadata_registry is not None:
        warn_deprecated(
            "preprocessing.pipeline.stage_metadata_registry",
            stacklevel=3,
        )
    if stage_contract_registry is not None and stage_metadata_registry is not None:
        raise DatasetBuildError(
            "dataset preprocessing registry arguments stage_contract_registry and "
            "stage_metadata_registry are aliases; only one may be passed"
        )
    if stage_metadata_registry is not None:
        return stage_metadata_registry
    return stage_contract_registry


def _normalize_report_rows(
    rows: Sequence[PreprocessingReportRow],
) -> tuple[PreprocessingReportRow, ...]:
    normalized: list[PreprocessingReportRow] = []
    for row in rows:
        normalized.append(validate_preprocessing_report_row(row))
    return tuple(normalized)


def _validate_quantitative_contract_before_execution(
    *,
    quantitative_state: QuantitativeContractState,
    stage_key: str,
    interpreted_contract: InterpretedPreprocessingStageContract,
) -> QuantitativeContractState:
    return interpreted_contract.quantitative_contract.validate_and_transition(
        quantitative_state,
        stage=stage_key,
        operation=interpreted_contract.operation,
        evidence=None,
    )


_INITIAL_PREPROCESSING_STATE_TABLES: frozenset[PreprocessingStateTableKey] = frozenset(
    {
        PreprocessingStateTableKey.DATASET_PHOSPHO,
        PreprocessingStateTableKey.DATASET_SITE_METADATA,
        PreprocessingStateTableKey.DATASET_SAMPLE_METADATA,
        PreprocessingStateTableKey.DATASET_TOTAL,
        PreprocessingStateTableKey.DATASET_IMPUTATION_OBSERVATION_MASK,
    }
)


def _validate_declared_table_dependencies(
    *,
    stage_key: str,
    consumed_input_tables: tuple[PreprocessingStateTableKey, ...],
    available_tables: set[PreprocessingStateTableKey],
) -> None:
    missing_tables = tuple(
        table for table in consumed_input_tables if table not in available_tables
    )
    if not missing_tables:
        return
    raise DatasetBuildError(
        "dataset preprocessing plan has invalid stage table dependencies: "
        f"stage={stage_key!r} consumes tables before they are available: "
        + ", ".join(table.value for table in missing_tables)
    )


def _initial_quantitative_meaning_evidence_mode(
    *,
    initial_quantitative_meaning: QuantitativeMeaning | None,
) -> QuantitativeMeaningEvidenceMode:
    if initial_quantitative_meaning is not None:
        return QuantitativeMeaningEvidenceMode.DECLARED_BY_CALLER
    return QuantitativeMeaningEvidenceMode.INFERRED_FROM_SCALE_CONTRACT


def _resolve_output_quantitative_meaning_evidence_mode(
    *,
    input_state: QuantitativeContractState,
    output_state: QuantitativeContractState,
    input_mode: QuantitativeMeaningEvidenceMode,
) -> QuantitativeMeaningEvidenceMode:
    if output_state.meaning is input_state.meaning:
        return input_mode
    return QuantitativeMeaningEvidenceMode.DERIVED_BY_PHOSPY_OPERATION


def _run_stage_pre_quantitative_contract_validation(
    *,
    stage: PreprocessingStage,
    state: PreprocessingState,
) -> None:
    stage.validate_before_quantitative_contract(state)


__all__ = ["PreprocessingPipeline", "_resolve_state_table"]
