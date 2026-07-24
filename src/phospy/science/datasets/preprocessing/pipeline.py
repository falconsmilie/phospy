"""Internal dataset preprocessing pipeline orchestration."""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import replace

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
    PreprocessingReportRow,
    PreprocessingStage,
    PreprocessingStageExecution,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.report_rows import (
    validate_preprocessing_report_row,
)
from phospy.science.datasets.preprocessing.stage_contract import (
    PreprocessingStageFactoryContext,
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

_LEGACY_STAGE_METADATA_REGISTRY_DEPRECATION_MESSAGE = (
    "PreprocessingPipeline(stage_metadata_registry=...) is deprecated because "
    "stage_metadata_registry is a legacy alias for stage_contract_registry; use "
    "stage_contract_registry instead. The legacy alias is planned for removal in "
    "PhosPy 1.8.0."
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
        self._stages_by_key = {stage.stage_key: stage for stage in stages}
        if len(self._stages_by_key) != len(stages):
            raise DatasetBuildError(
                "dataset preprocessing stage registry contains duplicate stage keys"
            )
        self._stage_contract_by_key = {
            metadata.stage_key: metadata for metadata in resolved_metadata_registry
        }

    def run(self, state: PreprocessingState) -> PreprocessingState:
        final_state, _ = self.run_with_trace(state)
        return final_state

    def run_with_trace(
        self,
        state: PreprocessingState,
    ) -> tuple[PreprocessingState, tuple[PreprocessingStageExecution, ...]]:
        current = state
        trace: list[PreprocessingStageExecution] = []
        report_rows: list[PreprocessingReportRow] = list(current.report_rows)
        event_validator = _TransformationEventSequenceValidator()
        diagnostics_defaults_resolver = _StageDiagnosticsDefaultsResolver()
        diagnostics_normalizer = _StageDiagnosticsNormalizer()
        fingerprint_service = _StageFingerprintService()
        trace_builder = _StageTraceBuilder()

        for stage_key in current.plan.stage_order:
            stage = self._resolve_stage(stage_key)
            contract = self._resolve_stage_contract(stage_key)
            previous = current
            interpreted_contract = contract.interpret(previous.plan)
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
            trace.append(
                trace_builder.run(
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
            )
        if report_rows:
            current = replace(current, report_rows=tuple(report_rows))
        return current, tuple(trace)

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


def _resolve_stage_contract_overrides(
    *,
    stage_contract_registry: tuple[PreprocessingStageMetadata, ...] | None,
    stage_metadata_registry: tuple[PreprocessingStageMetadata, ...] | None,
) -> tuple[PreprocessingStageMetadata, ...] | None:
    if stage_contract_registry is not None and stage_metadata_registry is not None:
        raise DatasetBuildError(
            "dataset preprocessing registry arguments stage_contract_registry and "
            "stage_metadata_registry are aliases; only one may be passed"
        )
    if stage_metadata_registry is not None:
        warnings.warn(
            _LEGACY_STAGE_METADATA_REGISTRY_DEPRECATION_MESSAGE,
            DeprecationWarning,
            stacklevel=3,
        )
        return stage_metadata_registry
    return stage_contract_registry


def _normalize_report_rows(
    rows: Sequence[PreprocessingReportRow],
) -> tuple[PreprocessingReportRow, ...]:
    normalized: list[PreprocessingReportRow] = []
    for row in rows:
        normalized.append(validate_preprocessing_report_row(row))
    return tuple(normalized)


__all__ = ["PreprocessingPipeline", "_resolve_state_table"]
