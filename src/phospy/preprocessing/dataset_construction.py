"""Published Preprocessing composition boundary for Dataset construction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.preparation_records.batch_correction import BatchCorrectionReport
from phospy.preprocessing.batch_correction.adapters.preprocessing import (
    SpsRuvStyleBatchCorrectionWorkflowRunner,
)
from phospy.preprocessing.batch_correction.contracts.output import (
    CorrectedPreprocessingOutput,
)
from phospy.preprocessing.batch_correction.contracts.stage import (
    BatchCorrectionAdequacyValidatorProtocol,
    BatchDesignMetadataValidatorProtocol,
    SpsRuvStyleBatchCorrectionRunner,
)
from phospy.preprocessing.batch_correction.metadata import (
    ResolvedBatchCorrectionMetadata,
    levels_in_sample_order,
)
from phospy.preprocessing.batch_correction.validation.design import (
    BatchCorrectionAdequacyValidator,
    BatchDesignMetadataValidator,
)
from phospy.preprocessing.batch_correction.validation.output import (
    validate_corrected_preprocessing_output,
)
from phospy.preprocessing.contracts import (
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingState,
    reject_external_corrected_output_after_downstream_preprocessing,
)
from phospy.preprocessing.contracts.configuration.total_protein import (
    DatasetProteinAwarePreparationConfig,
)
from phospy.preprocessing.contracts.stage import PreprocessingStageContract
from phospy.preprocessing.execution.corrected_output import (
    CorrectedPreprocessingOutputIntegrator,
)
from phospy.preprocessing.execution.pipeline import PreprocessingPipeline
from phospy.preprocessing.planning.registry import get_preprocessing_stage_metadata
from phospy.preprocessing.provenance.composition import PreprocessingProvenanceAdapter
from phospy.preprocessing.results.diagnostics import ProcessingTraceDiagnostics
from phospy.preprocessing.results.protein_aware import (
    ProteinAwarePreparationStageOutput,
)
from phospy.preprocessing.results.quantitative_evidence import (
    QuantitativeOperationEvidenceContext,
)
from phospy.preprocessing.results.report_rows import compose_stage_owned_report_tables
from phospy.preprocessing.stages.protein_aware import ProteinAwarePreparationStage
from phospy.preprocessing.validation.quantitative_evidence import (
    validate_quantitative_operation_evidence,
)
from phospy.quantitative.models import (
    IntensityScaleKind,
    IntensityScaleState,
    IntensityTransformationEvent,
    QuantitativeMeaning,
)
from phospy.quantitative.quantitative_contracts import QuantitativeContractState


class PreprocessingPipelineProtocol(Protocol):
    """Pipeline behaviour accepted by the published composition boundary."""

    def validate_quantitative_contracts(
        self,
        *,
        plan: PreprocessingPlan,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> QuantitativeContractState: ...

    def run_with_trace(
        self,
        state: PreprocessingState,
        *,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> tuple[PreprocessingState, tuple[PreprocessingStageExecution, ...]]: ...


class PreprocessingProvenanceAdapterProtocol(Protocol):
    """Passive provenance-table adapter used by Dataset construction."""

    def build_tables(
        self,
        *,
        plan: PreprocessingPlan,
        input_row_count: int,
        output_row_count: int,
        trace: tuple[PreprocessingStageExecution, ...],
    ) -> tuple[pd.DataFrame, pd.DataFrame]: ...


class CorrectedPreprocessingOutputIntegratorProtocol(Protocol):
    """Trusted corrected-output integration used by the composition boundary."""

    def run(
        self,
        *,
        state: PreprocessingState,
        correction_output: CorrectedPreprocessingOutput,
    ) -> tuple[PreprocessingState, PreprocessingStageExecution]: ...


class ProteinAwarePreparationStageProtocol(Protocol):
    """Protein-aware stage behaviour hidden behind the published boundary."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        total: pd.DataFrame | None,
        transformation_state: IntensityScaleState,
        config: DatasetProteinAwarePreparationConfig,
    ) -> ProteinAwarePreparationStageOutput | None: ...


@dataclass(frozen=True, slots=True)
class DatasetConstructionPreprocessingOutput:
    """Passive preprocessing output consumed by Dataset durable assembly."""

    phospho: pd.DataFrame
    site_metadata: pd.DataFrame
    sample_metadata: pd.DataFrame | None
    total: pd.DataFrame | None
    comparisons: pd.DataFrame | None
    imputation_observation_mask: pd.DataFrame | None
    comparison_group_stats: pd.DataFrame
    comparison_pair_stats: pd.DataFrame
    preprocessing_row_counts: pd.DataFrame
    preprocessing_operations: pd.DataFrame
    row_audit: pd.DataFrame
    preprocessing_trace: tuple[PreprocessingStageExecution, ...]
    duplicate_site_resolution: pd.DataFrame
    metadata_conflicts: pd.DataFrame
    batch_correction_metadata: ResolvedBatchCorrectionMetadata | None
    batch_correction_report: BatchCorrectionReport | None
    intensity_transformation_event: IntensityTransformationEvent | None


class DatasetConstructionPreprocessingRunner:
    """Execute preprocessing and publish only passive construction output."""

    def __init__(
        self,
        *,
        pipeline: PreprocessingPipelineProtocol | None = None,
        batch_correction_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
        provenance_adapter: PreprocessingProvenanceAdapterProtocol | None = None,
        correction_integrator: (
            CorrectedPreprocessingOutputIntegratorProtocol | None
        ) = None,
        batch_correction_metadata_validator: (
            BatchDesignMetadataValidatorProtocol | None
        ) = None,
        batch_correction_adequacy_validator: (
            BatchCorrectionAdequacyValidatorProtocol | None
        ) = None,
    ) -> None:
        self._pipeline = pipeline or PreprocessingPipeline(
            batch_correction_runner=batch_correction_runner,
            batch_correction_metadata_validator=batch_correction_metadata_validator,
            batch_correction_adequacy_validator=batch_correction_adequacy_validator,
        )
        self._provenance_adapter = (
            provenance_adapter or PreprocessingProvenanceAdapter()
        )
        self._correction_integrator = (
            correction_integrator or CorrectedPreprocessingOutputIntegrator()
        )

    def validate_quantitative_contracts(
        self,
        *,
        plan: PreprocessingPlan,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> QuantitativeContractState:
        return self._pipeline.validate_quantitative_contracts(
            plan=plan,
            initial_quantitative_scale_kind=initial_quantitative_scale_kind,
            initial_quantitative_meaning=initial_quantitative_meaning,
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        plan: PreprocessingPlan,
        corrected_preprocessing_output: CorrectedPreprocessingOutput | None = None,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> DatasetConstructionPreprocessingOutput:
        if corrected_preprocessing_output is not None and plan.batch_correction_enabled:
            raise PhosPyInputError(
                "dataset preprocessing received corrected_preprocessing_output "
                "while preprocessing_config.batch_correction also requests "
                "execution; correction must be applied exactly once"
            )
        if corrected_preprocessing_output is not None:
            reject_external_corrected_output_after_downstream_preprocessing(
                plan.stage_order
            )
        input_row_count = int(len(phospho.index))
        state, trace = self._pipeline.run_with_trace(
            PreprocessingState(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                plan=plan,
            ),
            initial_quantitative_scale_kind=initial_quantitative_scale_kind,
            initial_quantitative_meaning=initial_quantitative_meaning,
        )
        if corrected_preprocessing_output is not None:
            state, correction_trace = self._correction_integrator.run(
                state=state,
                correction_output=corrected_preprocessing_output,
            )
            trace = (*trace, correction_trace)
        row_counts, operations = self._provenance_adapter.build_tables(
            plan=plan,
            input_row_count=input_row_count,
            output_row_count=int(len(state.phospho.index)),
            trace=trace,
        )
        report_tables = compose_stage_owned_report_tables(state.report_rows)
        return DatasetConstructionPreprocessingOutput(
            phospho=state.phospho,
            site_metadata=state.site_metadata,
            sample_metadata=state.sample_metadata,
            total=state.total,
            comparisons=state.comparisons,
            imputation_observation_mask=state.imputation_observation_mask,
            comparison_group_stats=report_tables.comparison_group_stats,
            comparison_pair_stats=report_tables.comparison_pair_stats,
            preprocessing_row_counts=row_counts,
            preprocessing_operations=operations,
            row_audit=report_tables.row_audit,
            preprocessing_trace=trace,
            duplicate_site_resolution=report_tables.duplicate_site_resolution,
            metadata_conflicts=report_tables.metadata_conflicts,
            batch_correction_metadata=state.batch_correction_metadata,
            batch_correction_report=state.batch_correction_report,
            intensity_transformation_event=_resolve_intensity_transformation_event(
                trace
            ),
        )


class DatasetConstructionProteinAwareRunner:
    """Execute the Preprocessing-owned protein-aware stage."""

    def __init__(
        self,
        *,
        stage: ProteinAwarePreparationStageProtocol | None = None,
    ) -> None:
        self._stage = stage or ProteinAwarePreparationStage()

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        total: pd.DataFrame | None,
        intensity_scale_state: IntensityScaleState,
        plan: PreprocessingPlan,
    ) -> ProteinAwarePreparationStageOutput | None:
        return self._stage.run(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            transformation_state=intensity_scale_state,
            config=DatasetProteinAwarePreparationConfig(
                policy=plan.protein_aware_preparation_policy,
                protein_mapping_policy=plan.protein_aware_preparation_mapping_policy,
            ),
        )


def create_default_dataset_preprocessing_runner(
    *,
    batch_correction_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
) -> DatasetConstructionPreprocessingRunner:
    """Create the supported public builder's default Preprocessing composition."""

    return DatasetConstructionPreprocessingRunner(
        batch_correction_runner=(
            batch_correction_runner or SpsRuvStyleBatchCorrectionWorkflowRunner()
        ),
        batch_correction_metadata_validator=BatchDesignMetadataValidator(),
        batch_correction_adequacy_validator=BatchCorrectionAdequacyValidator(),
    )


def dataset_construction_stage_metadata(
    stage_key: str,
) -> PreprocessingStageContract:
    """Return registered passive stage metadata for durable Dataset assembly."""

    return get_preprocessing_stage_metadata(stage_key)


def validate_dataset_quantitative_operation_evidence(
    context: QuantitativeOperationEvidenceContext,
) -> None:
    """Validate passive quantitative evidence for durable Dataset assembly."""

    validate_quantitative_operation_evidence(context)


def parse_processing_trace_diagnostics(
    trace: tuple[PreprocessingStageExecution, ...] | None,
) -> ProcessingTraceDiagnostics:
    """Parse passive execution diagnostics for durable Dataset assembly."""

    return ProcessingTraceDiagnostics.from_trace(trace)


def _resolve_intensity_transformation_event(
    trace: tuple[PreprocessingStageExecution, ...],
) -> IntensityTransformationEvent | None:
    event: IntensityTransformationEvent | None = None
    for record in trace:
        if record.intensity_transformation_event is not None:
            event = record.intensity_transformation_event
    return event


__all__ = [
    "CorrectedPreprocessingOutputIntegratorProtocol",
    "CorrectedPreprocessingOutput",
    "DatasetConstructionPreprocessingOutput",
    "DatasetConstructionPreprocessingRunner",
    "DatasetConstructionProteinAwareRunner",
    "PreprocessingPipelineProtocol",
    "PreprocessingProvenanceAdapterProtocol",
    "ProcessingTraceDiagnostics",
    "ProteinAwarePreparationStageProtocol",
    "ProteinAwarePreparationStageOutput",
    "QuantitativeOperationEvidenceContext",
    "ResolvedBatchCorrectionMetadata",
    "create_default_dataset_preprocessing_runner",
    "dataset_construction_stage_metadata",
    "levels_in_sample_order",
    "parse_processing_trace_diagnostics",
    "validate_corrected_preprocessing_output",
    "validate_dataset_quantitative_operation_evidence",
]
