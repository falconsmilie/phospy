"""Builder-facing adapter for the internal preprocessing subsystem."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs import DATASET_BATCH_CORRECTION_METHOD_NONE
from phospy.contracts.configs.preprocessing.total_protein import (
    DatasetProteinAwarePreparationConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.science.datasets.preprocessing.correction_output import (
    CorrectedPreprocessingOutput,
    CorrectedPreprocessingOutputIntegrator,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.protein_aware_preparation import (
    ProteinAwarePreparationResult,
    ProteinAwarePreparationStage,
)
from phospy.science.datasets.preprocessing.provenance_adapter import (
    PreprocessingProvenanceAdapter,
)
from phospy.science.datasets.preprocessing.report_rows import (
    compose_stage_owned_report_tables,
)
from phospy.science.datasets.preprocessing.stages import (
    SpsRuvStyleBatchCorrectionRunner,
)
from phospy.science.datasets.preprocessing.state_builder import (
    DatasetProcessingStateBuilder,
)
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.transformations.models import (
    IntensityScaleState,
    IntensityTransformationEvent,
    QuantitativeMeaning,
)
from phospy.validation.datasets.preprocessing import (
    reject_external_corrected_output_after_downstream_preprocessing,
)


class DatasetPreprocessor:
    """Translate builder input into internal preprocessing state and run it."""

    def __init__(
        self,
        *,
        pipeline: PreprocessingPipeline | None = None,
        batch_correction_runner: SpsRuvStyleBatchCorrectionRunner | None = None,
        provenance_adapter: PreprocessingProvenanceAdapter | None = None,
        correction_integrator: CorrectedPreprocessingOutputIntegrator | None = None,
    ) -> None:
        self._pipeline = pipeline or PreprocessingPipeline(
            batch_correction_runner=batch_correction_runner
        )
        self._provenance_adapter = (
            provenance_adapter or PreprocessingProvenanceAdapter()
        )
        self._correction_integrator = (
            correction_integrator or CorrectedPreprocessingOutputIntegrator()
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
    ) -> PreprocessedDatasetBuildTables:
        if (
            corrected_preprocessing_output is not None
            and plan.batch_correction_method != DATASET_BATCH_CORRECTION_METHOD_NONE
        ):
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
        preprocessed_state, trace = self._pipeline.run_with_trace(
            PreprocessingState(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                plan=plan,
            )
        )
        if corrected_preprocessing_output is not None:
            preprocessed_state, correction_trace = self._correction_integrator.run(
                state=preprocessed_state,
                correction_output=corrected_preprocessing_output,
            )
            trace = (*trace, correction_trace)
        row_counts, operations = self._provenance_adapter.build_tables(
            plan=plan,
            input_row_count=input_row_count,
            output_row_count=int(len(preprocessed_state.phospho.index)),
            trace=trace,
        )
        report_tables = compose_stage_owned_report_tables(
            preprocessed_state.report_rows
        )
        intensity_transformation_event = _resolve_intensity_transformation_event(trace)
        return PreprocessedDatasetBuildTables(
            phospho=preprocessed_state.phospho,
            site_metadata=preprocessed_state.site_metadata,
            sample_metadata=preprocessed_state.sample_metadata,
            total=preprocessed_state.total,
            comparisons=preprocessed_state.comparisons,
            imputation_observation_mask=(
                preprocessed_state.imputation_observation_mask
            ),
            comparison_group_stats=report_tables.comparison_group_stats,
            comparison_pair_stats=report_tables.comparison_pair_stats,
            preprocessing_row_counts=row_counts,
            preprocessing_operations=operations,
            row_audit=report_tables.row_audit,
            preprocessing_trace=trace,
            duplicate_site_resolution=report_tables.duplicate_site_resolution,
            metadata_conflicts=report_tables.metadata_conflicts,
            batch_correction_metadata=preprocessed_state.batch_correction_metadata,
            batch_correction_report=preprocessed_state.batch_correction_report,
            intensity_transformation_event=intensity_transformation_event,
        )


class DatasetProteinAwarePreparationRunner:
    """Builder-facing adapter for protein-aware preparation after scale setup."""

    def __init__(
        self,
        *,
        stage: ProteinAwarePreparationStage | None = None,
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
    ) -> ProteinAwarePreparationResult | None:
        return self._stage.run(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            transformation_state=intensity_scale_state,
            config=DatasetProteinAwarePreparationConfig(
                policy=plan.protein_aware_preparation_policy,
                protein_mapping_policy=(plan.protein_aware_preparation_mapping_policy),
            ),
        )


_PROCESSING_STATE_BUILDER = DatasetProcessingStateBuilder()


def build_dataset_processing_state(
    *,
    plan: PreprocessingPlan,
    intensity_scale_state: IntensityScaleState,
    explicit_quantitative_meaning: QuantitativeMeaning | None = None,
    preprocessing_trace: tuple[PreprocessingStageExecution, ...] | None = None,
    final_phospho: pd.DataFrame | None = None,
    final_site_metadata: pd.DataFrame | None = None,
    final_sample_metadata: pd.DataFrame | None = None,
) -> DatasetProcessingState:
    """Build compact dataset processing state from the resolved preprocessing plan."""

    return _PROCESSING_STATE_BUILDER.build(
        plan=plan,
        intensity_scale_state=intensity_scale_state,
        explicit_quantitative_meaning=explicit_quantitative_meaning,
        preprocessing_trace=preprocessing_trace,
        final_phospho=final_phospho,
        final_site_metadata=final_site_metadata,
        final_sample_metadata=final_sample_metadata,
    )


def _resolve_intensity_transformation_event(
    trace: tuple[PreprocessingStageExecution, ...],
) -> IntensityTransformationEvent | None:
    event: IntensityTransformationEvent | None = None
    for record in trace:
        if record.intensity_transformation_event is None:
            continue
        event = record.intensity_transformation_event
    return event
