"""Dataset-facing adapter for the published Preprocessing boundary."""

from __future__ import annotations

import pandas as pd

from phospy.datasets.construction.contracts import PreprocessedDatasetBuildTables
from phospy.datasets.construction.protein_aware import (
    compose_protein_aware_preparation_result,
)
from phospy.datasets.construction.state_assembly import (
    DatasetProcessingStateBuilder,
)
from phospy.datasets.processing_state import DatasetProcessingState
from phospy.datasets.protein_aware import ProteinAwarePreparationResult
from phospy.preprocessing.contracts import (
    PreprocessingPlan,
    PreprocessingStageExecution,
)
from phospy.preprocessing.dataset_construction import (
    BatchCorrectionAdequacyValidatorProtocol,
    BatchDesignMetadataValidatorProtocol,
    CorrectedPreprocessingOutput,
    CorrectedPreprocessingOutputIntegratorProtocol,
    DatasetConstructionPreprocessingRunner,
    DatasetConstructionProteinAwareRunner,
    PreprocessingPipelineProtocol,
    PreprocessingProvenanceAdapterProtocol,
    ProteinAwarePreparationStageProtocol,
    SpsRuvStyleBatchCorrectionRunner,
)
from phospy.quantitative.models import (
    IntensityScaleKind,
    IntensityScaleState,
    QuantitativeMeaning,
)
from phospy.quantitative.quantitative_contracts import (
    QuantitativeContractState,
)


class DatasetPreprocessor:
    """Translate builder input into internal preprocessing state and run it."""

    def __init__(
        self,
        *,
        runner: DatasetConstructionPreprocessingRunner | None = None,
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
        self._runner = runner or DatasetConstructionPreprocessingRunner(
            pipeline=pipeline,
            batch_correction_runner=batch_correction_runner,
            provenance_adapter=provenance_adapter,
            correction_integrator=correction_integrator,
            batch_correction_metadata_validator=batch_correction_metadata_validator,
            batch_correction_adequacy_validator=batch_correction_adequacy_validator,
        )

    def validate_quantitative_contracts(
        self,
        *,
        plan: PreprocessingPlan,
        initial_quantitative_scale_kind: IntensityScaleKind | None = None,
        initial_quantitative_meaning: QuantitativeMeaning | None = None,
    ) -> QuantitativeContractState:
        return self._runner.validate_quantitative_contracts(
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
    ) -> PreprocessedDatasetBuildTables:
        output = self._runner.run(
            phospho=phospho,
            site_metadata=site_metadata,
            sample_metadata=sample_metadata,
            total=total,
            plan=plan,
            corrected_preprocessing_output=corrected_preprocessing_output,
            initial_quantitative_scale_kind=initial_quantitative_scale_kind,
            initial_quantitative_meaning=initial_quantitative_meaning,
        )
        return PreprocessedDatasetBuildTables(
            phospho=output.phospho,
            site_metadata=output.site_metadata,
            sample_metadata=output.sample_metadata,
            total=output.total,
            comparisons=output.comparisons,
            imputation_observation_mask=output.imputation_observation_mask,
            comparison_group_stats=output.comparison_group_stats,
            comparison_pair_stats=output.comparison_pair_stats,
            preprocessing_row_counts=output.preprocessing_row_counts,
            preprocessing_operations=output.preprocessing_operations,
            row_audit=output.row_audit,
            preprocessing_trace=output.preprocessing_trace,
            duplicate_site_resolution=output.duplicate_site_resolution,
            metadata_conflicts=output.metadata_conflicts,
            batch_correction_metadata=output.batch_correction_metadata,
            batch_correction_report=output.batch_correction_report,
            intensity_transformation_event=output.intensity_transformation_event,
        )


class DatasetProteinAwarePreparationRunner:
    """Builder-facing adapter for protein-aware preparation after scale setup."""

    def __init__(
        self,
        *,
        runner: DatasetConstructionProteinAwareRunner | None = None,
        stage: ProteinAwarePreparationStageProtocol | None = None,
    ) -> None:
        self._runner = runner or DatasetConstructionProteinAwareRunner(stage=stage)

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        total: pd.DataFrame | None,
        intensity_scale_state: IntensityScaleState,
        plan: PreprocessingPlan,
    ) -> ProteinAwarePreparationResult | None:
        stage_output = self._runner.run(
            phospho=phospho,
            site_metadata=site_metadata,
            total=total,
            intensity_scale_state=intensity_scale_state,
            plan=plan,
        )
        if stage_output is None:
            return None
        return compose_protein_aware_preparation_result(
            stage_output,
            provenance=None,
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
