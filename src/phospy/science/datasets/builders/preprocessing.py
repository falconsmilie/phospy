"""Builder-facing adapter for the internal preprocessing subsystem."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.configs import (
    DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH,
    DATASET_BATCH_CORRECTION_METHOD_NONE,
)
from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.science.datasets.preprocessing.batch_correction_metadata import (
    BatchCorrectionMetadataResolver,
)
from phospy.science.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingState,
)
from phospy.science.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.science.datasets.preprocessing.provenance_adapter import (
    PreprocessingProvenanceAdapter,
)
from phospy.science.datasets.preprocessing.report_rows import (
    compose_stage_owned_report_tables,
)
from phospy.science.datasets.preprocessing.state_builder import (
    DatasetProcessingStateBuilder,
)
from phospy.science.datasets.processing_state import DatasetProcessingState
from phospy.science.transformations.models import (
    IntensityScaleState,
    QuantitativeMeaning,
)
from phospy.validation.datasets.batch_correction import (
    BatchCorrectionAdequacyValidator,
)


class DatasetPreprocessor:
    """Translate builder input into internal preprocessing state and run it."""

    def __init__(
        self,
        *,
        pipeline: PreprocessingPipeline | None = None,
        provenance_adapter: PreprocessingProvenanceAdapter | None = None,
        batch_correction_metadata_resolver: BatchCorrectionMetadataResolver
        | None = None,
        batch_correction_adequacy_validator: BatchCorrectionAdequacyValidator
        | None = None,
    ) -> None:
        self._pipeline = pipeline or PreprocessingPipeline()
        self._provenance_adapter = (
            provenance_adapter or PreprocessingProvenanceAdapter()
        )
        self._batch_correction_metadata_resolver = (
            batch_correction_metadata_resolver or BatchCorrectionMetadataResolver()
        )
        self._batch_correction_adequacy_validator = (
            batch_correction_adequacy_validator or BatchCorrectionAdequacyValidator()
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        site_metadata: pd.DataFrame,
        sample_metadata: pd.DataFrame | None = None,
        total: pd.DataFrame | None = None,
        plan: PreprocessingPlan,
    ) -> PreprocessedDatasetBuildTables:
        input_row_count = int(len(phospho.index))
        batch_correction_metadata = None
        batch_correction_method = str(plan.batch_correction_method).strip()
        if batch_correction_method != DATASET_BATCH_CORRECTION_METHOD_NONE:
            if (
                batch_correction_method
                != DATASET_BATCH_CORRECTION_METHOD_LINEAR_RESIDUALIZE_BATCH
            ):
                raise PhosPyInputError(
                    "dataset preprocessing plan contains unsupported "
                    f"batch_correction_method={batch_correction_method!r}"
                )
            batch_correction_metadata = self._batch_correction_metadata_resolver.run(
                phospho=phospho,
                sample_metadata=sample_metadata,
                batch_column=plan.batch_correction_batch_column,
                condition_column=plan.batch_correction_condition_column,
            )
            self._batch_correction_adequacy_validator.run(
                batch_by_sample=batch_correction_metadata.batch_by_sample,
                condition_by_sample=batch_correction_metadata.condition_by_sample,
                sample_order=batch_correction_metadata.sample_order,
                preserve_condition_effects=(
                    plan.batch_correction_preserve_condition_effects
                ),
            )
        preprocessed_state, trace = self._pipeline.run_with_trace(
            PreprocessingState(
                phospho=phospho,
                site_metadata=site_metadata,
                sample_metadata=sample_metadata,
                total=total,
                plan=plan,
            )
        )
        row_counts, operations = self._provenance_adapter.build_tables(
            plan=plan,
            input_row_count=input_row_count,
            output_row_count=int(len(preprocessed_state.phospho.index)),
            trace=trace,
        )
        report_tables = compose_stage_owned_report_tables(
            preprocessed_state.report_rows
        )
        return PreprocessedDatasetBuildTables(
            phospho=preprocessed_state.phospho,
            site_metadata=preprocessed_state.site_metadata,
            sample_metadata=preprocessed_state.sample_metadata,
            total=preprocessed_state.total,
            comparisons=preprocessed_state.comparisons,
            comparison_group_stats=report_tables.comparison_group_stats,
            comparison_pair_stats=report_tables.comparison_pair_stats,
            preprocessing_row_counts=row_counts,
            preprocessing_operations=operations,
            row_audit=report_tables.row_audit,
            preprocessing_trace=trace,
            duplicate_site_resolution=report_tables.duplicate_site_resolution,
            metadata_conflicts=report_tables.metadata_conflicts,
            batch_correction_metadata=batch_correction_metadata,
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
