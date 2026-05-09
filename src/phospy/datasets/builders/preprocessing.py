"""Builder-facing adapter for the internal preprocessing subsystem."""

from __future__ import annotations

import pandas as pd

from phospy.datasets.builders.contracts import PreprocessedDatasetBuildTables
from phospy.datasets.preprocessing.models import (
    PreprocessingPlan,
    PreprocessingStageExecution,
    PreprocessingState,
)
from phospy.datasets.preprocessing.pipeline import PreprocessingPipeline
from phospy.datasets.preprocessing.provenance_adapter import (
    PreprocessingProvenanceAdapter,
)
from phospy.datasets.preprocessing.report_rows import (
    compose_stage_owned_report_tables,
)
from phospy.datasets.preprocessing.state_builder import DatasetProcessingStateBuilder
from phospy.datasets.processing_state import DatasetProcessingState
from phospy.transformations.models import IntensityScaleState, QuantitativeMeaning


class DatasetPreprocessor:
    """Translate builder input into internal preprocessing state and run it."""

    def __init__(
        self,
        *,
        pipeline: PreprocessingPipeline | None = None,
        provenance_adapter: PreprocessingProvenanceAdapter | None = None,
    ) -> None:
        self._pipeline = pipeline or PreprocessingPipeline()
        self._provenance_adapter = (
            provenance_adapter or PreprocessingProvenanceAdapter()
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
