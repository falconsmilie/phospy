from __future__ import annotations

from dataclasses import replace

from phospy.science.datasets.processing_state import (
    DatasetProcessingState,
    MissingDataDiagnosticsV1,
)


def imputed_processing_state(
    processing_state: DatasetProcessingState,
) -> DatasetProcessingState:
    """Return a valid imputed variant of a dataset processing state."""

    diagnostics = MissingDataDiagnosticsV1(
        missing_data_policy="impute_row_median",
        imputation_method_id="row_median",
        imputation_method_family="deterministic_row_statistic",
        input_missing_cell_count=1,
        output_missing_cell_count=0,
        imputed_cell_count=1,
        affected_row_count=1,
        affected_column_count=1,
        affected_row_ids=("row_a",),
        affected_column_ids=("sample_a",),
        imputed_row_ids=("row_a",),
        imputed_column_ids=("sample_a",),
        dropped_row_ids=(),
        method_parameters={"min_observed_values": 1},
        stage_order=("missing_data",),
        missingness_mask_hash="test-missingness-mask",
        imputation_mask_hash="test-imputation-mask",
        rows_not_imputable=(),
    )
    return replace(
        processing_state,
        missing_data=replace(
            processing_state.missing_data,
            policy="impute_row_median",
            min_observed_values=1,
            complete_matrix=True,
            imputed=True,
            diagnostics=diagnostics,
            has_missing_values=False,
            missing_value_count=0,
        ),
    )
