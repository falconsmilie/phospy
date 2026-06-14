"""Row-median missing-data policy implementation."""

from __future__ import annotations

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.models import PreprocessingState

from .models import RowImputationRecord, RowMedianPolicyOutcome


def run_row_median_policy(state: PreprocessingState) -> RowMedianPolicyOutcome:
    """Apply row-median policy numerical transformation."""

    min_observed_values = state.plan.missing_data_min_observed_values
    if not isinstance(min_observed_values, int):
        raise PhosPyInputError(
            "dataset build request "
            "preprocessing_config.missing_data.min_observed_values must be an "
            "int when missing_data.policy='impute_row_median'"
        )
    if min_observed_values > state.phospho.shape[1]:
        raise PhosPyInputError(
            "dataset build request "
            "preprocessing_config.missing_data.min_observed_values "
            f"({min_observed_values}) cannot exceed the number of phospho "
            f"samples ({state.phospho.shape[1]})"
        )

    observed_counts = state.phospho.notna().sum(axis=1)
    retained_mask = observed_counts >= min_observed_values
    dropped_observed_counts = observed_counts.loc[~retained_mask]
    filtered_phospho = state.phospho.loc[retained_mask]
    filtered_site_metadata = state.site_metadata.loc[filtered_phospho.index]

    row_medians = filtered_phospho.median(axis=1, skipna=True)
    imputed = filtered_phospho.T.fillna(row_medians).T
    imputed_mask = filtered_phospho.isna() & imputed.notna()

    imputed_rows = filtered_phospho.index[imputed_mask.any(axis=1)]
    imputed_cell_count = int(imputed_mask.to_numpy().sum())
    imputed_row_ids = tuple(str(row_id) for row_id in imputed_rows.tolist())
    imputed_column_ids = tuple(
        str(column_name)
        for column_name in filtered_phospho.columns[imputed_mask.any(axis=0)]
    )
    rows_not_imputable = tuple(
        str(row_id)
        for row_id in imputed.index[
            imputed.isna().any(axis=1) & filtered_phospho.isna().any(axis=1)
        ].tolist()
    )
    dropped_row_ids = tuple(
        str(row_id) for row_id in dropped_observed_counts.index.tolist()
    )
    output_missing_cell_count = int(imputed.isna().to_numpy().sum())
    row_medians_used = {
        str(row_id): float(row_medians.loc[row_id]) for row_id in imputed_rows
    }
    dropped_row_observed_values = tuple(
        (str(row_id), int(observed_value_count))
        for row_id, observed_value_count in dropped_observed_counts.items()
    )
    imputed_rows_audit = tuple(
        RowImputationRecord(
            row_id=str(row_id),
            imputed_columns=tuple(
                str(column_name)
                for column_name in filtered_phospho.columns[
                    imputed_mask.loc[row_id]
                ].tolist()
            ),
            imputed_cell_count=int(imputed_mask.loc[row_id].sum()),
        )
        for row_id in imputed_rows
    )
    return RowMedianPolicyOutcome(
        phospho=imputed,
        site_metadata=filtered_site_metadata,
        imputed_mask=imputed_mask,
        min_observed_values=int(min_observed_values),
        dropped_row_ids=dropped_row_ids,
        dropped_row_observed_values=dropped_row_observed_values,
        imputed_cell_count=imputed_cell_count,
        imputed_row_ids=imputed_row_ids,
        imputed_column_ids=imputed_column_ids,
        output_missing_cell_count=output_missing_cell_count,
        rows_not_imputable=rows_not_imputable,
        row_medians_used=row_medians_used,
        imputed_rows=imputed_rows_audit,
    )
