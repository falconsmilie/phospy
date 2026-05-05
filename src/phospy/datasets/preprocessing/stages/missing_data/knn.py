"""KNN missing-data policy implementation."""

from __future__ import annotations

import pandas as pd
from sklearn.impute import KNNImputer

from phospy.datasets.preprocessing.models import PreprocessingState
from phospy.errors.input import PhosPyInputError

from .diagnostics import label_preview
from .models import KnnPolicyOutcome, RowImputationRecord


def run_knn_policy(state: PreprocessingState) -> KnnPolicyOutcome:
    """Apply KNN policy numerical transformation."""

    k = state.plan.missing_data_k
    distance = state.plan.missing_data_distance
    max_missing_fraction_per_row = state.plan.missing_data_max_missing_fraction_per_row
    if k is None or distance is None or max_missing_fraction_per_row is None:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.missing_data.policy="
            "'impute_knn' requires k, distance, and "
            "max_missing_fraction_per_row"
        )
    k_value = int(k)
    distance_value = str(distance).strip()
    max_missing_fraction_value = float(max_missing_fraction_per_row)

    missing_fraction = state.phospho.isna().mean(axis=1)
    retained_mask = missing_fraction <= max_missing_fraction_value
    dropped_missing_fraction = missing_fraction.loc[~retained_mask]
    filtered_phospho = state.phospho.loc[retained_mask].copy(deep=True)
    filtered_site_metadata = state.site_metadata.loc[filtered_phospho.index]

    if filtered_phospho.empty:
        imputed = filtered_phospho.copy(deep=True)
    else:
        all_missing_columns = filtered_phospho.columns[
            filtered_phospho.notna().sum(axis=0).to_numpy(dtype=int, copy=False) == 0
        ]
        if len(all_missing_columns) > 0:
            raise PhosPyInputError(
                "dataset preprocessing stage 'missing_data' cannot apply "
                "missing_data.policy='impute_knn' because one or more columns "
                "have no observed values after row filtering. "
                f"affected column labels (preview): {label_preview(all_missing_columns.tolist())}. "
                "adjust missing_data.max_missing_fraction_per_row or input data."
            )
        imputer = KNNImputer(
            n_neighbors=k_value,
            metric="nan_euclidean",
        )
        imputed_values = imputer.fit_transform(filtered_phospho)
        if imputed_values.shape[1] != filtered_phospho.shape[1]:
            raise PhosPyInputError(
                "dataset preprocessing stage 'missing_data' cannot apply "
                "missing_data.policy='impute_knn' because the imputer could not "
                "retain all matrix columns during imputation. "
                "ensure every retained column has at least one observed value."
            )
        imputed = pd.DataFrame(
            imputed_values,
            index=filtered_phospho.index.copy(),
            columns=filtered_phospho.columns.copy(),
        )

    if filtered_phospho.empty:
        imputed_mask = filtered_phospho.isna() & filtered_phospho.notna()
    else:
        imputed_mask = filtered_phospho.isna() & imputed.notna()

    imputed_cell_count = int(imputed_mask.to_numpy().sum())
    imputed_row_ids = (
        tuple(
            str(row_id)
            for row_id in imputed.index[
                imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not imputed.empty
        else ()
    )
    imputed_column_ids = (
        tuple(
            str(column_name)
            for column_name in imputed.columns[
                imputed_mask.any(axis=0).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not imputed.empty
        else ()
    )
    unresolved_mask = imputed.isna() & filtered_phospho.isna()
    unresolved_row_ids = (
        tuple(
            str(row_id)
            for row_id in imputed.index[
                unresolved_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not imputed.empty
        else ()
    )
    unresolved_column_ids = (
        tuple(
            str(column_name)
            for column_name in imputed.columns[
                unresolved_mask.any(axis=0).to_numpy(dtype=bool, copy=False)
            ].tolist()
        )
        if not imputed.empty
        else ()
    )
    output_missing_cell_count = int(imputed.isna().to_numpy().sum())
    if output_missing_cell_count > 0:
        raise PhosPyInputError(
            "dataset preprocessing stage 'missing_data' could not complete "
            "missing_data.policy='impute_knn' because missing values remain after "
            "imputation. "
            f"remaining rows (preview): {label_preview(list(unresolved_row_ids))}. "
            f"remaining columns (preview): {label_preview(list(unresolved_column_ids))}. "
            "adjust missing_data.max_missing_fraction_per_row, k, or input data."
        )

    dropped_row_ids = tuple(
        str(row_id) for row_id in dropped_missing_fraction.index.tolist()
    )
    imputed_rows = (
        tuple(
            RowImputationRecord(
                row_id=str(row_id),
                imputed_columns=tuple(
                    str(column_name)
                    for column_name in imputed.columns[
                        imputed_mask.loc[row_id]
                    ].tolist()
                ),
                imputed_cell_count=int(imputed_mask.loc[row_id].sum()),
            )
            for row_id in imputed.index[
                imputed_mask.any(axis=1).to_numpy(dtype=bool, copy=False)
            ]
        )
        if not imputed.empty
        else ()
    )
    dropped_rows_missing_fraction = tuple(
        (str(row_id), float(missing_fraction_value))
        for row_id, missing_fraction_value in dropped_missing_fraction.items()
    )
    return KnnPolicyOutcome(
        phospho=imputed,
        site_metadata=filtered_site_metadata,
        k=k_value,
        distance=distance_value,
        max_missing_fraction_per_row=max_missing_fraction_value,
        dropped_row_ids=dropped_row_ids,
        dropped_rows_missing_fraction=dropped_rows_missing_fraction,
        imputed_cell_count=imputed_cell_count,
        imputed_row_ids=imputed_row_ids,
        imputed_column_ids=imputed_column_ids,
        output_missing_cell_count=output_missing_cell_count,
        rows_not_imputable=dropped_row_ids,
        imputed_rows=imputed_rows,
    )
