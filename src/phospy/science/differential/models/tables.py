"""Differential result table contracts and validation helpers."""
# pyright: reportMissingTypeStubs=false

from __future__ import annotations

from typing import cast

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.validation.common.dataframes import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_string_index,
    require_unique_columns,
    require_unique_index,
)
from phospy.validation.identity_contracts import (
    RESULT_IDENTITY_COLUMNS,
    RESULT_TABLE_IDENTITY_CONTRACT,
    enforce_phosphosite_identity_contract,
    enforce_required_identity_text_columns,
    enforce_result_identity_metadata_coherence,
)

_RESULT_STATISTIC_COLUMNS: tuple[str, ...] = ("logFC", "t", "P.Value", "adj.P.Val")
_PUBLIC_RESULT_IDENTITY_COLUMNS: tuple[str, ...] = RESULT_IDENTITY_COLUMNS

DIFFERENTIAL_RESULT_STATUS_COLUMN = "result_status"
DIFFERENTIAL_RESULT_STATUS_TESTED = "tested"
DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION = "withheld_high_imputation"
DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED = (
    "withheld_insufficient_observed_samples"
)
DIFFERENTIAL_RESULT_WITHHELD_STATUSES: tuple[str, ...] = (
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
)
DIFFERENTIAL_IMPUTATION_RESULT_COLUMNS: tuple[str, ...] = (
    "imputed_cell_count",
    "observed_cell_count",
    "imputed_fraction",
    "imputation_policy",
    "imputation_fraction_threshold",
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
)


def validate_result_table_contract(
    table: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    _validate_result_table_statistics(
        table=table,
        field_name=field_name,
        allow_imputation_withheld_status=True,
    )
    enforce_phosphosite_identity_contract(
        site_metadata=table,
        field_name=field_name,
        contract=RESULT_TABLE_IDENTITY_CONTRACT,
        error_type=PhosPyInputError,
        compare_raw_site_key_column_before_decode=True,
    )
    enforce_required_identity_text_columns(
        table=table,
        field_name=field_name,
        columns=_PUBLIC_RESULT_IDENTITY_COLUMNS,
        error_type=PhosPyInputError,
    )
    enforce_result_identity_metadata_coherence(
        table=table,
        field_name=field_name,
        context_label="Differential result identity metadata",
        error_type=PhosPyInputError,
    )


def validate_computation_result_table_contract(
    table: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    _validate_result_table_statistics(
        table=table,
        field_name=field_name,
        allow_imputation_withheld_status=False,
    )
    present_identity = [
        column for column in _PUBLIC_RESULT_IDENTITY_COLUMNS if column in table.columns
    ]
    if present_identity:
        joined = ", ".join(present_identity)
        raise PhosPyInputError(
            f"{field_name} must be stat-only and must not include identity columns: "
            f"{joined}"
        )


def _validate_result_table_statistics(
    *,
    table: pd.DataFrame,
    field_name: str,
    allow_imputation_withheld_status: bool,
) -> None:
    require_dataframe(
        table,
        field_name=field_name,
        allow_empty=False,
        error_type=PhosPyInputError,
    )
    require_non_empty_dataframe(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_string_index(
        table.index,
        field_name=f"{field_name}.index",
        error_type=PhosPyInputError,
    )
    require_unique_index(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    require_unique_columns(
        table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    missing = [
        column for column in _RESULT_STATISTIC_COLUMNS if column not in table.columns
    ]
    if missing:
        joined = ", ".join(missing)
        raise PhosPyInputError(f"{field_name} is missing required columns: {joined}")
    stat_table = cast(pd.DataFrame, table[list(_RESULT_STATISTIC_COLUMNS)])
    require_numeric_dataframe(
        stat_table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    if (
        allow_imputation_withheld_status
        and DIFFERENTIAL_RESULT_STATUS_COLUMN in table.columns
    ):
        _validate_imputation_status_statistics(
            table=table,
            stat_table=stat_table,
            field_name=field_name,
        )
    else:
        require_finite_numeric_dataframe(
            stat_table,
            field_name=field_name,
            error_type=PhosPyInputError,
            allow_missing=False,
        )
    _validate_unit_interval_column(
        table=table,
        column_name="P.Value",
        field_name=field_name,
    )
    _validate_unit_interval_column(
        table=table,
        column_name="adj.P.Val",
        field_name=field_name,
    )


def _validate_imputation_status_statistics(
    *,
    table: pd.DataFrame,
    stat_table: pd.DataFrame,
    field_name: str,
) -> None:
    status_values = table[DIFFERENTIAL_RESULT_STATUS_COLUMN].astype(str)
    allowed_statuses = {
        DIFFERENTIAL_RESULT_STATUS_TESTED,
        *DIFFERENTIAL_RESULT_WITHHELD_STATUSES,
    }
    unknown_statuses = sorted(set(status_values.tolist()) - allowed_statuses)
    if unknown_statuses:
        raise PhosPyInputError(
            f"{field_name}.{DIFFERENTIAL_RESULT_STATUS_COLUMN} contains unsupported "
            "values: " + ", ".join(repr(value) for value in unknown_statuses)
        )

    status_array = status_values.to_numpy(dtype=str)
    tested_mask: npt.NDArray[np.bool_] = np.asarray(
        status_array == DIFFERENTIAL_RESULT_STATUS_TESTED,
        dtype=bool,
    )
    withheld_mask: npt.NDArray[np.bool_] = np.isin(
        status_array,
        DIFFERENTIAL_RESULT_WITHHELD_STATUSES,
    )
    stat_values_float: npt.NDArray[np.float64] = np.asarray(
        stat_table.to_numpy(dtype=float),
        dtype=np.float64,
    )
    invalid_tested_mask = ~np.isfinite(stat_values_float[tested_mask, :])
    if bool(invalid_tested_mask.any()):
        raise PhosPyInputError(
            f"{field_name} rows with "
            f"{DIFFERENTIAL_RESULT_STATUS_COLUMN}="
            f"{DIFFERENTIAL_RESULT_STATUS_TESTED!r} must contain finite "
            "numeric logFC, t, P.Value, and adj.P.Val values"
        )

    withheld_row_positions = np.flatnonzero(withheld_mask)
    if not int(withheld_row_positions.size):
        return
    withheld_values: npt.NDArray[np.object_] = np.asarray(
        stat_table.to_numpy(dtype=object)[withheld_mask, :],
        dtype=object,
    )
    non_missing_mask: npt.NDArray[np.bool_] = np.asarray(
        ~pd.isna(withheld_values),
        dtype=bool,
    )
    if bool(non_missing_mask.any()):
        invalid_positions = np.argwhere(non_missing_mask)
        previews: list[str] = []
        for row_position, column_position in invalid_positions[:3]:
            source_row_position = int(withheld_row_positions[int(row_position)])
            previews.append(
                f"({stat_table.index[source_row_position]!r}, "
                f"{stat_table.columns[int(column_position)]!r})"
            )
        suffix = (
            ""
            if int(invalid_positions.shape[0]) <= 3
            else f", +{int(invalid_positions.shape[0] - 3)} more"
        )
        raise PhosPyInputError(
            f"{field_name} withheld imputation-policy rows must contain missing "
            "values for logFC, t, P.Value, and adj.P.Val; invalid values: "
            + ", ".join(previews)
            + suffix
        )


def _validate_unit_interval_column(
    *,
    table: pd.DataFrame,
    column_name: str,
    field_name: str,
) -> None:
    column = table[column_name]
    values = column.to_numpy(dtype=float)
    invalid_mask = (values < 0.0) | (values > 1.0)
    if not np.any(invalid_mask):
        return

    invalid_positions = np.flatnonzero(invalid_mask)
    preview: list[str] = []
    for position in invalid_positions[:3]:
        preview.append(f"({table.index[position]!r}, {values[position]:.6g})")
    suffix = (
        ""
        if invalid_positions.size <= 3
        else f", +{int(invalid_positions.size - 3)} more"
    )
    examples = ", ".join(preview)
    raise PhosPyInputError(
        f"{field_name}.{column_name} must be within [0, 1] for each feature; "
        f"invalid values: {examples}{suffix}; "
        f"invalid_entry_count={int(invalid_positions.size)}"
    )


__all__ = [
    "DIFFERENTIAL_IMPUTATION_RESULT_COLUMNS",
    "DIFFERENTIAL_RESULT_STATUS_COLUMN",
    "DIFFERENTIAL_RESULT_STATUS_TESTED",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED",
    "DIFFERENTIAL_RESULT_WITHHELD_STATUSES",
]
