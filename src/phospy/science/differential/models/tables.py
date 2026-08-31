"""Differential result table contracts and validation helpers."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.validation import (
    require_dataframe,
    require_finite_numeric_dataframe,
    require_non_empty_dataframe,
    require_numeric_dataframe,
    require_string_index,
    require_unique_columns,
    require_unique_index,
)
from phospy.science.sites.identity_contracts import RESULT_IDENTITY_COLUMNS
from phospy.science.sites.identity_rules.result_identity import (
    enforce_result_table_identity_contract,
)

_RESULT_STATISTIC_COLUMNS: tuple[str, ...] = ("logFC", "t", "P.Value", "adj.P.Val")
_PUBLIC_RESULT_IDENTITY_COLUMNS: tuple[str, ...] = RESULT_IDENTITY_COLUMNS

DIFFERENTIAL_RESULT_STATUS_COLUMN = "result_status"
DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN = "result_status_reason"
DIFFERENTIAL_RESULT_STATUS_TESTED = "tested"
DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT = "withheld_all_constant"
DIFFERENTIAL_RESULT_STATUS_WITHHELD_INVALID_NUMERIC_VALUES = (
    "withheld_invalid_numeric_values"
)
DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER = "withheld_other"
DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION = "withheld_high_imputation"
DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED = (
    "withheld_insufficient_observed_values"
)
DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE = (
    "withheld_protein_preparation_ineligible"
)
DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID = (
    "withheld_protein_covariate_invalid"
)
DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID = (
    "withheld_protein_augmented_design_invalid"
)
DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE = (
    "withheld_protein_contrast_non_estimable"
)
DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK = "protein_preparation_fallback"
DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED = "protein_preparation_excluded"
DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE = "protein_covariate_non_finite"
DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE = (
    "protein_covariate_zero_variance"
)
DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT = (
    "protein_augmented_design_rank_deficient"
)
DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED = (
    "protein_augmented_design_ill_conditioned"
)
DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF = (
    "protein_augmented_design_non_positive_residual_dof"
)
DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE = (
    "protein_contrast_non_estimable"
)
DIFFERENTIAL_RESULT_PROTEIN_AWARE_STATUS_REASONS: dict[str, tuple[str, ...]] = {
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE: (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED,
    ),
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID: (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE,
    ),
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID: (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED,
        DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF,
    ),
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE: (
        DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE,
    ),
}
_DIFFERENTIAL_RESULT_PROTEIN_AWARE_REASON_CODES = frozenset(
    reason
    for reasons in DIFFERENTIAL_RESULT_PROTEIN_AWARE_STATUS_REASONS.values()
    for reason in reasons
)
DIFFERENTIAL_RESULT_WITHHELD_STATUSES: tuple[str, ...] = (
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INVALID_NUMERIC_VALUES,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID,
    DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE,
)
DIFFERENTIAL_IMPUTATION_RESULT_COLUMNS: tuple[str, ...] = (
    "imputed_cell_count",
    "observed_cell_count",
    "imputed_fraction",
    "imputation_policy",
    "imputation_fraction_threshold",
    "contains_imputed_cells",
    "observed_only_fit",
    "residual_df_adjusted_for_imputation",
    "inferential_status",
    DIFFERENTIAL_RESULT_STATUS_COLUMN,
    DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN,
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
    enforce_result_table_identity_contract(
        table=table,
        field_name=field_name,
        error_type=PhosPyInputError,
        context_label="Differential result identity metadata",
        identity_columns=_PUBLIC_RESULT_IDENTITY_COLUMNS,
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
    stat_table = table[list(_RESULT_STATISTIC_COLUMNS)]
    require_numeric_dataframe(
        stat_table,
        field_name=field_name,
        error_type=PhosPyInputError,
    )
    if (
        allow_imputation_withheld_status
        and DIFFERENTIAL_RESULT_STATUS_COLUMN in table.columns
    ):
        _validate_status_statistics(
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


def _validate_status_statistics(
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
    if DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN not in table.columns:
        raise PhosPyInputError(
            f"{field_name} rows with {DIFFERENTIAL_RESULT_STATUS_COLUMN} must include "
            f"{DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN}"
        )
    validate_result_status_reason_contract(table, field_name=field_name)

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
    reason_column = table[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN]
    empty_reason_mask = np.asarray(
        [str(value).strip() == "" for value in reason_column.tolist()],
        dtype=bool,
    )
    invalid_reason_positions = np.flatnonzero(withheld_mask & empty_reason_mask)
    if int(invalid_reason_positions.size):
        invalid_labels = [
            str(table.index[int(position)]) for position in invalid_reason_positions[:3]
        ]
        suffix = (
            ""
            if int(invalid_reason_positions.size) <= 3
            else f", +{int(invalid_reason_positions.size - 3)} more"
        )
        raise PhosPyInputError(
            f"{field_name} withheld rows must include a non-empty "
            f"{DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN}; invalid rows: "
            + ", ".join(invalid_labels)
            + suffix
        )
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
            f"{field_name} withheld rows must contain missing "
            "values for logFC, t, P.Value, and adj.P.Val; invalid values: "
            + ", ".join(previews)
            + suffix
        )


def validate_result_status_reason_contract(
    table: pd.DataFrame,
    *,
    field_name: str,
) -> None:
    """Validate differential result status values and protein-aware reason codes."""

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
    if DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN not in table.columns:
        raise PhosPyInputError(
            f"{field_name} rows with {DIFFERENTIAL_RESULT_STATUS_COLUMN} must include "
            f"{DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN}"
        )
    _validate_protein_aware_status_reason_codes(
        table=table,
        status_values=status_values,
        field_name=field_name,
    )


def _validate_protein_aware_status_reason_codes(
    *,
    table: pd.DataFrame,
    status_values: pd.Series,
    field_name: str,
) -> None:
    reason_values = table[DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN].astype(str)
    invalid_rows: list[str] = []
    for row_label, status, reason in zip(
        table.index.tolist(),
        status_values.tolist(),
        reason_values.tolist(),
        strict=True,
    ):
        allowed_reasons = DIFFERENTIAL_RESULT_PROTEIN_AWARE_STATUS_REASONS.get(status)
        if allowed_reasons is None:
            if reason not in _DIFFERENTIAL_RESULT_PROTEIN_AWARE_REASON_CODES:
                continue
        elif reason in allowed_reasons:
            continue
        invalid_rows.append(f"{row_label!r}: status={status!r}, reason={reason!r}")
    if not invalid_rows:
        return
    preview = ", ".join(invalid_rows[:3])
    suffix = "" if len(invalid_rows) <= 3 else f", +{len(invalid_rows) - 3} more"
    raise PhosPyInputError(
        f"{field_name}.{DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN} contains "
        "unsupported protein-aware reason codes: "
        f"{preview}{suffix}"
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
    "DIFFERENTIAL_RESULT_PROTEIN_AWARE_STATUS_REASONS",
    "DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_ILL_CONDITIONED",
    "DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_NON_POSITIVE_RESIDUAL_DOF",
    "DIFFERENTIAL_RESULT_REASON_PROTEIN_AUGMENTED_DESIGN_RANK_DEFICIENT",
    "DIFFERENTIAL_RESULT_REASON_PROTEIN_CONTRAST_NON_ESTIMABLE",
    "DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_NON_FINITE",
    "DIFFERENTIAL_RESULT_REASON_PROTEIN_COVARIATE_ZERO_VARIANCE",
    "DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_EXCLUDED",
    "DIFFERENTIAL_RESULT_REASON_PROTEIN_PREPARATION_FALLBACK",
    "DIFFERENTIAL_RESULT_STATUS_COLUMN",
    "DIFFERENTIAL_RESULT_STATUS_REASON_COLUMN",
    "DIFFERENTIAL_RESULT_STATUS_TESTED",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_ALL_CONSTANT",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_HIGH_IMPUTATION",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_INSUFFICIENT_OBSERVED",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_INVALID_NUMERIC_VALUES",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_OTHER",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_AUGMENTED_DESIGN_INVALID",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_CONTRAST_NON_ESTIMABLE",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_COVARIATE_INVALID",
    "DIFFERENTIAL_RESULT_STATUS_WITHHELD_PROTEIN_PREPARATION_INELIGIBLE",
    "DIFFERENTIAL_RESULT_WITHHELD_STATUSES",
    "validate_result_status_reason_contract",
]
