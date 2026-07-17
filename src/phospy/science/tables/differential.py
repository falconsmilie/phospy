"""Reporting helpers for differential result tables."""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd

from phospy.science.tables._differential_validation import (
    require_boolean,
    require_column_name,
    require_differential_result_columns,
    require_na_position,
    require_non_negative_threshold,
    require_numeric_result_column,
    require_probability_threshold,
)

ADJUSTED_P_VALUE_COLUMN = "adj.P.Val"
RAW_P_VALUE_COLUMN = "P.Value"
LOG_FOLD_CHANGE_COLUMN = "logFC"


def filter_differential_results(
    table: pd.DataFrame,
    *,
    adjusted_p_value_max: float | None = None,
    p_value_max: float | None = None,
    min_abs_effect_size: float | None = None,
    adjusted_p_value_column: str = ADJUSTED_P_VALUE_COLUMN,
    p_value_column: str = RAW_P_VALUE_COLUMN,
    effect_size_column: str = LOG_FOLD_CHANGE_COLUMN,
) -> pd.DataFrame:
    """Filter an existing differential result table without changing it."""

    frame = pd.DataFrame(
        require_differential_result_columns(
            table,
            columns=(),
            field_name="differential result table",
        ),
        copy=True,
    )
    mask = pd.Series(True, index=frame.index, dtype=bool)

    if adjusted_p_value_max is not None:
        threshold = require_probability_threshold(
            adjusted_p_value_max,
            field_name="adjusted_p_value_max",
        )
        column_name = require_column_name(
            adjusted_p_value_column,
            field_name="adjusted_p_value_column",
        )
        values = require_numeric_result_column(
            frame,
            column_name=column_name,
            field_name="differential result table",
        )
        mask &= values <= threshold

    if p_value_max is not None:
        threshold = require_probability_threshold(
            p_value_max,
            field_name="p_value_max",
        )
        column_name = require_column_name(
            p_value_column,
            field_name="p_value_column",
        )
        values = require_numeric_result_column(
            frame,
            column_name=column_name,
            field_name="differential result table",
        )
        mask &= values <= threshold

    if min_abs_effect_size is not None:
        threshold = require_non_negative_threshold(
            min_abs_effect_size,
            field_name="min_abs_effect_size",
        )
        column_name = require_column_name(
            effect_size_column,
            field_name="effect_size_column",
        )
        values = require_numeric_result_column(
            frame,
            column_name=column_name,
            field_name="differential result table",
        )
        mask &= values.abs() >= threshold

    return pd.DataFrame(frame.loc[mask, :], copy=True)


def rank_differential_results(
    table: pd.DataFrame,
    *,
    by: str,
    ascending: bool = True,
    absolute: bool = False,
    na_position: Literal["first", "last"] = "last",
) -> pd.DataFrame:
    """Rank an existing differential result table by one numeric column."""

    frame = pd.DataFrame(
        require_differential_result_columns(
            table,
            columns=(),
            field_name="differential result table",
        ),
        copy=True,
    )
    column_name = require_column_name(by, field_name="by")
    ascending_value = require_boolean(ascending, field_name="ascending")
    absolute_value = require_boolean(absolute, field_name="absolute")
    na_position_value = require_na_position(na_position, field_name="na_position")
    values = require_numeric_result_column(
        frame,
        column_name=column_name,
        field_name="differential result table",
    )
    sort_key = values.abs() if absolute_value else values
    sort_frame = pd.DataFrame(
        {
            "sort_key": sort_key.to_numpy(dtype="float64", copy=False),
            "original_order": np.arange(int(frame.shape[0])),
        }
    )
    sort_frame = sort_frame.sort_values(
        by=["sort_key", "original_order"],
        ascending=[ascending_value, True],
        na_position=na_position_value,
        kind="mergesort",
    )
    return pd.DataFrame(frame.iloc[sort_frame.index.to_numpy()], copy=True)


__all__ = [
    "ADJUSTED_P_VALUE_COLUMN",
    "LOG_FOLD_CHANGE_COLUMN",
    "RAW_P_VALUE_COLUMN",
    "filter_differential_results",
    "rank_differential_results",
]
