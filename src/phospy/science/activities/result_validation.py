"""Activity result-table validation and export helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phospy._deprecations import warn_deprecated
from phospy.errors.validation import PhosPyValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.frames.ownership import export_dataframe
from phospy.science.activities.semantics import (
    ActivityInputSemantics,
    ActivityProfileAxis,
    ActivityProfileMetadata,
    ActivityQuantitativeSemantics,
)
from phospy.science.tables.activity import ActivityMatrix


def _resolve_activity_result_semantics(
    *,
    activity_matrix: pd.DataFrame,
    input_semantics: ActivityInputSemantics | None,
    profile_metadata: ActivityProfileMetadata | None,
) -> tuple[ActivityInputSemantics, ActivityProfileMetadata]:
    profile_ids = tuple(str(column) for column in activity_matrix.columns)
    if input_semantics is None and profile_metadata is None:
        warn_deprecated(
            "activities.result.missing_semantics",
            stacklevel=3,
        )
        input_semantics = ActivityInputSemantics(
            profile_axis=ActivityProfileAxis.SAMPLE,
            quantitative_semantics=(
                ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE
            ),
        )
        profile_metadata = ActivityProfileMetadata(
            axis=ActivityProfileAxis.SAMPLE,
            profile_ids=profile_ids,
            sample_ids=profile_ids,
        )
    elif input_semantics is None or profile_metadata is None:
        raise WorkflowBoundaryError(
            "activity_result.input_semantics and activity_result.profile_metadata "
            "must be provided together"
        )
    if not isinstance(input_semantics, ActivityInputSemantics):
        raise WorkflowBoundaryError(
            "activity_result.input_semantics must be ActivityInputSemantics"
        )
    if not isinstance(profile_metadata, ActivityProfileMetadata):
        raise WorkflowBoundaryError(
            "activity_result.profile_metadata must be ActivityProfileMetadata"
        )
    if profile_metadata.axis is not input_semantics.profile_axis:
        raise WorkflowBoundaryError(
            "activity_result.profile_metadata.axis must match "
            "activity_result.input_semantics.profile_axis"
        )
    observed_profile_ids = tuple(str(value) for value in profile_metadata.profile_ids)
    if observed_profile_ids != profile_ids:
        raise WorkflowBoundaryError(
            "activity_result.profile_metadata.profile_ids must exactly match "
            "activity_result.activity_matrix columns; "
            f"expected={profile_ids!r}, got={observed_profile_ids!r}"
        )
    return input_semantics, profile_metadata


def _validate_activity_statistics_profile_contract(
    *,
    statistics_table: pd.DataFrame,
    input_semantics: ActivityInputSemantics,
    profile_metadata: ActivityProfileMetadata,
) -> None:
    profile_values = tuple(
        str(value) for value in statistics_table.loc[:, "profile_id"].tolist()
    )
    declared_profile_ids = tuple(str(value) for value in profile_metadata.profile_ids)
    declared_profile_set = set(declared_profile_ids)
    unknown_profile_ids = sorted(set(profile_values).difference(declared_profile_set))
    if unknown_profile_ids:
        raise PhosPyValidationError(
            "activity_result.statistics_table.profile_id values must occur in "
            "activity_result.profile_metadata.profile_ids; "
            f"unknown_profile_ids={tuple(unknown_profile_ids)!r}"
        )

    namespace_ids = _activity_statistics_profile_namespace(
        profile_metadata=profile_metadata,
    )
    namespace_set = set(namespace_ids)
    outside_namespace = sorted(set(profile_values).difference(namespace_set))
    if outside_namespace:
        axis = _activity_profile_axis_value(profile_metadata.axis)
        raise PhosPyValidationError(
            "activity_result.statistics_table.profile_id values must use the "
            f"declared {axis} profile namespace; "
            f"unknown_profile_ids={tuple(outside_namespace)!r}"
        )

    if input_semantics.has_real_condition_contract:
        if "condition" not in statistics_table.columns:
            return
        condition_values = tuple(
            str(value) for value in statistics_table.loc[:, "condition"].tolist()
        )
        mismatched = sorted(
            {
                profile_id
                for profile_id, condition_id in zip(
                    profile_values,
                    condition_values,
                    strict=True,
                )
                if condition_id != profile_id
            }
        )
        if mismatched:
            raise PhosPyValidationError(
                "activity_result.statistics_table.condition must equal profile_id "
                "for condition-summary activity results; "
                f"mismatched_profile_ids={tuple(mismatched)!r}"
            )
        return

    if "condition" in statistics_table.columns:
        axis = _activity_profile_axis_value(input_semantics.profile_axis)
        raise PhosPyValidationError(
            "activity_result.statistics_table.condition is reserved for "
            "condition-summary activity results; "
            f"profile_axis={axis!r} does not define a biological condition contract"
        )


def _activity_statistics_profile_namespace(
    *,
    profile_metadata: ActivityProfileMetadata,
) -> tuple[str, ...]:
    axis = profile_metadata.axis
    if axis is ActivityProfileAxis.SAMPLE:
        return tuple(str(value) for value in profile_metadata.sample_ids)
    if axis is ActivityProfileAxis.CONDITION_SUMMARY:
        return tuple(str(value) for value in profile_metadata.condition_ids)
    if axis is ActivityProfileAxis.CONTRAST:
        return tuple(str(value) for value in profile_metadata.contrast_ids)
    return tuple(str(value) for value in profile_metadata.profile_ids)


def _activity_profile_axis_value(axis: ActivityProfileAxis | str) -> str:
    if isinstance(axis, ActivityProfileAxis):
        return axis.value
    return str(axis)


def _activity_profile_axis_name(
    input_semantics: ActivityInputSemantics,
) -> str:
    if input_semantics.has_real_condition_contract:
        return "condition"
    return "profile_id"


def _apply_activity_profile_axis_name(
    frame: pd.DataFrame,
    *,
    input_semantics: ActivityInputSemantics,
) -> pd.DataFrame:
    if input_semantics.profile_axis is ActivityProfileAxis.SAMPLE:
        return frame
    renamed = frame.copy(deep=False)
    renamed.columns = renamed.columns.copy()
    renamed.columns.name = _activity_profile_axis_name(input_semantics)
    return renamed


def _apply_optional_activity_profile_axis_name(
    frame: pd.DataFrame | None,
    *,
    input_semantics: ActivityInputSemantics,
) -> pd.DataFrame | None:
    if frame is None:
        return None
    return _apply_activity_profile_axis_name(
        frame,
        input_semantics=input_semantics,
    )


def _validate_optional_activity_matrix(
    matrix: pd.DataFrame | None,
    *,
    field_name: str,
    assume_owned: bool,
) -> pd.DataFrame | None:
    if matrix is None:
        return None
    return ActivityMatrix(
        frame=matrix,
        field_name=field_name,
        _assume_owned=assume_owned,
    ).frame


def _validate_optional_probability_matrix(
    matrix: pd.DataFrame | None,
    *,
    field_name: str,
    assume_owned: bool,
) -> pd.DataFrame | None:
    matrix = _validate_optional_activity_matrix(
        matrix,
        field_name=field_name,
        assume_owned=assume_owned,
    )
    if matrix is None:
        return None
    values = matrix.to_numpy(dtype="float64", copy=False)
    finite_mask = np.isfinite(values)
    finite_values = values[finite_mask]
    if ((finite_values < 0.0) | (finite_values > 1.0)).any():
        raise PhosPyValidationError(
            f"{field_name} must be between 0.0 and 1.0 when present"
        )
    return matrix


def _empty_activity_matrix() -> pd.DataFrame:
    return pd.DataFrame(dtype=float)


def _empty_count_matrix() -> pd.DataFrame:
    return pd.DataFrame(dtype="int64")


def _empty_count_series(name: str) -> pd.Series:
    series = pd.Series(dtype="int64", name=name)
    series.index.name = "kinase"
    return series


def _empty_target_table() -> pd.DataFrame:
    return pd.DataFrame(columns=["site_id", "kinase", "score"])


def _export_public_target_table(table: pd.DataFrame) -> pd.DataFrame:
    exported = export_dataframe(table)
    if {"site_key", "display_id"}.issubset(exported.columns):
        return exported
    legacy_columns = ["site_id", "kinase", "score"]
    if not all(column in exported.columns for column in legacy_columns):
        return exported
    return exported.loc[:, legacy_columns]


__all__ = [
    "_apply_activity_profile_axis_name",
    "_apply_optional_activity_profile_axis_name",
    "_empty_activity_matrix",
    "_empty_count_matrix",
    "_empty_count_series",
    "_empty_target_table",
    "_export_public_target_table",
    "_resolve_activity_result_semantics",
    "_validate_activity_statistics_profile_contract",
    "_validate_optional_activity_matrix",
    "_validate_optional_probability_matrix",
]
