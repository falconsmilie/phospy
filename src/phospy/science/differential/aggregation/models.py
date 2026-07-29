"""Typed peptide-to-site differential estimate-combination models.

The preferred PhosPy peptide-to-site lane resolves peptide evidence to a
site-level sample-intensity matrix before running the core differential model.
The models in this module cover the narrower post-hoc lane for caller-supplied
peptide-level differential estimates. That lane is supported only when the
input uncertainty is typed and dependence assumptions are explicit.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import cast

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.frames.ownership import export_dataframe, own_dataframe
from phospy.frames.validation import (
    require_columns,
    require_dataframe,
    require_non_empty_dataframe,
    require_unique_columns,
)
from phospy.provenance.scientific_policy_models import ScientificPolicyRecord
from phospy.science.configs.differential import (
    MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG,
    SUPPORTED_MULTIPLE_TESTING_METHODS,
    MultipleTestingMethod,
)

PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS = "supported_typed_estimate_combination_v1"

PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SAMPLE_INTENSITY = (
    "sample_intensity_resolution_before_differential_model"
)
PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE = "single_estimate_passthrough"
PEPTIDE_TO_SITE_AGGREGATION_LEVEL_POSTHOC = "posthoc_independent_estimate_combination"

PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH = (
    "single_estimate_passthrough"
)
PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P = "stouffer_signed_p_independent"
PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE = (
    "fixed_effect_inverse_variance_independent"
)
SUPPORTED_PEPTIDE_TO_SITE_UNCERTAINTY_METHODS: tuple[str, ...] = (
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH,
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P,
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE,
)

PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES = "independent_sources"
PEPTIDE_TO_SITE_DEPENDENCE_POLICY_SAME_EXPERIMENT_CORRELATED = (
    "same_experiment_correlated"
)
SUPPORTED_PEPTIDE_TO_SITE_ESTIMATE_DEPENDENCE_POLICIES: tuple[str, ...] = (
    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES,
    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_SAME_EXPERIMENT_CORRELATED,
)
SUPPORTED_PEPTIDE_TO_SITE_CONFIG_DEPENDENCE_POLICIES: tuple[str, ...] = (
    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES,
)

PEPTIDE_TO_SITE_MAPPING_POLICY_EXPLICIT_SITE_ID = "explicit_site_id"
PEPTIDE_TO_SITE_MAPPING_POLICY_KEEP_JOINT = "keep_joint"
PEPTIDE_TO_SITE_MAPPING_POLICY_SPLIT_EQUAL_WEIGHT = "split_equal_weight"
PEPTIDE_TO_SITE_MAPPING_POLICY_EXCLUDE_FROM_STATISTICAL_MODEL = (
    "exclude_from_statistical_model"
)
SUPPORTED_PEPTIDE_TO_SITE_MAPPING_POLICIES: tuple[str, ...] = (
    PEPTIDE_TO_SITE_MAPPING_POLICY_EXPLICIT_SITE_ID,
    PEPTIDE_TO_SITE_MAPPING_POLICY_KEEP_JOINT,
    PEPTIDE_TO_SITE_MAPPING_POLICY_SPLIT_EQUAL_WEIGHT,
    PEPTIDE_TO_SITE_MAPPING_POLICY_EXCLUDE_FROM_STATISTICAL_MODEL,
)

PEPTIDE_DIFFERENTIAL_ESTIMATE_REQUIRED_COLUMNS: tuple[str, ...] = (
    "site_id",
    "peptide_id",
    "effect",
    "standard_error",
    "statistic",
    "p_value",
    "residual_degrees_of_freedom",
    "moderated_degrees_of_freedom",
    "source_experiment_id",
    "dependence_policy",
    "peptide_to_site_mapping_policy",
)
PEPTIDE_DIFFERENTIAL_ESTIMATE_OPTIONAL_COLUMNS: tuple[str, ...] = (
    "mapping_weight",
    "mapping_uncertainty",
)

PEPTIDE_TO_SITE_AGGREGATION_RESULT_REQUIRED_COLUMNS: tuple[str, ...] = (
    "logFC",
    "standard_error",
    "uncertainty_statistic",
    "P.Value",
    "adj.P.Val",
    "residual_degrees_of_freedom",
    "moderated_degrees_of_freedom",
    "n_peptide_observations",
    "n_peptides_used",
    "source_experiment_ids",
    "peptide_to_site_mapping_policy",
    "multi_site_estimate_count",
    "aggregation_level",
    "dependence_assumption",
    "uncertainty_method",
    "correction_method",
    "p_value_method",
    "statistic_distribution",
)


@dataclass(frozen=True, slots=True, init=False)
class PeptideDifferentialEstimateTable:
    """Typed peptide-level differential estimates for post-hoc site combination.

    Required columns record the original finite-degree-of-freedom uncertainty.
    ``p_value`` is the original two-sided p-value from the peptide-level model;
    single-estimate site outputs pass it through directly rather than converting
    the statistic to a normal approximation.
    """

    _frame: pd.DataFrame

    def __init__(
        self,
        frame: pd.DataFrame,
    ) -> None:
        frame = own_dataframe(
            frame,
            field_name="peptide_differential_estimate_table",
            error_type=PhosPyInputError,
        )
        require_dataframe(
            frame,
            field_name="peptide_differential_estimate_table",
            allow_empty=False,
            error_type=PhosPyInputError,
        )
        require_non_empty_dataframe(
            frame,
            field_name="peptide_differential_estimate_table",
            error_type=PhosPyInputError,
        )
        require_unique_columns(
            frame,
            field_name="peptide_differential_estimate_table",
            error_type=PhosPyInputError,
        )
        require_columns(
            frame,
            field_name="peptide_differential_estimate_table",
            required_columns=PEPTIDE_DIFFERENTIAL_ESTIMATE_REQUIRED_COLUMNS,
            error_type=PhosPyInputError,
        )

        selected_columns = [
            *PEPTIDE_DIFFERENTIAL_ESTIMATE_REQUIRED_COLUMNS,
            *(
                column
                for column in PEPTIDE_DIFFERENTIAL_ESTIMATE_OPTIONAL_COLUMNS
                if column in frame.columns
            ),
        ]
        canonical = _copy_selected_columns(frame, selected_columns)
        for column_name in (
            "site_id",
            "peptide_id",
            "source_experiment_id",
        ):
            _set_column_values(
                canonical,
                column_name,
                tuple(
                    _canonical_non_empty_string(
                        value,
                        field_name=(
                            f"peptide_differential_estimate_table.{column_name}"
                        ),
                    )
                    for value in _column_values(canonical, column_name)
                ),
            )
        _set_column_values(
            canonical,
            "peptide_to_site_mapping_policy",
            tuple(
                _canonical_mapping_policy(value)
                for value in _column_values(
                    canonical,
                    "peptide_to_site_mapping_policy",
                )
            ),
        )
        _set_column_values(
            canonical,
            "dependence_policy",
            tuple(
                _canonical_estimate_dependence_policy(value)
                for value in _column_values(canonical, "dependence_policy")
            ),
        )
        for column_name in ("effect", "statistic"):
            _set_column_values(
                canonical,
                column_name,
                tuple(
                    _canonical_finite_float(
                        value,
                        field_name=(
                            f"peptide_differential_estimate_table.{column_name}"
                        ),
                    )
                    for value in _column_values(canonical, column_name)
                ),
            )
        _set_column_values(
            canonical,
            "standard_error",
            tuple(
                _canonical_positive_finite_float(
                    value,
                    field_name="peptide_differential_estimate_table.standard_error",
                )
                for value in _column_values(canonical, "standard_error")
            ),
        )
        _set_column_values(
            canonical,
            "p_value",
            tuple(
                _canonical_p_value(
                    value,
                    field_name="peptide_differential_estimate_table.p_value",
                )
                for value in _column_values(canonical, "p_value")
            ),
        )
        for column_name in (
            "residual_degrees_of_freedom",
            "moderated_degrees_of_freedom",
        ):
            _set_column_values(
                canonical,
                column_name,
                tuple(
                    _canonical_positive_finite_float(
                        value,
                        field_name=(
                            f"peptide_differential_estimate_table.{column_name}"
                        ),
                    )
                    for value in _column_values(canonical, column_name)
                ),
            )
        if "mapping_weight" in canonical.columns:
            _set_column_values(
                canonical,
                "mapping_weight",
                tuple(
                    _canonical_positive_finite_float(
                        value,
                        field_name=(
                            "peptide_differential_estimate_table.mapping_weight"
                        ),
                    )
                    for value in _column_values(canonical, "mapping_weight")
                ),
            )
        if "mapping_uncertainty" in canonical.columns:
            _set_column_values(
                canonical,
                "mapping_uncertainty",
                tuple(
                    _canonical_bool(
                        value,
                        field_name=(
                            "peptide_differential_estimate_table.mapping_uncertainty"
                        ),
                    )
                    for value in _column_values(canonical, "mapping_uncertainty")
                ),
            )
        object.__setattr__(self, "_frame", canonical)

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self._frame)

    @property
    def site_ids(self) -> tuple[str, ...]:
        return tuple(
            dict.fromkeys(
                str(value) for value in _column_values(self._frame, "site_id")
            )
        )


@dataclass(frozen=True, slots=True)
class PeptideToSiteAggregationConfig:
    """Configuration for typed peptide-to-site estimate combination."""

    uncertainty_method: str = PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P
    min_estimates_per_site: int = 1
    dependence_policy: str = PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES
    multiple_testing_method: MultipleTestingMethod = (
        MULTIPLE_TESTING_METHOD_BENJAMINI_HOCHBERG
    )

    def __post_init__(self) -> None:
        if self.uncertainty_method not in SUPPORTED_PEPTIDE_TO_SITE_UNCERTAINTY_METHODS:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_PEPTIDE_TO_SITE_UNCERTAINTY_METHODS
            )
            raise PhosPyInputError(
                "peptide_to_site_aggregation.uncertainty_method must be one of: "
                f"{supported}"
            )
        min_estimates = cast(object, self.min_estimates_per_site)
        if not isinstance(min_estimates, int) or min_estimates < 1:
            raise PhosPyInputError(
                "peptide_to_site_aggregation.min_estimates_per_site must be an int >= 1"
            )
        if (
            self.dependence_policy
            not in SUPPORTED_PEPTIDE_TO_SITE_CONFIG_DEPENDENCE_POLICIES
        ):
            supported = ", ".join(
                repr(value)
                for value in SUPPORTED_PEPTIDE_TO_SITE_CONFIG_DEPENDENCE_POLICIES
            )
            raise PhosPyInputError(
                "peptide_to_site_aggregation.dependence_policy must be one of: "
                f"{supported}"
            )
        if self.multiple_testing_method not in SUPPORTED_MULTIPLE_TESTING_METHODS:
            supported = ", ".join(
                repr(value) for value in SUPPORTED_MULTIPLE_TESTING_METHODS
            )
            raise PhosPyInputError(
                "peptide_to_site_aggregation.multiple_testing_method must be one of: "
                f"{supported}"
            )
        object.__setattr__(
            self,
            "uncertainty_method",
            str(self.uncertainty_method),
        )
        object.__setattr__(
            self,
            "dependence_policy",
            str(self.dependence_policy),
        )


@dataclass(frozen=True, slots=True, init=False)
class PeptideToSiteAggregationResult:
    """Site-level result from typed peptide estimate combination."""

    contrast_name: str
    table: pd.DataFrame
    warnings: tuple[str, ...]
    provenance: Mapping[str, object]
    scientific_policies: tuple[ScientificPolicyRecord, ...]

    def __init__(
        self,
        *,
        contrast_name: str,
        table: pd.DataFrame,
        warnings: Sequence[str] = (),
        provenance: Mapping[str, object] | None = None,
        scientific_policies: tuple[ScientificPolicyRecord, ...] = (),
    ) -> None:
        contrast_name_value = _canonical_non_empty_string(
            contrast_name,
            field_name="peptide_to_site_aggregation_result.contrast_name",
        )
        table = own_dataframe(
            table,
            field_name="peptide_to_site_aggregation_result.table",
            error_type=PhosPyInputError,
        )
        require_dataframe(
            table,
            field_name="peptide_to_site_aggregation_result.table",
            allow_empty=False,
            error_type=PhosPyInputError,
        )
        require_non_empty_dataframe(
            table,
            field_name="peptide_to_site_aggregation_result.table",
            error_type=PhosPyInputError,
        )
        require_unique_columns(
            table,
            field_name="peptide_to_site_aggregation_result.table",
            error_type=PhosPyInputError,
        )
        require_columns(
            table,
            field_name="peptide_to_site_aggregation_result.table",
            required_columns=PEPTIDE_TO_SITE_AGGREGATION_RESULT_REQUIRED_COLUMNS,
            error_type=PhosPyInputError,
        )
        _validate_result_p_values(table)
        warnings_tuple = tuple(str(value) for value in warnings)
        for policy in scientific_policies:
            if not isinstance(cast(object, policy), ScientificPolicyRecord):
                raise PhosPyInputError(
                    "peptide_to_site_aggregation_result.scientific_policies must "
                    "contain ScientificPolicyRecord values"
                )
        object.__setattr__(self, "contrast_name", contrast_name_value)
        object.__setattr__(self, "table", table)
        object.__setattr__(self, "warnings", warnings_tuple)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {}
                if provenance is None
                else {str(key): value for key, value in provenance.items()}
            ),
        )
        object.__setattr__(self, "scientific_policies", tuple(scientific_policies))

    def to_dataframe(self) -> pd.DataFrame:
        return export_dataframe(self.table)


def _copy_selected_columns(
    frame: pd.DataFrame,
    columns: Sequence[str],
) -> pd.DataFrame:
    selected = cast(pd.DataFrame, frame.loc[:, list(columns)])  # pyright: ignore[reportUnknownMemberType, reportArgumentType, reportCallIssue] - pandas-stubs cannot express generic DataFrame.loc column selection.
    copied = cast(pd.DataFrame, selected.copy(deep=True))  # pyright: ignore[reportUnknownMemberType] - pandas-stubs loses concrete DataFrame type for copy.
    return copied


def _column_values(frame: pd.DataFrame, column_name: str) -> tuple[object, ...]:
    column: pd.Series = frame.loc[:, column_name]  # pyright: ignore[reportUnknownMemberType, reportArgumentType, reportCallIssue] - pandas-stubs cannot express scalar-column DataFrame.loc selection.
    array = column.to_numpy(dtype=object)  # pyright: ignore[reportUnknownMemberType] - pandas-stubs loses Series.to_numpy return type on unparameterized Series.
    return tuple(cast(Sequence[object], array))


def _set_column_values(
    frame: pd.DataFrame,
    column_name: str,
    values: Sequence[object],
) -> None:
    frame.loc[:, column_name] = list(values)  # pyright: ignore[reportUnknownMemberType, reportArgumentType, reportCallIssue] - pandas-stubs cannot express runtime-valid DataFrame.loc assignment.


def _canonical_estimate_dependence_policy(value: object) -> str:
    text = _canonical_non_empty_string(
        value,
        field_name="peptide_differential_estimate_table.dependence_policy",
    )
    if text not in SUPPORTED_PEPTIDE_TO_SITE_ESTIMATE_DEPENDENCE_POLICIES:
        supported = ", ".join(
            repr(policy)
            for policy in SUPPORTED_PEPTIDE_TO_SITE_ESTIMATE_DEPENDENCE_POLICIES
        )
        raise PhosPyInputError(
            "peptide_differential_estimate_table.dependence_policy must be one of: "
            f"{supported}"
        )
    return text


def _canonical_mapping_policy(value: object) -> str:
    text = _canonical_non_empty_string(
        value,
        field_name="peptide_differential_estimate_table.peptide_to_site_mapping_policy",
    )
    if text not in SUPPORTED_PEPTIDE_TO_SITE_MAPPING_POLICIES:
        supported = ", ".join(
            repr(policy) for policy in SUPPORTED_PEPTIDE_TO_SITE_MAPPING_POLICIES
        )
        raise PhosPyInputError(
            "peptide_differential_estimate_table."
            "peptide_to_site_mapping_policy must be one of: "
            f"{supported}"
        )
    return text


def _canonical_non_empty_string(value: object, *, field_name: str) -> str:
    if _is_missing(value):
        raise PhosPyInputError(f"{field_name} must not contain missing values")
    if not isinstance(value, str):
        raise PhosPyInputError(f"{field_name} must contain non-empty string values")
    stripped = value.strip()
    if stripped == "":
        raise PhosPyInputError(f"{field_name} must contain non-empty string values")
    return stripped


def _canonical_bool(value: object, *, field_name: str) -> bool:
    if isinstance(value, bool):
        return value
    if type(value).__name__ == "bool_":
        return bool(value)
    raise PhosPyInputError(f"{field_name} must contain bool values")


def _canonical_finite_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric")
    resolved = float(value)
    if not math.isfinite(resolved):
        raise PhosPyInputError(f"{field_name} must be finite")
    return resolved


def _canonical_positive_finite_float(value: object, *, field_name: str) -> float:
    resolved = _canonical_finite_float(value, field_name=field_name)
    if resolved <= 0.0:
        raise PhosPyInputError(f"{field_name} must be > 0.0")
    return resolved


def _canonical_p_value(value: object, *, field_name: str) -> float:
    resolved = _canonical_finite_float(value, field_name=field_name)
    if resolved < 0.0 or resolved > 1.0:
        raise PhosPyInputError(f"{field_name} must be within [0.0, 1.0]")
    return resolved


def _is_missing(value: object) -> bool:
    try:
        return bool(pd.Series([value], dtype="object").isna().iloc[0])
    except (TypeError, ValueError):
        return False


def _validate_result_p_values(table: pd.DataFrame) -> None:
    for column_name in ("P.Value", "adj.P.Val"):
        for value in _column_values(table, column_name):
            if isinstance(value, bool) or not isinstance(value, int | float):
                raise PhosPyInputError(
                    "peptide_to_site_aggregation_result.table "
                    f"{column_name} values must be numeric p-values or NaN"
                )
            number = float(value)
            if math.isnan(number):
                continue
            if not math.isfinite(number) or number < 0.0 or number > 1.0:
                raise PhosPyInputError(
                    "peptide_to_site_aggregation_result.table "
                    f"{column_name} values must be within [0.0, 1.0]"
                )


__all__ = [
    "PEPTIDE_DIFFERENTIAL_ESTIMATE_OPTIONAL_COLUMNS",
    "PEPTIDE_DIFFERENTIAL_ESTIMATE_REQUIRED_COLUMNS",
    "PEPTIDE_TO_SITE_AGGREGATION_LEVEL_POSTHOC",
    "PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SAMPLE_INTENSITY",
    "PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE",
    "PEPTIDE_TO_SITE_AGGREGATION_RESULT_REQUIRED_COLUMNS",
    "PEPTIDE_TO_SITE_AGGREGATION_SUPPORT_STATUS",
    "PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES",
    "PEPTIDE_TO_SITE_DEPENDENCE_POLICY_SAME_EXPERIMENT_CORRELATED",
    "PEPTIDE_TO_SITE_MAPPING_POLICY_EXCLUDE_FROM_STATISTICAL_MODEL",
    "PEPTIDE_TO_SITE_MAPPING_POLICY_EXPLICIT_SITE_ID",
    "PEPTIDE_TO_SITE_MAPPING_POLICY_KEEP_JOINT",
    "PEPTIDE_TO_SITE_MAPPING_POLICY_SPLIT_EQUAL_WEIGHT",
    "PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE",
    "PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH",
    "PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P",
    "SUPPORTED_PEPTIDE_TO_SITE_CONFIG_DEPENDENCE_POLICIES",
    "SUPPORTED_PEPTIDE_TO_SITE_ESTIMATE_DEPENDENCE_POLICIES",
    "SUPPORTED_PEPTIDE_TO_SITE_MAPPING_POLICIES",
    "SUPPORTED_PEPTIDE_TO_SITE_UNCERTAINTY_METHODS",
    "PeptideDifferentialEstimateTable",
    "PeptideToSiteAggregationConfig",
    "PeptideToSiteAggregationResult",
]
