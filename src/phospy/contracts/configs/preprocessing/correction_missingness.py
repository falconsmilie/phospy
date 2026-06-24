"""Missingness and temporary-imputation contracts for future correction.

These contracts describe intent and provenance shape only. They do not impute
values, run SPS/RUV-style correction, or convert temporary imputed values into
observed biological evidence.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import TypeAlias

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.policies import PolicyEnum
from phospy.validation.common.config_values import require_non_empty_string
from phospy.validation.common.numbers import require_optional_int_at_least

JsonScalar: TypeAlias = str | int | float | bool | None


class TemporaryImputationMethod(PolicyEnum):
    """Temporary imputation labels for future correction mechanics."""

    NONE = "none"
    ROW_MEDIAN_TEMPORARY = "row_median_temporary"
    MINPROB_TEMPORARY = "minprob_temporary"
    KNN_TEMPORARY = "knn_temporary"
    UNSUPPORTED = "unsupported"


class OriginallyMissingCellTracking(PolicyEnum):
    """How originally missing cells are tracked across correction mechanics."""

    NONE = "none"
    OBSERVATION_MASK = "observation_mask"
    EXISTING_IMPUTATION_PROVENANCE = "existing_imputation_provenance"
    UNSUPPORTED = "unsupported"


class CorrectedMissingCellAction(PolicyEnum):
    """Required handling for cells that were missing before temporary imputation."""

    RESTORE_MISSING = "restore_missing"
    MASK_MISSING = "mask_missing"
    FLAG_MISSING = "flag_missing"
    WITHHOLD_MISSING = "withhold_missing"
    UNSUPPORTED = "unsupported"


class RowSampleEligibilityImpact(PolicyEnum):
    """Correction eligibility impact of original missingness."""

    NO_CHANGE = "no_change"
    EXCLUDE_ROWS_WITH_ORIGINALLY_MISSING_VALUES = (
        "exclude_rows_with_originally_missing_values"
    )
    EXCLUDE_SAMPLES_WITH_ORIGINALLY_MISSING_VALUES = (
        "exclude_samples_with_originally_missing_values"
    )
    EXCLUDE_ROWS_WITH_INSUFFICIENT_OBSERVED_VALUES = (
        "exclude_rows_with_insufficient_observed_values"
    )
    REQUIRE_COMPLETE_CASES = "require_complete_cases"
    UNSUPPORTED = "unsupported"


TEMPORARY_IMPUTATION_METHODS = frozenset(TemporaryImputationMethod)
ORIGINALLY_MISSING_CELL_TRACKING_POLICIES = frozenset(OriginallyMissingCellTracking)
CORRECTED_MISSING_CELL_ACTIONS = frozenset(CorrectedMissingCellAction)
ROW_SAMPLE_ELIGIBILITY_IMPACTS = frozenset(RowSampleEligibilityImpact)

_REQUIRED_METHOD_PARAMETERS: Mapping[TemporaryImputationMethod, tuple[str, ...]] = {
    TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY: ("min_observed_values",),
    TemporaryImputationMethod.MINPROB_TEMPORARY: ("q", "width"),
    TemporaryImputationMethod.KNN_TEMPORARY: (
        "k",
        "distance",
        "max_missing_fraction_per_row",
    ),
}
_RANDOM_TEMPORARY_IMPUTATION_METHODS = frozenset(
    {TemporaryImputationMethod.MINPROB_TEMPORARY}
)
_UNSAFE_MISSINGNESS_ELIGIBILITY_IMPACTS = frozenset(
    {
        RowSampleEligibilityImpact.EXCLUDE_SAMPLES_WITH_ORIGINALLY_MISSING_VALUES,
        RowSampleEligibilityImpact.REQUIRE_COMPLETE_CASES,
        RowSampleEligibilityImpact.UNSUPPORTED,
    }
)


@dataclass(frozen=True, slots=True)
class ObservationMask:
    """Coordinate mask for cells that were missing in the original matrix."""

    feature_ids: tuple[str, ...]
    sample_ids: tuple[str, ...]
    originally_missing_cells: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        feature_ids = _require_non_empty_string_tuple(
            self.feature_ids,
            field_name="observation mask.feature_ids",
        )
        sample_ids = _require_non_empty_string_tuple(
            self.sample_ids,
            field_name="observation mask.sample_ids",
        )
        originally_missing_cells = _normalise_cell_coordinates(
            self.originally_missing_cells,
            feature_ids=feature_ids,
            sample_ids=sample_ids,
        )
        object.__setattr__(self, "feature_ids", feature_ids)
        object.__setattr__(self, "sample_ids", sample_ids)
        object.__setattr__(
            self,
            "originally_missing_cells",
            originally_missing_cells,
        )

    def is_originally_missing(self, feature_id: str, sample_id: str) -> bool:
        """Return whether one matrix coordinate was originally missing."""

        _require_known_coordinate(
            feature_id,
            sample_id,
            feature_ids=self.feature_ids,
            sample_ids=self.sample_ids,
        )
        return (feature_id, sample_id) in frozenset(self.originally_missing_cells)

    def is_originally_observed(self, feature_id: str, sample_id: str) -> bool:
        """Return whether one matrix coordinate was originally observed."""

        _require_known_coordinate(
            feature_id,
            sample_id,
            feature_ids=self.feature_ids,
            sample_ids=self.sample_ids,
        )
        return not self.is_originally_missing(feature_id, sample_id)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible payload for future provenance recording."""

        return {
            "feature_ids": list(self.feature_ids),
            "sample_ids": list(self.sample_ids),
            "originally_missing_cells": [
                {"feature_id": feature_id, "sample_id": sample_id}
                for feature_id, sample_id in self.originally_missing_cells
            ],
            "originally_missing_cell_count": len(self.originally_missing_cells),
        }


@dataclass(frozen=True, slots=True)
class TemporaryImputationPolicy:
    """Temporary-imputation intent for correction mechanics only."""

    allowed: bool = False
    method: TemporaryImputationMethod = TemporaryImputationMethod.NONE
    method_parameters: tuple[tuple[str, JsonScalar], ...] = ()
    random_seed: int | None = None
    supported: bool = True
    unsupported_reason: str | None = None
    imputed_values_are_observed_evidence: bool = False

    def __post_init__(self) -> None:
        method = TemporaryImputationMethod.parse(
            self.method,
            field_name="temporary imputation policy.method",
        )
        parameters = _normalise_method_parameters(self.method_parameters)
        random_seed = require_optional_int_at_least(
            self.random_seed,
            field_name="temporary imputation policy.random_seed",
            minimum=0,
            error_type=PhosPyInputError,
        )
        _validate_supported_state(
            supported=self.supported,
            unsupported_reason=self.unsupported_reason,
            field_prefix="temporary imputation policy",
        )
        if not isinstance(self.allowed, bool):
            raise PhosPyInputError("temporary imputation policy.allowed must be a bool")
        if not isinstance(self.imputed_values_are_observed_evidence, bool):
            raise PhosPyInputError(
                "temporary imputation policy.imputed_values_are_observed_evidence "
                "must be a bool"
            )
        if self.imputed_values_are_observed_evidence:
            raise PhosPyInputError(
                "temporary imputation policy must not treat temporary imputed "
                "values as observed biological evidence"
            )
        if method is TemporaryImputationMethod.UNSUPPORTED and self.supported:
            raise PhosPyInputError(
                "temporary imputation policy.supported must be False when "
                "method='unsupported'"
            )
        if method is TemporaryImputationMethod.NONE and self.allowed:
            raise PhosPyInputError(
                "temporary imputation policy.method must not be 'none' when "
                "allowed=True"
            )
        if method is TemporaryImputationMethod.NONE and (
            parameters or random_seed is not None
        ):
            raise PhosPyInputError(
                "temporary imputation policy.method_parameters and random_seed "
                "must be empty when method='none'"
            )
        if method is not TemporaryImputationMethod.NONE and not self.allowed:
            raise PhosPyInputError(
                "temporary imputation policy.allowed must be True when a "
                "temporary imputation method is requested"
            )
        if method is not TemporaryImputationMethod.UNSUPPORTED:
            _validate_method_parameters(method=method, parameters=parameters)
        if method in _RANDOM_TEMPORARY_IMPUTATION_METHODS and random_seed is None:
            raise PhosPyInputError(
                "temporary imputation policy.random_seed is required when "
                f"method='{method.value}' uses randomness"
            )
        object.__setattr__(self, "method", method)
        object.__setattr__(self, "method_parameters", parameters)
        object.__setattr__(self, "random_seed", random_seed)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible payload for future provenance recording."""

        return {
            "allowed": self.allowed,
            "method": self.method.value,
            "method_parameters": dict(self.method_parameters),
            "random_seed": self.random_seed,
            "supported": self.supported,
            "unsupported_reason": self.unsupported_reason,
            "imputed_values_are_observed_evidence": (
                self.imputed_values_are_observed_evidence
            ),
        }


@dataclass(frozen=True, slots=True)
class CorrectionMaskPolicy:
    """How future correction must handle originally missing cells."""

    corrected_missing_cell_action: CorrectedMissingCellAction = (
        CorrectedMissingCellAction.RESTORE_MISSING
    )
    preserve_observation_mask: bool = True
    supported: bool = True
    unsupported_reason: str | None = None

    def __post_init__(self) -> None:
        action = CorrectedMissingCellAction.parse(
            self.corrected_missing_cell_action,
            field_name="correction mask policy.corrected_missing_cell_action",
        )
        if not isinstance(self.preserve_observation_mask, bool):
            raise PhosPyInputError(
                "correction mask policy.preserve_observation_mask must be a bool"
            )
        if not self.preserve_observation_mask:
            raise PhosPyInputError(
                "correction mask policy must preserve the observation mask so "
                "originally missing cells remain distinguishable"
            )
        _validate_supported_state(
            supported=self.supported,
            unsupported_reason=self.unsupported_reason,
            field_prefix="correction mask policy",
        )
        if action is CorrectedMissingCellAction.UNSUPPORTED and self.supported:
            raise PhosPyInputError(
                "correction mask policy.supported must be False when "
                "corrected_missing_cell_action='unsupported'"
            )
        object.__setattr__(self, "corrected_missing_cell_action", action)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible payload for future provenance recording."""

        return {
            "corrected_missing_cell_action": (self.corrected_missing_cell_action.value),
            "preserve_observation_mask": self.preserve_observation_mask,
            "supported": self.supported,
            "unsupported_reason": self.unsupported_reason,
        }


@dataclass(frozen=True, slots=True)
class CorrectionMissingnessPolicy:
    """Missingness contract for future SPS/RUV-style correction requests."""

    temporary_imputation: TemporaryImputationPolicy = field(
        default_factory=TemporaryImputationPolicy
    )
    originally_missing_cells_tracked_by: OriginallyMissingCellTracking = (
        OriginallyMissingCellTracking.NONE
    )
    correction_mask_policy: CorrectionMaskPolicy = field(
        default_factory=CorrectionMaskPolicy
    )
    row_sample_eligibility_impact: RowSampleEligibilityImpact = (
        RowSampleEligibilityImpact.NO_CHANGE
    )
    observation_mask: ObservationMask | None = None
    supported: bool = True
    unsupported_reason: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.temporary_imputation, TemporaryImputationPolicy):
            raise PhosPyInputError(
                "correction missingness policy.temporary_imputation must be a "
                "TemporaryImputationPolicy"
            )
        tracking = OriginallyMissingCellTracking.parse(
            self.originally_missing_cells_tracked_by,
            field_name=(
                "correction missingness policy.originally_missing_cells_tracked_by"
            ),
        )
        if not isinstance(self.correction_mask_policy, CorrectionMaskPolicy):
            raise PhosPyInputError(
                "correction missingness policy.correction_mask_policy must be a "
                "CorrectionMaskPolicy"
            )
        eligibility = RowSampleEligibilityImpact.parse(
            self.row_sample_eligibility_impact,
            field_name=("correction missingness policy.row_sample_eligibility_impact"),
        )
        if self.observation_mask is not None and not isinstance(
            self.observation_mask,
            ObservationMask,
        ):
            raise PhosPyInputError(
                "correction missingness policy.observation_mask must be an "
                "ObservationMask when provided"
            )
        _validate_supported_state(
            supported=self.supported,
            unsupported_reason=self.unsupported_reason,
            field_prefix="correction missingness policy",
        )
        if self.temporary_imputation.allowed and tracking is (
            OriginallyMissingCellTracking.NONE
        ):
            raise PhosPyInputError(
                "correction missingness policy must track originally missing "
                "cells when temporary imputation is allowed"
            )
        if (
            tracking is OriginallyMissingCellTracking.UNSUPPORTED
            or eligibility is RowSampleEligibilityImpact.UNSUPPORTED
        ) and self.supported:
            raise PhosPyInputError(
                "correction missingness policy.supported must be False when "
                "tracking or row/sample eligibility impact is 'unsupported'"
            )
        if self.observation_mask is not None and tracking is not (
            OriginallyMissingCellTracking.OBSERVATION_MASK
        ):
            raise PhosPyInputError(
                "correction missingness policy.observation_mask requires "
                "originally_missing_cells_tracked_by='observation_mask'"
            )
        object.__setattr__(self, "originally_missing_cells_tracked_by", tracking)
        object.__setattr__(self, "row_sample_eligibility_impact", eligibility)

    def to_payload(self) -> dict[str, object]:
        """Return a JSON-compatible payload for future provenance recording."""

        return {
            "temporary_imputation": self.temporary_imputation.to_payload(),
            "originally_missing_cells_tracked_by": (
                self.originally_missing_cells_tracked_by.value
            ),
            "correction_mask_policy": self.correction_mask_policy.to_payload(),
            "row_sample_eligibility_impact": (self.row_sample_eligibility_impact.value),
            "observation_mask": (
                None
                if self.observation_mask is None
                else self.observation_mask.to_payload()
            ),
            "supported": self.supported,
            "unsupported_reason": self.unsupported_reason,
        }


class CorrectionMissingnessCompatibilityValidator:
    """Validate missing-value policy compatibility before correction execution."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        policy: CorrectionMissingnessPolicy | None = None,
        allow_complete_observed_data: bool = True,
        context: str = "correction missingness validation",
    ) -> None:
        if not isinstance(phospho, pd.DataFrame):
            raise PhosPyInputError(f"{context} requires a phospho DataFrame")
        feature_ids = _normalise_axis_labels(
            phospho.index,
            field_name=f"{context} phospho.index",
        )
        sample_ids = _normalise_axis_labels(
            phospho.columns,
            field_name=f"{context} phospho.columns",
        )
        missing_cells = _matrix_missing_cells(
            phospho,
            feature_ids=feature_ids,
            sample_ids=sample_ids,
        )
        has_missing = bool(missing_cells)

        if policy is None:
            if has_missing:
                raise PhosPyInputError(
                    f"{context} found missing values but no explicit supported "
                    "correction missingness policy was provided"
                )
            if not allow_complete_observed_data:
                raise PhosPyInputError(
                    f"{context} requires an explicit correction missingness policy"
                )
            return

        if not isinstance(policy, CorrectionMissingnessPolicy):
            raise PhosPyInputError(
                f"{context} policy must be a CorrectionMissingnessPolicy"
            )
        _require_supported_execution_policy(policy=policy, context=context)
        if policy.temporary_imputation.allowed and policy.observation_mask is None:
            raise PhosPyInputError(
                f"{context} temporary imputation requires an observation mask"
            )

        if not has_missing:
            if not allow_complete_observed_data:
                raise PhosPyInputError(
                    f"{context} complete observed data are not allowed by this "
                    "validation boundary"
                )
            _validate_complete_matrix_mask(policy=policy, context=context)
            return

        if policy.observation_mask is None:
            raise PhosPyInputError(
                f"{context} found missing values but policy does not provide an "
                "observation mask for originally missing cells"
            )
        if (
            policy.originally_missing_cells_tracked_by
            is not OriginallyMissingCellTracking.OBSERVATION_MASK
        ):
            raise PhosPyInputError(
                f"{context} requires originally_missing_cells_tracked_by="
                "'observation_mask' when missing values are present"
            )
        _validate_observation_mask_matches_matrix(
            mask=policy.observation_mask,
            feature_ids=feature_ids,
            sample_ids=sample_ids,
            missing_cells=missing_cells,
            context=context,
        )
        if (
            policy.row_sample_eligibility_impact
            in _UNSAFE_MISSINGNESS_ELIGIBILITY_IMPACTS
        ):
            raise PhosPyInputError(
                f"{context} row/sample eligibility impact "
                f"{policy.row_sample_eligibility_impact.value!r} makes correction "
                "unsafe for matrices with originally missing values"
            )


def _require_non_empty_string_tuple(
    value: object,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise PhosPyInputError(f"{field_name} must be a non-empty sequence")
    resolved = tuple(
        require_non_empty_string(
            item,
            field_name=f"{field_name}[]",
            error_type=PhosPyInputError,
        )
        for item in value
    )
    if not resolved:
        raise PhosPyInputError(f"{field_name} must be a non-empty sequence")
    if len(set(resolved)) != len(resolved):
        raise PhosPyInputError(f"{field_name} must not contain duplicates")
    return resolved


def _normalise_cell_coordinates(
    value: object,
    *,
    feature_ids: tuple[str, ...],
    sample_ids: tuple[str, ...],
) -> tuple[tuple[str, str], ...]:
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise PhosPyInputError(
            "observation mask.originally_missing_cells must be a sequence of "
            "(feature_id, sample_id) pairs"
        )
    coordinates: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for cell in value:
        if not isinstance(cell, Sequence) or isinstance(cell, str) or len(cell) != 2:
            raise PhosPyInputError(
                "observation mask.originally_missing_cells must contain only "
                "(feature_id, sample_id) pairs"
            )
        feature_id = require_non_empty_string(
            cell[0],
            field_name="observation mask.originally_missing_cells[].feature_id",
            error_type=PhosPyInputError,
        )
        sample_id = require_non_empty_string(
            cell[1],
            field_name="observation mask.originally_missing_cells[].sample_id",
            error_type=PhosPyInputError,
        )
        _require_known_coordinate(
            feature_id,
            sample_id,
            feature_ids=feature_ids,
            sample_ids=sample_ids,
        )
        coordinate = (feature_id, sample_id)
        if coordinate in seen:
            raise PhosPyInputError(
                "observation mask.originally_missing_cells must not contain "
                "duplicate coordinates"
            )
        seen.add(coordinate)
        coordinates.append(coordinate)
    return tuple(coordinates)


def _require_known_coordinate(
    feature_id: str,
    sample_id: str,
    *,
    feature_ids: tuple[str, ...],
    sample_ids: tuple[str, ...],
) -> None:
    if feature_id not in feature_ids:
        raise PhosPyInputError(
            "observation mask coordinate feature_id must be present in "
            f"feature_ids; got {feature_id!r}"
        )
    if sample_id not in sample_ids:
        raise PhosPyInputError(
            "observation mask coordinate sample_id must be present in "
            f"sample_ids; got {sample_id!r}"
        )


def _normalise_method_parameters(value: object) -> tuple[tuple[str, JsonScalar], ...]:
    if value is None:
        return ()
    if isinstance(value, Mapping):
        items = tuple(value.items())
    elif isinstance(value, Sequence) and not isinstance(value, str):
        items = tuple(value)
    else:
        raise PhosPyInputError(
            "temporary imputation policy.method_parameters must be a mapping or "
            "sequence of (name, value) pairs"
        )
    resolved: list[tuple[str, JsonScalar]] = []
    seen: set[str] = set()
    for item in items:
        if not isinstance(item, Sequence) or isinstance(item, str) or len(item) != 2:
            raise PhosPyInputError(
                "temporary imputation policy.method_parameters must contain only "
                "(name, value) pairs"
            )
        name = require_non_empty_string(
            item[0],
            field_name="temporary imputation policy.method_parameters[].name",
            error_type=PhosPyInputError,
        )
        if name in seen:
            raise PhosPyInputError(
                "temporary imputation policy.method_parameters must not contain "
                "duplicate names"
            )
        seen.add(name)
        parameter_value = _require_json_scalar(item[1])
        resolved.append((name, parameter_value))
    return tuple(resolved)


def _require_json_scalar(value: object) -> JsonScalar:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise PhosPyInputError(
                "temporary imputation policy.method_parameters values must be "
                "finite JSON scalar values"
            )
        return value
    raise PhosPyInputError(
        "temporary imputation policy.method_parameters values must be JSON "
        "scalar values"
    )


def _validate_method_parameters(
    *,
    method: TemporaryImputationMethod,
    parameters: tuple[tuple[str, JsonScalar], ...],
) -> None:
    required = _REQUIRED_METHOD_PARAMETERS.get(method, ())
    if not required:
        if parameters:
            raise PhosPyInputError(
                "temporary imputation policy.method_parameters must be empty "
                f"when method='{method.value}'"
            )
        return

    values = dict(parameters)
    missing = tuple(name for name in required if name not in values)
    if missing:
        raise PhosPyInputError(
            "temporary imputation policy.method_parameters missing required "
            f"parameter(s) for method='{method.value}': {', '.join(missing)}"
        )
    unexpected = tuple(name for name in values if name not in required)
    if unexpected:
        raise PhosPyInputError(
            "temporary imputation policy.method_parameters contains unsupported "
            f"parameter(s) for method='{method.value}': {', '.join(unexpected)}"
        )

    if method is TemporaryImputationMethod.ROW_MEDIAN_TEMPORARY:
        _require_int_parameter_at_least(
            values["min_observed_values"],
            name="min_observed_values",
            minimum=1,
            method=method,
        )
    elif method is TemporaryImputationMethod.MINPROB_TEMPORARY:
        _require_float_parameter_between(
            values["q"],
            name="q",
            lower=0.0,
            upper=0.5,
            include_upper=False,
            method=method,
        )
        _require_float_parameter_between(
            values["width"],
            name="width",
            lower=0.0,
            upper=1.0,
            include_upper=True,
            method=method,
        )
    elif method is TemporaryImputationMethod.KNN_TEMPORARY:
        _require_int_parameter_at_least(
            values["k"],
            name="k",
            minimum=1,
            method=method,
        )
        if values["distance"] != "nan_euclidean":
            raise PhosPyInputError(
                "temporary imputation policy.method_parameters['distance'] must "
                f"be 'nan_euclidean' when method='{method.value}'"
            )
        _require_float_parameter_between(
            values["max_missing_fraction_per_row"],
            name="max_missing_fraction_per_row",
            lower=0.0,
            upper=1.0,
            include_upper=True,
            method=method,
        )


def _require_int_parameter_at_least(
    value: JsonScalar,
    *,
    name: str,
    minimum: int,
    method: TemporaryImputationMethod,
) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise PhosPyInputError(
            f"temporary imputation policy.method_parameters['{name}'] must be "
            f"an int when method='{method.value}'"
        )
    if value < minimum:
        raise PhosPyInputError(
            f"temporary imputation policy.method_parameters['{name}'] must be "
            f"greater than or equal to {minimum} when method='{method.value}'"
        )


def _require_float_parameter_between(
    value: JsonScalar,
    *,
    name: str,
    lower: float,
    upper: float,
    include_upper: bool,
    method: TemporaryImputationMethod,
) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PhosPyInputError(
            f"temporary imputation policy.method_parameters['{name}'] must be "
            f"a float when method='{method.value}'"
        )
    numeric = float(value)
    upper_ok = numeric <= upper if include_upper else numeric < upper
    if not math.isfinite(numeric) or not (lower < numeric and upper_ok):
        operator = "<=" if include_upper else "<"
        raise PhosPyInputError(
            f"temporary imputation policy.method_parameters['{name}'] must "
            f"satisfy {lower:g} < value {operator} {upper:g} when "
            f"method='{method.value}'"
        )


def _require_supported_execution_policy(
    *,
    policy: CorrectionMissingnessPolicy,
    context: str,
) -> None:
    if not policy.supported:
        raise PhosPyInputError(
            f"{context} policy is unsupported: {policy.unsupported_reason}"
        )
    if not policy.temporary_imputation.supported:
        raise PhosPyInputError(
            f"{context} temporary imputation policy is unsupported: "
            f"{policy.temporary_imputation.unsupported_reason}"
        )
    if not policy.correction_mask_policy.supported:
        raise PhosPyInputError(
            f"{context} correction mask policy is unsupported: "
            f"{policy.correction_mask_policy.unsupported_reason}"
        )
    if (
        policy.temporary_imputation.method is TemporaryImputationMethod.UNSUPPORTED
        or policy.originally_missing_cells_tracked_by
        is OriginallyMissingCellTracking.UNSUPPORTED
        or policy.correction_mask_policy.corrected_missing_cell_action
        is CorrectedMissingCellAction.UNSUPPORTED
        or policy.row_sample_eligibility_impact
        is RowSampleEligibilityImpact.UNSUPPORTED
    ):
        raise PhosPyInputError(
            f"{context} contains unsupported missingness or imputation policies"
        )
    if not policy.correction_mask_policy.preserve_observation_mask:
        raise PhosPyInputError(
            f"{context} correction mask policy must preserve observation masks"
        )


def _validate_complete_matrix_mask(
    *,
    policy: CorrectionMissingnessPolicy,
    context: str,
) -> None:
    if (
        policy.observation_mask is not None
        and policy.observation_mask.originally_missing_cells
    ):
        raise PhosPyInputError(
            f"{context} observation mask declares originally missing cells for a "
            "complete observed matrix"
        )


def _validate_observation_mask_matches_matrix(
    *,
    mask: ObservationMask,
    feature_ids: tuple[str, ...],
    sample_ids: tuple[str, ...],
    missing_cells: tuple[tuple[str, str], ...],
    context: str,
) -> None:
    if mask.feature_ids != feature_ids:
        raise PhosPyInputError(
            f"{context} observation mask feature_ids must match phospho.index"
        )
    if mask.sample_ids != sample_ids:
        raise PhosPyInputError(
            f"{context} observation mask sample_ids must match phospho.columns"
        )
    if frozenset(mask.originally_missing_cells) != frozenset(missing_cells):
        raise PhosPyInputError(
            f"{context} observation mask originally_missing_cells must match the "
            "matrix missing-value coordinates"
        )


def _matrix_missing_cells(
    phospho: pd.DataFrame,
    *,
    feature_ids: tuple[str, ...],
    sample_ids: tuple[str, ...],
) -> tuple[tuple[str, str], ...]:
    missing = phospho.isna().to_numpy()
    cells: list[tuple[str, str]] = []
    row_count, column_count = missing.shape
    for row_index in range(row_count):
        for column_index in range(column_count):
            if missing[row_index, column_index]:
                cells.append((feature_ids[row_index], sample_ids[column_index]))
    return tuple(cells)


def _normalise_axis_labels(index: pd.Index, *, field_name: str) -> tuple[str, ...]:
    labels: list[str] = []
    for position, value in enumerate(index.tolist()):
        if bool(pd.Series((value,), dtype="object").isna().iat[0]):
            raise PhosPyInputError(
                f"{field_name} must not contain missing labels; position {position}"
            )
        label = str(value).strip()
        if not label:
            raise PhosPyInputError(
                f"{field_name} must not contain blank labels; position {position}"
            )
        labels.append(label)
    if len(set(labels)) != len(labels):
        duplicates = tuple(label for label in labels if labels.count(label) > 1)
        raise PhosPyInputError(
            f"{field_name} must not contain duplicate labels: "
            f"{', '.join(dict.fromkeys(duplicates))}"
        )
    return tuple(labels)


def _validate_supported_state(
    *,
    supported: object,
    unsupported_reason: object | None,
    field_prefix: str,
) -> None:
    if not isinstance(supported, bool):
        raise PhosPyInputError(f"{field_prefix}.supported must be a bool")
    if supported:
        if unsupported_reason is not None:
            raise PhosPyInputError(
                f"{field_prefix}.unsupported_reason must be None when supported=True"
            )
        return
    require_non_empty_string(
        unsupported_reason,
        field_name=f"{field_prefix}.unsupported_reason",
        error_type=PhosPyInputError,
        when_provided=True,
    )


__all__ = [
    "CORRECTED_MISSING_CELL_ACTIONS",
    "ORIGINALLY_MISSING_CELL_TRACKING_POLICIES",
    "ROW_SAMPLE_ELIGIBILITY_IMPACTS",
    "TEMPORARY_IMPUTATION_METHODS",
    "CorrectedMissingCellAction",
    "CorrectionMaskPolicy",
    "CorrectionMissingnessCompatibilityValidator",
    "CorrectionMissingnessPolicy",
    "ObservationMask",
    "OriginallyMissingCellTracking",
    "RowSampleEligibilityImpact",
    "TemporaryImputationMethod",
    "TemporaryImputationPolicy",
]
