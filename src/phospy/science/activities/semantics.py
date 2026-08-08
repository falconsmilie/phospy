"""Typed kinase-activity input semantics.

This module defines the interpretation boundary for matrices consumed by
activity methods. It deliberately stores profile-axis and quantitative meaning
as structured objects instead of inferring them from matrix column labels.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from enum import Enum
from typing import cast

import pandas as pd

from phospy._deprecations import warn_deprecated
from phospy.errors.validation import PhosPyValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.frames.comparison import dataframe_equals
from phospy.frames.ownership import export_dataframe
from phospy.science.tables.activity import ActivityMatrix


class ActivityProfileAxis(str, Enum):
    """Axis represented by activity input matrix columns."""

    SAMPLE = "sample"
    CONDITION_SUMMARY = "condition_summary"
    CONTRAST = "contrast"
    EFFECT = "effect"


class ActivityQuantitativeSemantics(str, Enum):
    """Quantitative interpretation of values supplied to an activity method."""

    SAMPLE_LEVEL_ABUNDANCE = "sample_level_abundance"
    CONDITION_SUMMARY_ABUNDANCE = "condition_summary_abundance"
    CONTRAST_LOG_FOLD_CHANGE = "contrast_log_fold_change"
    STANDARDISED_EFFECT = "standardised_effect"


_ACTIVITY_INPUT_SEMANTIC_AXIS_BY_QUANTITY = {
    ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE: ActivityProfileAxis.SAMPLE,
    ActivityQuantitativeSemantics.CONDITION_SUMMARY_ABUNDANCE: (
        ActivityProfileAxis.CONDITION_SUMMARY
    ),
    ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE: (
        ActivityProfileAxis.CONTRAST
    ),
    ActivityQuantitativeSemantics.STANDARDISED_EFFECT: ActivityProfileAxis.EFFECT,
}
_CONTRAST_OR_EFFECT_QUANTITIES = frozenset(
    {
        ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE,
        ActivityQuantitativeSemantics.STANDARDISED_EFFECT,
    }
)


@dataclass(frozen=True, slots=True)
class ActivityAggregationRecord:
    """One condition-summary profile and the source profiles aggregated into it."""

    profile_id: str
    source_profile_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "profile_id",
            _require_non_empty_text(
                self.profile_id,
                field_name="activity_aggregation_record.profile_id",
            ),
        )
        object.__setattr__(
            self,
            "source_profile_ids",
            _normalize_non_empty_text_tuple(
                self.source_profile_ids,
                field_name="activity_aggregation_record.source_profile_ids",
            ),
        )

    def to_payload(self) -> dict[str, object]:
        return {
            "profile_id": self.profile_id,
            "source_profile_ids": list(self.source_profile_ids),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityAggregationRecord:
        return cls(
            profile_id=_payload_text(
                payload,
                "profile_id",
                field_name="activity_aggregation_record.profile_id",
            ),
            source_profile_ids=tuple(
                _payload_text_sequence(
                    payload,
                    "source_profile_ids",
                    field_name="activity_aggregation_record.source_profile_ids",
                )
            ),
        )


@dataclass(frozen=True, slots=True)
class ActivityAggregationMetadata:
    """Typed aggregation metadata required for condition-summary abundance."""

    aggregation_method: str
    records: tuple[ActivityAggregationRecord, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "aggregation_method",
            _require_non_empty_text(
                self.aggregation_method,
                field_name="activity_aggregation_metadata.aggregation_method",
            ),
        )
        records = tuple(self.records)
        if not records:
            raise WorkflowBoundaryError(
                "activity_aggregation_metadata.records must contain at least one "
                "ActivityAggregationRecord"
            )
        for record in records:
            if not isinstance(record, ActivityAggregationRecord):
                raise WorkflowBoundaryError(
                    "activity_aggregation_metadata.records must contain "
                    "ActivityAggregationRecord values"
                )
        profile_ids = [record.profile_id for record in records]
        if len(profile_ids) != len(set(profile_ids)):
            raise WorkflowBoundaryError(
                "activity_aggregation_metadata.records profile_id values must be unique"
            )
        object.__setattr__(self, "records", records)

    @property
    def profile_ids(self) -> tuple[str, ...]:
        return tuple(record.profile_id for record in self.records)

    def to_payload(self) -> dict[str, object]:
        return {
            "aggregation_method": self.aggregation_method,
            "records": [record.to_payload() for record in self.records],
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityAggregationMetadata:
        raw_records = payload.get("records")
        if not isinstance(raw_records, Sequence) or isinstance(
            raw_records,
            (str, bytes, bytearray),
        ):
            raise ValueError("activity_aggregation_metadata.records must be a list")
        records: list[ActivityAggregationRecord] = []
        for position, raw_record in enumerate(raw_records):
            if not isinstance(raw_record, Mapping):
                raise ValueError(
                    "activity_aggregation_metadata.records "
                    f"[{position}] must be an object"
                )
            records.append(ActivityAggregationRecord.from_payload(raw_record))
        return cls(
            aggregation_method=_payload_text(
                payload,
                "aggregation_method",
                field_name="activity_aggregation_metadata.aggregation_method",
            ),
            records=tuple(records),
        )


@dataclass(frozen=True, slots=True)
class ActivityInputSemantics:
    """Typed activity input axis and quantitative semantics."""

    profile_axis: ActivityProfileAxis | str
    quantitative_semantics: ActivityQuantitativeSemantics | str

    def __post_init__(self) -> None:
        axis = _normalize_profile_axis(self.profile_axis)
        quantity = _normalize_quantitative_semantics(self.quantitative_semantics)
        expected_axis = _ACTIVITY_INPUT_SEMANTIC_AXIS_BY_QUANTITY[quantity]
        if axis is not expected_axis:
            raise WorkflowBoundaryError(
                "activity input semantics are inconsistent: quantitative semantics "
                f"{quantity.value!r} requires profile_axis={expected_axis.value!r}, "
                f"got {axis.value!r}"
            )
        object.__setattr__(self, "profile_axis", axis)
        object.__setattr__(self, "quantitative_semantics", quantity)

    @property
    def is_abundance(self) -> bool:
        return self.quantitative_semantics in {
            ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE,
            ActivityQuantitativeSemantics.CONDITION_SUMMARY_ABUNDANCE,
        }

    @property
    def requires_contrast_or_effect_input(self) -> bool:
        return self.quantitative_semantics in _CONTRAST_OR_EFFECT_QUANTITIES

    @property
    def has_real_condition_contract(self) -> bool:
        return (
            self.profile_axis is ActivityProfileAxis.CONDITION_SUMMARY
            and self.quantitative_semantics
            is ActivityQuantitativeSemantics.CONDITION_SUMMARY_ABUNDANCE
        )

    def to_payload(self) -> dict[str, str]:
        axis = cast(ActivityProfileAxis, self.profile_axis)
        quantity = cast(ActivityQuantitativeSemantics, self.quantitative_semantics)
        return {
            "profile_axis": axis.value,
            "quantitative_semantics": quantity.value,
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityInputSemantics:
        return cls(
            profile_axis=_payload_text(
                payload,
                "profile_axis",
                field_name="activity_input_semantics.profile_axis",
            ),
            quantitative_semantics=_payload_text(
                payload,
                "quantitative_semantics",
                field_name="activity_input_semantics.quantitative_semantics",
            ),
        )


@dataclass(frozen=True, slots=True)
class ActivityProfileMetadata:
    """Typed labels for activity input/output profiles."""

    profile_ids: tuple[str, ...]
    axis: ActivityProfileAxis | str
    sample_ids: tuple[str, ...] = ()
    condition_ids: tuple[str, ...] = ()
    contrast_ids: tuple[str, ...] = ()
    aggregation_metadata: ActivityAggregationMetadata | None = None

    def __post_init__(self) -> None:
        axis = _normalize_profile_axis(self.axis)
        profile_ids = _normalize_optional_text_tuple(
            self.profile_ids,
            field_name="activity_profile_metadata.profile_ids",
        )
        if len(profile_ids) != len(set(profile_ids)):
            raise WorkflowBoundaryError(
                "activity_profile_metadata.profile_ids must be unique"
            )
        sample_ids = _normalize_optional_text_tuple(
            self.sample_ids,
            field_name="activity_profile_metadata.sample_ids",
        )
        condition_ids = _normalize_optional_text_tuple(
            self.condition_ids,
            field_name="activity_profile_metadata.condition_ids",
        )
        contrast_ids = _normalize_optional_text_tuple(
            self.contrast_ids,
            field_name="activity_profile_metadata.contrast_ids",
        )
        aggregation_metadata = self.aggregation_metadata
        if aggregation_metadata is not None and not isinstance(
            aggregation_metadata,
            ActivityAggregationMetadata,
        ):
            raise WorkflowBoundaryError(
                "activity_profile_metadata.aggregation_metadata must be "
                "ActivityAggregationMetadata or None"
            )
        if axis is ActivityProfileAxis.SAMPLE:
            sample_ids = profile_ids if not sample_ids else sample_ids
            _require_matching_labels(
                observed=sample_ids,
                expected=profile_ids,
                field_name="activity_profile_metadata.sample_ids",
            )
        elif axis is ActivityProfileAxis.CONDITION_SUMMARY:
            condition_ids = profile_ids if not condition_ids else condition_ids
            _require_matching_labels(
                observed=condition_ids,
                expected=profile_ids,
                field_name="activity_profile_metadata.condition_ids",
            )
            if aggregation_metadata is None:
                raise WorkflowBoundaryError(
                    "condition-summary activity input requires explicit "
                    "ActivityAggregationMetadata"
                )
            _require_matching_labels(
                observed=aggregation_metadata.profile_ids,
                expected=profile_ids,
                field_name="activity_profile_metadata.aggregation_metadata.profile_ids",
            )
        elif axis is ActivityProfileAxis.CONTRAST:
            contrast_ids = profile_ids if not contrast_ids else contrast_ids
            _require_matching_labels(
                observed=contrast_ids,
                expected=profile_ids,
                field_name="activity_profile_metadata.contrast_ids",
            )
        elif axis is ActivityProfileAxis.EFFECT:
            pass
        object.__setattr__(self, "axis", axis)
        object.__setattr__(self, "profile_ids", profile_ids)
        object.__setattr__(self, "sample_ids", sample_ids)
        object.__setattr__(self, "condition_ids", condition_ids)
        object.__setattr__(self, "contrast_ids", contrast_ids)
        object.__setattr__(self, "aggregation_metadata", aggregation_metadata)

    def to_payload(self) -> dict[str, object]:
        axis = cast(ActivityProfileAxis, self.axis)
        return {
            "axis": axis.value,
            "profile_ids": list(self.profile_ids),
            "sample_ids": list(self.sample_ids),
            "condition_ids": list(self.condition_ids),
            "contrast_ids": list(self.contrast_ids),
            "aggregation_metadata": (
                None
                if self.aggregation_metadata is None
                else self.aggregation_metadata.to_payload()
            ),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
    ) -> ActivityProfileMetadata:
        raw_aggregation = payload.get("aggregation_metadata")
        aggregation = None
        if raw_aggregation is not None:
            if not isinstance(raw_aggregation, Mapping):
                raise ValueError(
                    "activity_profile_metadata.aggregation_metadata must be an object"
                )
            aggregation = ActivityAggregationMetadata.from_payload(raw_aggregation)
        return cls(
            axis=_payload_text(
                payload,
                "axis",
                field_name="activity_profile_metadata.axis",
            ),
            profile_ids=tuple(
                _payload_text_sequence(
                    payload,
                    "profile_ids",
                    field_name="activity_profile_metadata.profile_ids",
                )
            ),
            sample_ids=tuple(
                _payload_text_sequence(
                    payload,
                    "sample_ids",
                    field_name="activity_profile_metadata.sample_ids",
                    allow_missing=True,
                )
            ),
            condition_ids=tuple(
                _payload_text_sequence(
                    payload,
                    "condition_ids",
                    field_name="activity_profile_metadata.condition_ids",
                    allow_missing=True,
                )
            ),
            contrast_ids=tuple(
                _payload_text_sequence(
                    payload,
                    "contrast_ids",
                    field_name="activity_profile_metadata.contrast_ids",
                    allow_missing=True,
                )
            ),
            aggregation_metadata=aggregation,
        )


@dataclass(frozen=True, slots=True, eq=False)
class ActivityInputMatrix:
    """Activity input matrix paired with explicit semantics.

    Python equality and hashing are identity-based. Use
    :meth:`scientifically_equals` for explicit matrix-content comparison.
    """

    __hash__ = object.__hash__

    frame: pd.DataFrame
    semantics: ActivityInputSemantics
    profile_metadata: ActivityProfileMetadata
    field_name: str = field(default="activity_input.matrix", repr=False, compare=False)
    _assume_owned: bool = field(default=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.semantics, ActivityInputSemantics):
            raise WorkflowBoundaryError(
                "activity_input.semantics must be ActivityInputSemantics"
            )
        if not isinstance(self.profile_metadata, ActivityProfileMetadata):
            raise WorkflowBoundaryError(
                "activity_input.profile_metadata must be ActivityProfileMetadata"
            )
        frame = ActivityMatrix(
            frame=self.frame,
            field_name=self.field_name,
            _assume_owned=bool(self._assume_owned),
        ).frame
        profile_ids = tuple(str(column) for column in frame.columns)
        _require_matching_labels(
            observed=self.profile_metadata.profile_ids,
            expected=profile_ids,
            field_name="activity_input.profile_metadata.profile_ids",
        )
        if self.profile_metadata.axis is not self.semantics.profile_axis:
            raise WorkflowBoundaryError(
                "activity_input.profile_metadata.axis must match "
                "activity_input.semantics.profile_axis"
            )
        if (
            self.semantics.profile_axis is ActivityProfileAxis.CONDITION_SUMMARY
            and self.profile_metadata.aggregation_metadata is None
        ):
            raise WorkflowBoundaryError(
                "condition-summary activity input requires explicit aggregation "
                "metadata"
            )
        object.__setattr__(self, "frame", frame)

    @property
    def matrix(self) -> pd.DataFrame:
        """Return a defensive matrix snapshot."""

        return export_dataframe(self.frame)

    def scientifically_equals(self, other: object) -> bool:
        """Return ``True`` when another activity input has the same content."""

        if not isinstance(other, ActivityInputMatrix):
            return False
        return (
            dataframe_equals(self.frame, other.frame)
            and self.semantics == other.semantics
            and self.profile_metadata == other.profile_metadata
        )

    @classmethod
    def sample_level_abundance(
        cls,
        frame: pd.DataFrame,
        *,
        field_name: str = "activity_input.sample_level_abundance",
        _assume_owned: bool = False,
    ) -> ActivityInputMatrix:
        profile_ids = tuple(str(column) for column in frame.columns)
        return cls(
            frame=frame,
            semantics=ActivityInputSemantics(
                profile_axis=ActivityProfileAxis.SAMPLE,
                quantitative_semantics=(
                    ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE
                ),
            ),
            profile_metadata=ActivityProfileMetadata(
                axis=ActivityProfileAxis.SAMPLE,
                profile_ids=profile_ids,
                sample_ids=profile_ids,
            ),
            field_name=field_name,
            _assume_owned=_assume_owned,
        )

    @classmethod
    def condition_summary_abundance(
        cls,
        frame: pd.DataFrame,
        *,
        aggregation_metadata: ActivityAggregationMetadata,
        field_name: str = "activity_input.condition_summary_abundance",
        _assume_owned: bool = False,
    ) -> ActivityInputMatrix:
        profile_ids = tuple(str(column) for column in frame.columns)
        return cls(
            frame=frame,
            semantics=ActivityInputSemantics(
                profile_axis=ActivityProfileAxis.CONDITION_SUMMARY,
                quantitative_semantics=(
                    ActivityQuantitativeSemantics.CONDITION_SUMMARY_ABUNDANCE
                ),
            ),
            profile_metadata=ActivityProfileMetadata(
                axis=ActivityProfileAxis.CONDITION_SUMMARY,
                profile_ids=profile_ids,
                condition_ids=profile_ids,
                aggregation_metadata=aggregation_metadata,
            ),
            field_name=field_name,
            _assume_owned=_assume_owned,
        )

    @classmethod
    def contrast_log_fold_change(
        cls,
        frame: pd.DataFrame,
        *,
        field_name: str = "activity_input.contrast_log_fold_change",
        _assume_owned: bool = False,
    ) -> ActivityInputMatrix:
        profile_ids = tuple(str(column) for column in frame.columns)
        return cls(
            frame=frame,
            semantics=ActivityInputSemantics(
                profile_axis=ActivityProfileAxis.CONTRAST,
                quantitative_semantics=(
                    ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE
                ),
            ),
            profile_metadata=ActivityProfileMetadata(
                axis=ActivityProfileAxis.CONTRAST,
                profile_ids=profile_ids,
                contrast_ids=profile_ids,
            ),
            field_name=field_name,
            _assume_owned=_assume_owned,
        )

    @classmethod
    def standardised_effect(
        cls,
        frame: pd.DataFrame,
        *,
        field_name: str = "activity_input.standardised_effect",
        _assume_owned: bool = False,
    ) -> ActivityInputMatrix:
        profile_ids = tuple(str(column) for column in frame.columns)
        return cls(
            frame=frame,
            semantics=ActivityInputSemantics(
                profile_axis=ActivityProfileAxis.EFFECT,
                quantitative_semantics=(
                    ActivityQuantitativeSemantics.STANDARDISED_EFFECT
                ),
            ),
            profile_metadata=ActivityProfileMetadata(
                axis=ActivityProfileAxis.EFFECT,
                profile_ids=profile_ids,
            ),
            field_name=field_name,
            _assume_owned=_assume_owned,
        )

    @classmethod
    def from_payload(
        cls,
        *,
        frame: pd.DataFrame,
        semantics_payload: Mapping[str, object],
        profile_metadata_payload: Mapping[str, object],
        field_name: str = "activity_input.matrix",
    ) -> ActivityInputMatrix:
        return cls(
            frame=frame,
            semantics=ActivityInputSemantics.from_payload(semantics_payload),
            profile_metadata=ActivityProfileMetadata.from_payload(
                profile_metadata_payload
            ),
            field_name=field_name,
        )


def normalize_activity_input_matrix(
    value: ActivityInputMatrix | pd.DataFrame,
    *,
    field_name: str,
    legacy_dataframe_semantics: ActivityInputSemantics | None = None,
    legacy_dataframe_warning: str | None = None,
) -> ActivityInputMatrix:
    """Return a typed activity input matrix from current or legacy inputs."""

    if isinstance(value, ActivityInputMatrix):
        return value
    if not isinstance(value, pd.DataFrame):
        raise WorkflowBoundaryError(f"{field_name} must be ActivityInputMatrix")
    if legacy_dataframe_semantics is None:
        raise WorkflowBoundaryError(
            f"{field_name} must be ActivityInputMatrix with explicit semantics"
        )
    if legacy_dataframe_warning is not None:
        warn_deprecated(
            "activities.ssgsea.effect_matrix_dataframe",
            stacklevel=3,
        )
    quantity = cast(
        ActivityQuantitativeSemantics,
        legacy_dataframe_semantics.quantitative_semantics,
    )
    if quantity is ActivityQuantitativeSemantics.SAMPLE_LEVEL_ABUNDANCE:
        return ActivityInputMatrix.sample_level_abundance(
            value,
            field_name=field_name,
        )
    if quantity is ActivityQuantitativeSemantics.CONTRAST_LOG_FOLD_CHANGE:
        return ActivityInputMatrix.contrast_log_fold_change(
            value,
            field_name=field_name,
        )
    if quantity is ActivityQuantitativeSemantics.STANDARDISED_EFFECT:
        return ActivityInputMatrix.standardised_effect(
            value,
            field_name=field_name,
        )
    raise WorkflowBoundaryError(
        f"{field_name} cannot infer condition-summary aggregation metadata from a "
        "legacy DataFrame; provide ActivityInputMatrix.condition_summary_abundance"
    )


def require_abundance_activity_input(
    activity_input: ActivityInputMatrix,
    *,
    field_name: str,
) -> ActivityInputMatrix:
    if activity_input.semantics.is_abundance:
        return activity_input
    quantity = cast(
        ActivityQuantitativeSemantics,
        activity_input.semantics.quantitative_semantics,
    )
    raise WorkflowBoundaryError(
        f"{field_name} requires sample-level or condition-summary abundance input; "
        f"got {quantity.value!r}"
    )


def require_contrast_or_effect_activity_input(
    activity_input: ActivityInputMatrix,
    *,
    field_name: str,
) -> ActivityInputMatrix:
    if activity_input.semantics.requires_contrast_or_effect_input:
        return activity_input
    quantity = cast(
        ActivityQuantitativeSemantics,
        activity_input.semantics.quantitative_semantics,
    )
    raise WorkflowBoundaryError(
        f"{field_name} requires explicit contrast/effect input; got {quantity.value!r}"
    )


def _normalize_profile_axis(value: ActivityProfileAxis | str) -> ActivityProfileAxis:
    if isinstance(value, ActivityProfileAxis):
        return value
    try:
        return ActivityProfileAxis(str(value))
    except ValueError as exc:
        allowed = ", ".join(axis.value for axis in ActivityProfileAxis)
        raise WorkflowBoundaryError(
            f"activity profile axis must be one of: {allowed}"
        ) from exc


def _normalize_quantitative_semantics(
    value: ActivityQuantitativeSemantics | str,
) -> ActivityQuantitativeSemantics:
    if isinstance(value, ActivityQuantitativeSemantics):
        return value
    normalized = str(value).strip()
    if normalized == "standardized_effect":
        normalized = ActivityQuantitativeSemantics.STANDARDISED_EFFECT.value
    try:
        return ActivityQuantitativeSemantics(normalized)
    except ValueError as exc:
        allowed = ", ".join(
            quantity.value for quantity in ActivityQuantitativeSemantics
        )
        raise WorkflowBoundaryError(
            f"activity quantitative semantics must be one of: {allowed}"
        ) from exc


def _normalize_non_empty_text_tuple(
    values: Sequence[object],
    *,
    field_name: str,
) -> tuple[str, ...]:
    normalized = _normalize_optional_text_tuple(values, field_name=field_name)
    if not normalized:
        raise WorkflowBoundaryError(f"{field_name} must contain at least one value")
    return normalized


def _normalize_optional_text_tuple(
    values: Sequence[object],
    *,
    field_name: str,
) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise WorkflowBoundaryError(f"{field_name} must be a sequence of strings")
    try:
        raw_values = tuple(values)
    except TypeError as exc:
        raise WorkflowBoundaryError(
            f"{field_name} must be a sequence of strings"
        ) from exc
    return tuple(
        _require_non_empty_text(value, field_name=f"{field_name}[{position}]")
        for position, value in enumerate(raw_values)
    )


def _require_matching_labels(
    *,
    observed: Sequence[object],
    expected: Sequence[object],
    field_name: str,
) -> None:
    observed_values = tuple(str(value) for value in observed)
    expected_values = tuple(str(value) for value in expected)
    if observed_values == expected_values:
        return
    raise WorkflowBoundaryError(
        f"{field_name} must exactly match activity matrix profile_ids; "
        f"expected={expected_values!r}, got={observed_values!r}"
    )


def _require_non_empty_text(value: object, *, field_name: str) -> str:
    if not isinstance(value, str):
        raise WorkflowBoundaryError(f"{field_name} must be a non-empty string")
    normalized = value.strip()
    if not normalized:
        raise WorkflowBoundaryError(f"{field_name} must be a non-empty string")
    return normalized


def _payload_text(
    payload: Mapping[str, object],
    key: str,
    *,
    field_name: str,
) -> str:
    value = payload.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")
    return value.strip()


def _payload_text_sequence(
    payload: Mapping[str, object],
    key: str,
    *,
    field_name: str,
    allow_missing: bool = False,
) -> tuple[str, ...]:
    value = payload.get(key)
    if value is None and allow_missing:
        return ()
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes, bytearray)):
        raise ValueError(f"{field_name} must be a list of strings")
    result: list[str] = []
    for position, item in enumerate(value):
        if not isinstance(item, str) or not item.strip():
            raise ValueError(f"{field_name}[{position}] must be a non-empty string")
        result.append(item.strip())
    return tuple(result)


def coerce_semantics_payload(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> ActivityInputSemantics:
    try:
        return ActivityInputSemantics.from_payload(value)
    except (ValueError, PhosPyValidationError, WorkflowBoundaryError) as exc:
        raise WorkflowBoundaryError(f"{field_name} is invalid: {exc}") from exc


def coerce_profile_metadata_payload(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> ActivityProfileMetadata:
    try:
        return ActivityProfileMetadata.from_payload(value)
    except (ValueError, PhosPyValidationError, WorkflowBoundaryError) as exc:
        raise WorkflowBoundaryError(f"{field_name} is invalid: {exc}") from exc


__all__ = [
    "ActivityAggregationMetadata",
    "ActivityAggregationRecord",
    "ActivityInputMatrix",
    "ActivityInputSemantics",
    "ActivityProfileAxis",
    "ActivityProfileMetadata",
    "ActivityQuantitativeSemantics",
    "coerce_profile_metadata_payload",
    "coerce_semantics_payload",
    "normalize_activity_input_matrix",
    "require_abundance_activity_input",
    "require_contrast_or_effect_activity_input",
]
