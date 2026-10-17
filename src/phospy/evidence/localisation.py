"""Localisation-confidence normalisation shared by phosphosite importers."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError

LOCALISATION_CONFIDENCE_SCALE_PROBABILITY = "probability"
LOCALISATION_CONFIDENCE_SCALE_PERCENT = "percent"
SUPPORTED_LOCALISATION_CONFIDENCE_SCALES: tuple[str, ...] = (
    LOCALISATION_CONFIDENCE_SCALE_PROBABILITY,
    LOCALISATION_CONFIDENCE_SCALE_PERCENT,
)
LOCALISATION_CONFIDENCE_OUTPUT_COLUMN = "localisation_confidence"
_EXAMPLE_LIMIT = 5


@dataclass(frozen=True, slots=True)
class LocalisationConfidenceNormalisationReport:
    """Diagnostics from normalising one localisation-confidence column."""

    source_column: str
    output_column: str
    scale: str
    row_count: int
    missing_count: int
    invalid_count: int
    invalid_examples: tuple[str, ...]

    def to_payload(self) -> dict[str, object]:
        return {
            "source_column": self.source_column,
            "output_column": self.output_column,
            "scale": self.scale,
            "row_count": int(self.row_count),
            "missing_count": int(self.missing_count),
            "invalid_count": int(self.invalid_count),
            "invalid_examples": list(self.invalid_examples),
        }


@dataclass(frozen=True, slots=True)
class LocalisationProbabilityAssessment:
    """Parsed localisation-probability assessment for one metadata column."""

    normalized: pd.Series
    missing_mask: pd.Series
    invalid_mask: pd.Series
    invalid_examples: tuple[str, ...]

    @property
    def missing_count(self) -> int:
        return int(self.missing_mask.sum())

    @property
    def invalid_count(self) -> int:
        return int(self.invalid_mask.sum())


def assess_localisation_probability_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[Exception],
    column_name: str = "localisation_probability",
) -> LocalisationProbabilityAssessment | None:
    """Parse optional localisation probability values with diagnostics."""

    if column_name not in site_metadata.columns:
        return None
    values = site_metadata[column_name]
    values_index = pd.Index(values.index)
    missing_mask = values.isna()
    blank_string_mask = values.map(
        lambda value: isinstance(value, str) and value.strip() == ""
    )
    bool_mask = values.map(lambda value: isinstance(value, (bool, np.bool_)))
    parse_exempt_mask = missing_mask | blank_string_mask | bool_mask
    numeric_values = pd.to_numeric(values.mask(parse_exempt_mask), errors="coerce")
    finite_mask = pd.Series(
        np.isfinite(numeric_values.to_numpy(dtype=float, copy=False, na_value=np.nan)),
        index=values_index,
    )
    valid_numeric_mask = (
        ~parse_exempt_mask
        & numeric_values.notna()
        & finite_mask
        & numeric_values.ge(0.0)
        & numeric_values.le(1.0)
    )
    missing_mask = pd.Series(
        (missing_mask | blank_string_mask).to_numpy(dtype=bool, copy=False),
        index=values_index,
        dtype="boolean",
    )
    invalid_mask = pd.Series(
        ((~missing_mask.astype(bool)) & (~valid_numeric_mask)).to_numpy(
            dtype=bool,
            copy=False,
        ),
        index=values_index,
        dtype="boolean",
    )
    normalized = pd.Series(pd.NA, index=values_index, dtype="Float64")
    if bool(valid_numeric_mask.any()):
        normalized.loc[valid_numeric_mask] = numeric_values.loc[
            valid_numeric_mask
        ].astype(float)
    invalid_examples: list[str] = []
    invalid_positions = np.flatnonzero(invalid_mask.to_numpy(dtype=bool, copy=False))
    for position in invalid_positions[:_EXAMPLE_LIMIT]:
        site_id = values.index[int(position)]
        raw_value = values.at[site_id]
        parsed = _parse_confidence_value(raw_value, scale="probability")
        invalid_examples.append(f"{site_id!r}:{raw_value!r}:{parsed}")

    if bool((missing_mask & invalid_mask).any()):
        raise error_type(
            f"{field_name}.{column_name} localisation parsing produced inconsistent "
            "missing/invalid masks"
        )
    return LocalisationProbabilityAssessment(
        normalized=normalized,
        missing_mask=missing_mask,
        invalid_mask=invalid_mask,
        invalid_examples=tuple(invalid_examples),
    )


def validate_localisation_probability_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[Exception],
    column_name: str = "localisation_probability",
) -> None:
    """Validate optional localisation probability values when present."""

    assessment = assess_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=column_name,
    )
    if assessment is None or assessment.invalid_count == 0:
        return
    raise error_type(
        f"{field_name}.{column_name} must contain values in [0.0, 1.0] or missing; "
        f"invalid_row_count={assessment.invalid_count}; "
        f"examples={_summarise_examples(list(assessment.invalid_examples), limit=3)}"
    )


def validate_localisation_confidence_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[Exception],
    column_name: str = "localisation_confidence",
) -> None:
    """Validate optional localisation-confidence values when present."""

    validate_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=column_name,
    )


def normalise_localisation_confidence_series(
    values: pd.Series,
    *,
    source_column: str,
    scale: str,
    output_column: str = LOCALISATION_CONFIDENCE_OUTPUT_COLUMN,
) -> tuple[pd.Series, LocalisationConfidenceNormalisationReport, tuple[str, ...]]:
    """Normalise localisation confidence to probabilities in ``[0.0, 1.0]``.

    Invalid and missing values are preserved as missing candidate values and
    reported. Rows are never dropped by this helper; downstream dataset
    validation remains builder-owned.
    """

    _validate_scale(scale)
    values_index = pd.Index(values.index)
    missing_mask = values.isna()
    blank_string_mask = values.map(
        lambda value: isinstance(value, str) and value.strip() == ""
    )
    bool_mask = values.map(lambda value: isinstance(value, (bool, np.bool_)))
    parse_exempt_mask = missing_mask | blank_string_mask | bool_mask
    numeric_values = pd.to_numeric(values.mask(parse_exempt_mask), errors="coerce")
    finite_mask = pd.Series(
        np.isfinite(numeric_values.to_numpy(dtype=float, copy=False, na_value=np.nan)),
        index=values_index,
    )
    upper_bound = 100.0 if scale == LOCALISATION_CONFIDENCE_SCALE_PERCENT else 1.0
    valid_numeric_mask = (
        ~parse_exempt_mask
        & numeric_values.notna()
        & finite_mask
        & numeric_values.ge(0.0)
        & numeric_values.le(upper_bound)
    )
    missing_mask = missing_mask | blank_string_mask
    invalid_mask = (~missing_mask) & (~valid_numeric_mask)
    normalized = pd.Series(pd.NA, index=values_index, dtype="Float64")
    if bool(valid_numeric_mask.any()):
        normalised_values = numeric_values.loc[valid_numeric_mask].astype(float)
        if scale == LOCALISATION_CONFIDENCE_SCALE_PERCENT:
            normalised_values = normalised_values / 100.0
        normalized.loc[valid_numeric_mask] = normalised_values
    invalid_examples: list[str] = []
    for row_id in values.index[invalid_mask].tolist()[:_EXAMPLE_LIMIT]:
        raw_value = values.at[row_id]
        parsed = _parse_confidence_value(raw_value, scale=scale)
        invalid_examples.append(f"{row_id!r}:{raw_value!r}:{parsed}")
    missing_count = int(missing_mask.sum())
    invalid_count = int(invalid_mask.sum())

    report = LocalisationConfidenceNormalisationReport(
        source_column=source_column,
        output_column=output_column,
        scale=scale,
        row_count=int(values.shape[0]),
        missing_count=missing_count,
        invalid_count=invalid_count,
        invalid_examples=tuple(invalid_examples),
    )
    warnings: list[str] = []
    if invalid_count:
        warnings.append(
            "localisation confidence contained invalid values; invalid values were "
            "preserved as missing candidate values and reported in diagnostics"
        )
    if missing_count:
        warnings.append(
            "localisation confidence contained missing values; rows were retained "
            "for builder-owned localisation policy handling"
        )
    return normalized, report, tuple(warnings)


def _parse_confidence_value(value: object, *, scale: str) -> float | str | None:
    if _is_missing(value):
        return None
    if isinstance(value, bool):
        return "bool_not_allowed"
    if isinstance(value, str):
        token = value.strip()
        if token == "":
            return None
        try:
            numeric_value = float(token)
        except ValueError:
            return "not_numeric"
    elif isinstance(value, (int, float)):
        numeric_value = float(value)
    else:
        return "unsupported_type"
    if not math.isfinite(numeric_value):
        return "not_finite"
    if scale == LOCALISATION_CONFIDENCE_SCALE_PERCENT:
        if numeric_value < 0.0 or numeric_value > 100.0:
            return "out_of_percent_range"
        return numeric_value / 100.0
    if numeric_value < 0.0 or numeric_value > 1.0:
        return "out_of_probability_range"
    return numeric_value


def _validate_scale(scale: str) -> None:
    if scale in SUPPORTED_LOCALISATION_CONFIDENCE_SCALES:
        return
    supported = ", ".join(
        repr(value) for value in SUPPORTED_LOCALISATION_CONFIDENCE_SCALES
    )
    raise PhosPyInputError(f"localisation confidence scale must be one of: {supported}")


def _summarise_examples(values: list[str], *, limit: int = _EXAMPLE_LIMIT) -> str:
    suffix = "" if len(values) <= limit else " ..."
    return "[" + ", ".join(values[:limit]) + suffix + "]"


def _is_missing(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value: object = value
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value: object = value
        return str(temporal_value) == "NaT"
    return False


__all__ = [
    "LOCALISATION_CONFIDENCE_OUTPUT_COLUMN",
    "LOCALISATION_CONFIDENCE_SCALE_PERCENT",
    "LOCALISATION_CONFIDENCE_SCALE_PROBABILITY",
    "SUPPORTED_LOCALISATION_CONFIDENCE_SCALES",
    "LocalisationConfidenceNormalisationReport",
    "LocalisationProbabilityAssessment",
    "assess_localisation_probability_column",
    "normalise_localisation_confidence_series",
    "validate_localisation_confidence_column",
    "validate_localisation_probability_column",
]
