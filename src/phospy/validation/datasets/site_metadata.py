"""Phosphosite metadata validation and workflow policy enforcement."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, TypeVar, cast

import numpy as np
import pandas as pd

from phospy.science.sites.metadata_validation import (
    enforce_site_identity_rows,
    validate_site_identity_metadata,
    validate_site_sequence_column,
)
from phospy.science.sites.sequence_context import (
    SequenceContextContract,
    enforce_centred_site_sequence_context,
    enforce_site_sequence_context_contract,
)

if TYPE_CHECKING:
    from phospy.contracts.configs.localisation import LocalisationRequirement

ErrorType = TypeVar("ErrorType", bound=Exception)
_EXAMPLE_LIMIT = 5


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


def validate_localisation_probability_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "localisation_probability",
) -> None:
    """Validate optional localisation probability values when the column exists."""

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


def enforce_localisation_requirement(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    workflow_name: str,
    requirement: LocalisationRequirement,
    error_type: type[ErrorType],
    column_name: str = "localisation_confidence",
) -> None:
    """Enforce workflow-level localisation policy using row-context diagnostics."""

    if not requirement.requires_probability_column:
        return
    resolved_column_name = _resolve_localisation_column_name(
        site_metadata=site_metadata,
        requested_column_name=column_name,
    )
    if resolved_column_name is None:
        site_examples = _site_id_examples(site_metadata.index)
        raise error_type(
            f"{workflow_name} requires localisation metadata policy={requirement.policy}; "
            f"missing required column={field_name}.{column_name}; "
            f"affected_rows={int(site_metadata.shape[0])}; "
            f"example_site_ids={site_examples}"
        )

    assessment = assess_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=resolved_column_name,
    )
    if assessment is None:  # pragma: no cover - defensive guard
        return
    if assessment.invalid_count > 0:
        invalid_sites = _site_id_examples(
            _index_by_boolean_mask(site_metadata.index, assessment.invalid_mask)
        )
        raise error_type(
            f"{workflow_name} requires localisation metadata policy={requirement.policy}; "
            f"invalid values in {field_name}.{resolved_column_name}; "
            f"affected_rows={assessment.invalid_count}; "
            f"example_site_ids={invalid_sites}; "
            f"example_values={_summarise_examples(list(assessment.invalid_examples), limit=3)}"
        )
    if requirement.require_present and assessment.missing_count > 0:
        missing_sites = _site_id_examples(
            _index_by_boolean_mask(site_metadata.index, assessment.missing_mask)
        )
        raise error_type(
            f"{workflow_name} requires localisation metadata policy={requirement.policy}; "
            f"missing values in {field_name}.{resolved_column_name}; "
            f"affected_rows={assessment.missing_count}; "
            f"example_site_ids={missing_sites}"
        )
    if requirement.minimum_probability is None:
        return
    below_threshold = assessment.normalized.notna() & (
        assessment.normalized.astype("float64") < requirement.minimum_probability
    )
    below_threshold_count = int(below_threshold.sum())
    if below_threshold_count <= 0:
        return
    below_threshold_sites = _site_id_examples(
        _index_by_boolean_mask(site_metadata.index, below_threshold)
    )
    threshold = float(requirement.minimum_probability)
    raise error_type(
        f"{workflow_name} requires localisation metadata policy={requirement.policy}; "
        f"{field_name}.{resolved_column_name} must be >= {threshold:.3f}; "
        f"affected_rows={below_threshold_count}; "
        f"example_site_ids={below_threshold_sites}"
    )


def enforce_required_non_empty_string_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    workflow_name: str,
    column_name: str,
    error_type: type[ErrorType],
) -> None:
    """Require one site-metadata column to be non-missing, non-empty strings."""

    if column_name not in site_metadata.columns:
        site_examples = _site_id_examples(site_metadata.index)
        raise error_type(
            f"{field_name} is missing required columns: {column_name}; "
            f"{workflow_name} requires {field_name}.{column_name}; "
            f"missing required column={field_name}.{column_name}; "
            f"affected_rows={int(site_metadata.shape[0])}; "
            f"example_site_ids={site_examples}"
        )
    column = site_metadata[column_name]
    invalid_mask = pd.Series(False, index=pd.Index(column.index), dtype="boolean")
    for site_id, raw_value in column.items():
        if _is_missing(raw_value):
            invalid_mask.at[site_id] = True
            continue
        if not isinstance(raw_value, str):
            invalid_mask.at[site_id] = True
            continue
        if raw_value.strip() == "":
            invalid_mask.at[site_id] = True
    invalid_count = int(invalid_mask.sum())
    if invalid_count == 0:
        return
    invalid_sites = _site_id_examples(
        _index_by_boolean_mask(site_metadata.index, invalid_mask)
    )
    raise error_type(
        f"{workflow_name} requires {field_name}.{column_name} to contain non-empty "
        f"string values; affected_rows={invalid_count}; "
        f"example_site_ids={invalid_sites}"
    )


def assess_localisation_probability_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
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
        parsed = _parse_localisation_probability(raw_value)
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


def validate_localisation_confidence_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "localisation_confidence",
) -> None:
    """Validate optional localisation-confidence values when the column exists."""

    validate_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=column_name,
    )


def assess_localisation_confidence_column(
    *,
    site_metadata: pd.DataFrame,
    field_name: str,
    error_type: type[ErrorType],
    column_name: str = "localisation_confidence",
) -> LocalisationProbabilityAssessment | None:
    """Parse optional localisation-confidence values with diagnostics."""

    return assess_localisation_probability_column(
        site_metadata=site_metadata,
        field_name=field_name,
        error_type=error_type,
        column_name=column_name,
    )


def _parse_localisation_probability(value: object) -> float | str | None:
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
    if numeric_value < 0.0 or numeric_value > 1.0:
        return "out_of_range"
    return float(numeric_value)


def _resolve_localisation_column_name(
    *,
    site_metadata: pd.DataFrame,
    requested_column_name: str,
) -> str | None:
    if requested_column_name in site_metadata.columns:
        return requested_column_name
    if (
        requested_column_name == "localisation_confidence"
        and "localisation_probability" in site_metadata.columns
    ):
        return "localisation_probability"
    return None


def _site_id_examples(index: pd.Index, *, limit: int = _EXAMPLE_LIMIT) -> str:
    labels = [str(value) for value in index.tolist()]
    if not labels:
        return "(none)"
    preview = ", ".join(repr(label) for label in labels[:limit])
    suffix = "" if len(labels) <= limit else f", +{len(labels) - limit} more"
    return f"[{preview}{suffix}]"


def _index_by_boolean_mask(index: pd.Index, mask: pd.Series) -> pd.Index:
    labels = index.tolist()
    mask_values = cast(list[object], mask.tolist())
    selected = [
        label
        for label, include in zip(labels, mask_values, strict=True)
        if bool(include)
    ]
    return pd.Index(selected)


def _summarise_examples(values: list[str], *, limit: int = _EXAMPLE_LIMIT) -> str:
    if not values:
        return "(none)"
    preview = ", ".join(values[:limit])
    suffix = "" if len(values) <= limit else f", +{len(values) - limit} more"
    return f"[{preview}{suffix}]"


def _is_missing(value: object) -> bool:
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, np.floating):
        scalar_value = cast(object, value)
        return str(scalar_value).lower() == "nan"
    if isinstance(value, (np.datetime64, np.timedelta64)):
        temporal_value = cast(object, value)
        return str(temporal_value) == "NaT"
    return False


__all__ = [
    "LocalisationProbabilityAssessment",
    "SequenceContextContract",
    "assess_localisation_confidence_column",
    "assess_localisation_probability_column",
    "enforce_centred_site_sequence_context",
    "enforce_site_sequence_context_contract",
    "enforce_site_identity_rows",
    "enforce_required_non_empty_string_column",
    "enforce_localisation_requirement",
    "validate_localisation_confidence_column",
    "validate_site_sequence_column",
    "validate_localisation_probability_column",
    "validate_site_identity_metadata",
]
