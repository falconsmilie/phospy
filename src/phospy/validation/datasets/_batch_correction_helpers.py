"""Shared private helpers for batch-correction dataset validation."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.input import PhosPyInputError

_FloatArray = npt.NDArray[np.float64]

NOT_PROVIDED_VALUES = frozenset({"not_provided", "not provided"})
MISSING_ENVIRONMENT_VALUES = NOT_PROVIDED_VALUES | frozenset({"unknown"})
MISSING = object()


def require_non_empty_mapping(
    value: object,
    *,
    field_name: str,
) -> None:
    if not isinstance(value, Mapping) or not value:
        raise PhosPyInputError(f"{field_name} must be a non-empty object")


def reject_not_provided_required_mapping(
    value: Mapping[str, object],
    *,
    field_name: str,
) -> None:
    for key, item in value.items():
        if is_not_provided(item):
            raise PhosPyInputError(
                f"{field_name}.{key} must not be recorded as not_provided for "
                "applied SPS/RUV-style corrected output"
            )


def has_non_missing_text(value: object) -> bool:
    return not is_missing_required_text(value) and not is_not_provided(value)


def is_missing_environment_text(value: object) -> bool:
    return (
        is_missing_required_text(value)
        or str(value).strip().lower() in MISSING_ENVIRONMENT_VALUES
    )


def is_missing_required_text(value: object) -> bool:
    return value is None or str(value).strip() == ""


def is_not_provided(value: object) -> bool:
    return str(value).strip().lower() in NOT_PROVIDED_VALUES


def normalize_sample_order(
    sample_order: Sequence[str],
    *,
    context: str = "linear_residualize_batch",
) -> tuple[str, ...]:
    samples = tuple(str(sample).strip() for sample in sample_order)
    blank_positions = [
        position for position, sample in enumerate(samples) if sample == ""
    ]
    if blank_positions:
        raise PhosPyInputError(
            f"{context} sample_order contains blank sample labels at "
            f"positions {format_positions(blank_positions)}"
        )
    if len(set(samples)) != len(samples):
        duplicates = list(
            dict.fromkeys(sample for sample in samples if samples.count(sample) > 1)
        )
        raise PhosPyInputError(
            f"{context} sample_order contains duplicate sample "
            f"labels: {format_labels(duplicates)}"
        )
    return samples


def normalize_sample_index(index: pd.Index, *, field_name: str) -> pd.Index:
    normalized: list[str] = []
    missing_positions: list[int] = []
    blank_positions: list[int] = []
    for position, value in enumerate(index.tolist()):
        if is_missing_value(value):
            missing_positions.append(position)
            continue
        label = str(value).strip()
        if label == "":
            blank_positions.append(position)
            continue
        normalized.append(label)
    if missing_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain missing sample labels; positions: "
            f"{format_positions(missing_positions)}"
        )
    if blank_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain blank sample labels; positions: "
            f"{format_positions(blank_positions)}"
        )
    return pd.Index(normalized, name=index.name)


def require_unique_sample_index(index: pd.Index, *, field_name: str) -> None:
    if index.is_unique:
        return
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in list(index):
        label = str(value)
        if label in seen and label not in duplicates:
            duplicates.append(label)
        seen.add(label)
    raise PhosPyInputError(
        f"{field_name} contains duplicate sample labels: {format_labels(duplicates)}"
    )


def require_unique_metadata_columns(
    sample_metadata: pd.DataFrame,
    *,
    context: str,
) -> None:
    if sample_metadata.columns.is_unique:
        return
    seen: set[str] = set()
    duplicate_columns: list[str] = []
    for value in list(sample_metadata.columns):
        column = str(value)
        if column in seen and column not in duplicate_columns:
            duplicate_columns.append(column)
        seen.add(column)
    raise PhosPyInputError(
        f"{context} sample_metadata contains duplicated columns: "
        f"{format_labels(duplicate_columns)}"
    )


def require_metadata_column(
    sample_metadata: pd.DataFrame,
    *,
    column: str,
    context: str,
) -> None:
    if column == "":
        raise PhosPyInputError(f"{context} metadata column names must be non-empty")
    matches = [
        candidate for candidate in sample_metadata.columns if candidate == column
    ]
    if len(matches) == 1:
        return
    if len(matches) > 1:
        raise PhosPyInputError(
            f"{context} sample_metadata contains duplicated column {column!r}"
        )
    raise PhosPyInputError(
        f"{context} sample_metadata is missing required column {column!r}"
    )


def resolve_column_labels(
    sample_metadata: pd.DataFrame,
    *,
    sample_order: tuple[str, ...],
    column: str,
    label_kind: str,
    context: str,
) -> dict[str, str]:
    labels: dict[str, str] = {}
    missing_samples: list[str] = []
    blank_samples: list[str] = []
    for sample in sample_order:
        value = sample_metadata.at[sample, column]
        if is_missing_value(value):
            missing_samples.append(sample)
            continue
        label = str(value).strip()
        if label == "":
            blank_samples.append(sample)
            continue
        labels[sample] = label
    if missing_samples:
        raise PhosPyInputError(
            f"{context} sample_metadata column {column!r} contains missing "
            f"{label_kind} labels for samples: {format_labels(missing_samples)}"
        )
    if blank_samples:
        raise PhosPyInputError(
            f"{context} sample_metadata column {column!r} contains blank "
            f"{label_kind} labels for samples: {format_labels(blank_samples)}"
        )
    return labels


def resolve_labels(
    labels_by_sample: Mapping[str, object],
    *,
    sample_order: tuple[str, ...],
    label_kind: str,
    context: str = "linear_residualize_batch",
) -> tuple[str, ...]:
    labels: list[str] = []
    missing_samples: list[str] = []
    blank_samples: list[str] = []
    for sample in sample_order:
        if sample not in labels_by_sample:
            missing_samples.append(sample)
            continue
        value = labels_by_sample[sample]
        if is_missing_value(value):
            missing_samples.append(sample)
            continue
        label = str(value).strip()
        if label == "":
            blank_samples.append(sample)
            continue
        labels.append(label)

    if missing_samples:
        raise PhosPyInputError(
            f"{context} requires {label_kind} labels for every "
            f"sample; missing {label_kind} labels for samples: "
            f"{format_labels(missing_samples)}"
        )
    if blank_samples:
        raise PhosPyInputError(
            f"{context} requires {label_kind} labels for every "
            f"sample; blank {label_kind} labels for samples: "
            f"{format_labels(blank_samples)}"
        )
    return tuple(labels)


def levels_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    levels: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        levels.append(label)
    return tuple(levels)


def unique_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    return levels_in_order(labels)


def duplicates_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    return tuple(label for label in levels_in_order(labels) if counts[label] > 1)


def singleton_levels(labels: Sequence[str]) -> tuple[str, ...]:
    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    return tuple(label for label in levels_in_order(labels) if counts[label] == 1)


def treatment_coded_design(
    labels: Sequence[str],
    *,
    include_intercept: bool,
) -> _FloatArray:
    levels = levels_in_order(labels)
    row_width = (1 if include_intercept else 0) + max(len(levels) - 1, 0)
    if row_width == 0:
        return np.empty((len(labels), 0), dtype=float)

    rows: list[list[float]] = []
    for label in labels:
        row: list[float] = []
        if include_intercept:
            row.append(1.0)
        row.extend(1.0 if label == level else 0.0 for level in levels[1:])
        rows.append(row)
    return np.asarray(rows, dtype=float)


def matrix_rank(matrix: _FloatArray) -> int:
    if matrix.size == 0:
        return 0
    return int(np.linalg.matrix_rank(matrix))


def condition_is_determined_by_batch(
    *,
    batch_labels: Sequence[str],
    condition_labels: Sequence[str],
) -> bool:
    conditions_by_batch: dict[str, set[str]] = {}
    for batch, condition in zip(batch_labels, condition_labels, strict=True):
        conditions_by_batch.setdefault(batch, set()).add(condition)
    return all(len(conditions) == 1 for conditions in conditions_by_batch.values())


def batch_is_determined_by_condition(
    *,
    batch_labels: Sequence[str],
    condition_labels: Sequence[str],
) -> bool:
    batches_by_condition: dict[str, set[str]] = {}
    for batch, condition in zip(batch_labels, condition_labels, strict=True):
        batches_by_condition.setdefault(condition, set()).add(batch)
    return all(len(batches) == 1 for batches in batches_by_condition.values())


def same_label_partition(
    left_labels: Sequence[str],
    right_labels: Sequence[str],
) -> bool:
    left = tuple(left_labels)
    right = tuple(right_labels)
    if len(left) != len(right) or not left:
        return False
    for first_index in range(len(left)):
        for second_index in range(len(left)):
            same_left = left[first_index] == left[second_index]
            same_right = right[first_index] == right[second_index]
            if same_left != same_right:
                return False
    return True


def is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def format_positions(positions: Sequence[int]) -> str:
    preview = ", ".join(str(position) for position in tuple(positions)[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def format_labels(labels: Sequence[str]) -> str:
    preview = ", ".join(repr(value) for value in tuple(labels)[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"
