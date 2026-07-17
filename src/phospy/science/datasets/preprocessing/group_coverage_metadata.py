"""Group-aware coverage-filter metadata resolution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError


@dataclass(frozen=True, slots=True)
class ResolvedGroupCoverageFilterMetadata:
    """Sample groups aligned to phospho matrix sample order."""

    group_column: str
    sample_order: tuple[str, ...]
    group_by_sample: Mapping[str, str]
    sample_order_by_group: Mapping[str, tuple[str, ...]]

    @property
    def group_labels(self) -> tuple[str, ...]:
        return tuple(self.sample_order_by_group.keys())


class GroupCoverageFilterMetadataValidator:
    """Validate and resolve sample groups for coverage filtering."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        group_column: str | None,
        min_groups_passing_threshold: int,
    ) -> ResolvedGroupCoverageFilterMetadata:
        if sample_metadata is None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.group_coverage_filter "
                "requires sample_metadata input data"
            )
        if group_column is None or not str(group_column).strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.group_coverage_filter."
                "group_column must be a non-empty string"
            )
        resolved_group_column = str(group_column).strip()
        sample_order = _normalize_label_index(
            phospho.columns,
            field_name="dataset build request phospho.columns",
        )
        metadata_index = _normalize_label_index(
            sample_metadata.index,
            field_name="dataset build request sample_metadata.index",
        )
        _require_unique_labels(
            sample_order,
            field_name="dataset build request phospho.columns",
        )
        _require_unique_labels(
            metadata_index,
            field_name="dataset build request sample_metadata.index",
        )
        _require_column(sample_metadata, column=resolved_group_column)
        _require_no_missing_or_extra_metadata_rows(
            metadata_index=metadata_index,
            sample_order=sample_order,
        )

        group_by_sample = _resolve_groups_by_sample(
            sample_metadata=sample_metadata,
            metadata_index=metadata_index,
            column=resolved_group_column,
            sample_order=sample_order,
        )
        sample_order_by_group = _group_samples(
            group_by_sample=group_by_sample,
            sample_order=tuple(sample_order.tolist()),
        )
        if min_groups_passing_threshold > len(sample_order_by_group):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.group_coverage_filter."
                "min_groups_passing_threshold cannot exceed the number of observed "
                f"groups; threshold={int(min_groups_passing_threshold)}, "
                f"observed_groups={int(len(sample_order_by_group))}"
            )
        return ResolvedGroupCoverageFilterMetadata(
            group_column=resolved_group_column,
            sample_order=tuple(sample_order.tolist()),
            group_by_sample=group_by_sample,
            sample_order_by_group=sample_order_by_group,
        )


def require_numeric_group_coverage_matrix(phospho: pd.DataFrame) -> None:
    """Require numeric, non-boolean phospho columns for finite-value counting."""

    invalid_columns = [
        str(column)
        for column in phospho.columns
        if (
            not pd.api.types.is_numeric_dtype(phospho[column])
            or pd.api.types.is_bool_dtype(phospho[column])
        )
    ]
    if invalid_columns:
        preview = ", ".join(repr(value) for value in invalid_columns[:5])
        suffix = "" if len(invalid_columns) <= 5 else " ..."
        raise PhosPyInputError(
            "dataset build request preprocessing_config.group_coverage_filter "
            "requires numeric phospho sample columns for finite-value counting. "
            f"Non-numeric columns: {preview}{suffix}"
        )


def _normalize_label_index(index: pd.Index, *, field_name: str) -> pd.Index:
    normalized: list[str] = []
    missing_positions: list[int] = []
    blank_positions: list[int] = []
    for position, value in enumerate(index.tolist()):
        if _is_missing_value(value):
            missing_positions.append(position)
            continue
        label = str(value).strip()
        if label == "":
            blank_positions.append(position)
            continue
        normalized.append(label)

    if missing_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain missing labels for group-aware "
            f"coverage filtering; found positions: "
            f"{_format_positions(missing_positions)}"
        )
    if blank_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain blank labels for group-aware coverage "
            f"filtering; found positions: {_format_positions(blank_positions)}"
        )
    return pd.Index(normalized, name=index.name)


def _require_unique_labels(index: pd.Index, *, field_name: str) -> None:
    if index.is_unique:
        return
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in index.tolist():
        label = str(value)
        if label in seen and label not in duplicates:
            duplicates.append(label)
        seen.add(label)
    preview = ", ".join(repr(value) for value in duplicates[:5])
    suffix = "" if len(duplicates) <= 5 else " ..."
    raise PhosPyInputError(
        f"{field_name} contains duplicate sample labels for group-aware coverage "
        f"filtering: {preview}{suffix}"
    )


def _require_column(sample_metadata: pd.DataFrame, *, column: str) -> None:
    matches = [
        candidate for candidate in sample_metadata.columns if candidate == column
    ]
    if len(matches) == 1:
        return
    if len(matches) > 1:
        raise PhosPyInputError(
            "dataset build request preprocessing_config.group_coverage_filter."
            f"group_column resolves to duplicate sample_metadata column {column!r}"
        )
    raise PhosPyInputError(
        "dataset build request preprocessing_config.group_coverage_filter."
        f"group_column references missing sample_metadata column {column!r}"
    )


def _require_no_missing_or_extra_metadata_rows(
    *,
    metadata_index: pd.Index,
    sample_order: pd.Index,
) -> None:
    metadata_samples = set(metadata_index.tolist())
    matrix_samples = set(sample_order.tolist())
    missing_samples = [
        sample for sample in sample_order.tolist() if sample not in metadata_samples
    ]
    extra_samples = [
        sample for sample in metadata_index.tolist() if sample not in matrix_samples
    ]
    if not missing_samples and not extra_samples:
        return

    details: list[str] = []
    if missing_samples:
        details.append(
            "missing rows for coverage-filter samples: "
            f"{_format_label_preview(missing_samples)}"
        )
    if extra_samples:
        details.append(
            "rows not present in phospho columns for coverage-filter metadata "
            f"resolution: {_format_label_preview(extra_samples)}"
        )
    raise PhosPyInputError(
        "dataset build request sample_metadata has incompatible sample rows for "
        f"group-aware coverage filtering; {'; '.join(details)}"
    )


def _resolve_groups_by_sample(
    *,
    sample_metadata: pd.DataFrame,
    metadata_index: pd.Index,
    column: str,
    sample_order: pd.Index,
) -> dict[str, str]:
    groups_by_metadata_sample: dict[str, str] = {}
    missing_samples: list[str] = []
    blank_samples: list[str] = []
    column_values = sample_metadata[column].tolist()
    for sample, value in zip(
        metadata_index.tolist(),
        column_values,
        strict=True,
    ):
        if _is_missing_value(value):
            missing_samples.append(str(sample))
            continue
        label = str(value).strip()
        if label == "":
            blank_samples.append(str(sample))
            continue
        groups_by_metadata_sample[str(sample)] = label

    if missing_samples:
        raise PhosPyInputError(
            f"dataset build request sample_metadata column {column!r} contains "
            "missing group labels for samples: "
            f"{_format_label_preview(missing_samples)}"
        )
    if blank_samples:
        raise PhosPyInputError(
            f"dataset build request sample_metadata column {column!r} contains "
            f"blank group labels for samples: {_format_label_preview(blank_samples)}"
        )
    return {
        str(sample): groups_by_metadata_sample[str(sample)]
        for sample in sample_order.tolist()
    }


def _group_samples(
    *,
    group_by_sample: Mapping[str, str],
    sample_order: tuple[str, ...],
) -> dict[str, tuple[str, ...]]:
    grouped: dict[str, list[str]] = {}
    for sample in sample_order:
        group = group_by_sample[sample]
        grouped.setdefault(group, []).append(sample)
    return {group: tuple(samples) for group, samples in grouped.items()}


def _is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _format_positions(positions: list[int]) -> str:
    preview = ", ".join(str(position) for position in positions[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def _format_label_preview(labels: list[str]) -> str:
    preview = ", ".join(repr(value) for value in labels[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"


__all__ = [
    "GroupCoverageFilterMetadataValidator",
    "ResolvedGroupCoverageFilterMetadata",
    "require_numeric_group_coverage_matrix",
]
