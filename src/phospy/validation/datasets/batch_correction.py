"""Batch-correction metadata and design adequacy validation."""
# pyright: reportUnnecessaryIsInstance=false, reportUnknownMemberType=false
# Runtime boundary guards are intentionally retained for untyped external callers.

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import cast

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError


@dataclass(frozen=True, slots=True)
class ResolvedBatchDesignMetadata:
    """Batch, condition, and optional replicate labels in matrix sample order."""

    batch_by_sample: Mapping[str, str]
    condition_by_sample: Mapping[str, str]
    sample_order: tuple[str, ...]
    replicate_by_sample: Mapping[str, str] | None = None

    @property
    def batch_labels(self) -> tuple[str, ...]:
        return tuple(self.batch_by_sample[sample] for sample in self.sample_order)

    @property
    def condition_labels(self) -> tuple[str, ...]:
        return tuple(self.condition_by_sample[sample] for sample in self.sample_order)

    @property
    def replicate_labels(self) -> tuple[str, ...] | None:
        if self.replicate_by_sample is None:
            return None
        return tuple(self.replicate_by_sample[sample] for sample in self.sample_order)


class SampleMetadataAlignmentValidator:
    """Validate sample metadata without repairing or dropping rows."""

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        required_columns: Sequence[str],
        context: str = "batch-correction",
    ) -> tuple[str, ...]:
        if sample_metadata is None:
            raise PhosPyInputError(
                f"{context} validation requires sample_metadata input data"
            )
        if not isinstance(phospho, pd.DataFrame):
            raise PhosPyInputError(f"{context} validation requires phospho DataFrame")
        if not isinstance(sample_metadata, pd.DataFrame):
            raise PhosPyInputError(
                f"{context} validation requires sample_metadata DataFrame"
            )

        sample_order = _normalize_sample_index(
            phospho.columns,
            field_name=f"{context} phospho.columns",
        )
        metadata_index = _normalize_sample_index(
            sample_metadata.index,
            field_name=f"{context} sample_metadata.index",
        )
        _require_unique_sample_index(
            sample_order,
            field_name=f"{context} phospho.columns",
        )
        _require_unique_sample_index(
            metadata_index,
            field_name=f"{context} sample_metadata.index",
        )
        _require_unique_metadata_columns(sample_metadata, context=context)

        metadata_samples = set(metadata_index.tolist())
        matrix_samples = set(sample_order.tolist())
        missing_samples = [
            sample for sample in sample_order.tolist() if sample not in metadata_samples
        ]
        extra_samples = [
            sample for sample in metadata_index.tolist() if sample not in matrix_samples
        ]
        if missing_samples or extra_samples:
            details: list[str] = []
            if missing_samples:
                details.append(
                    "missing sample_metadata rows for matrix columns: "
                    f"{_format_labels(missing_samples)}"
                )
            if extra_samples:
                details.append(
                    "sample_metadata rows not present in matrix columns: "
                    f"{_format_labels(extra_samples)}"
                )
            raise PhosPyInputError(
                f"{context} sample_metadata is misaligned with matrix columns; "
                + "; ".join(details)
            )

        for column in required_columns:
            _require_metadata_column(
                sample_metadata,
                column=str(column).strip(),
                context=context,
            )
        return tuple(str(sample) for sample in sample_order.tolist())


class BatchStructureValidator:
    """Validate batch labels in aligned sample metadata."""

    def run(
        self,
        *,
        sample_metadata: pd.DataFrame,
        sample_order: Sequence[str],
        batch_column: str,
        context: str = "batch-correction",
    ) -> Mapping[str, str]:
        column = str(batch_column).strip()
        _require_metadata_column(sample_metadata, column=column, context=context)
        return _resolve_column_labels(
            sample_metadata,
            sample_order=tuple(sample_order),
            column=column,
            label_kind="batch",
            context=context,
        )


class ConditionStructureValidator:
    """Validate one or more condition columns in aligned sample metadata."""

    def run(
        self,
        *,
        sample_metadata: pd.DataFrame,
        sample_order: Sequence[str],
        condition_columns: Sequence[str],
        context: str = "batch-correction",
    ) -> Mapping[str, str]:
        if isinstance(condition_columns, str) or not isinstance(
            condition_columns, Sequence
        ):
            raise PhosPyInputError(
                f"{context} validation condition_columns must be a non-empty "
                "sequence of column names"
            )
        columns = tuple(str(column).strip() for column in condition_columns)
        if not columns or any(column == "" for column in columns):
            raise PhosPyInputError(
                f"{context} validation condition_columns must contain non-empty "
                "column names"
            )
        if len(set(columns)) != len(columns):
            raise PhosPyInputError(
                f"{context} validation condition_columns must not contain duplicates"
            )
        for column in columns:
            _require_metadata_column(sample_metadata, column=column, context=context)

        sample_tuple = tuple(sample_order)
        resolved_by_column = tuple(
            _resolve_column_labels(
                sample_metadata,
                sample_order=sample_tuple,
                column=column,
                label_kind="condition",
                context=context,
            )
            for column in columns
        )
        if len(columns) == 1:
            return resolved_by_column[0]
        return {
            sample: "|".join(
                f"{column}={labels[sample]}"
                for column, labels in zip(columns, resolved_by_column, strict=True)
            )
            for sample in sample_tuple
        }


class ReplicateStructureValidator:
    """Validate optional replicate labels when a request requires them."""

    def run(
        self,
        *,
        sample_metadata: pd.DataFrame,
        sample_order: Sequence[str],
        replicate_column: str | None,
        required: bool,
        context: str = "batch-correction",
    ) -> Mapping[str, str] | None:
        if replicate_column is None:
            if required:
                raise PhosPyInputError(
                    f"{context} validation requires replicate_column metadata"
                )
            return None
        column = str(replicate_column).strip()
        if column == "":
            raise PhosPyInputError(
                f"{context} validation replicate_column must be non-empty when provided"
            )
        _require_metadata_column(sample_metadata, column=column, context=context)
        return _resolve_column_labels(
            sample_metadata,
            sample_order=tuple(sample_order),
            column=column,
            label_kind="replicate",
            context=context,
        )


class BatchDesignMetadataValidator:
    """Coordinate sample, batch, condition, and replicate metadata validation."""

    def __init__(
        self,
        *,
        sample_alignment_validator: SampleMetadataAlignmentValidator | None = None,
        batch_structure_validator: BatchStructureValidator | None = None,
        condition_structure_validator: ConditionStructureValidator | None = None,
        replicate_structure_validator: ReplicateStructureValidator | None = None,
    ) -> None:
        self._sample_alignment_validator = (
            sample_alignment_validator or SampleMetadataAlignmentValidator()
        )
        self._batch_structure_validator = (
            batch_structure_validator or BatchStructureValidator()
        )
        self._condition_structure_validator = (
            condition_structure_validator or ConditionStructureValidator()
        )
        self._replicate_structure_validator = (
            replicate_structure_validator or ReplicateStructureValidator()
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        batch_column: str,
        condition_columns: Sequence[str],
        replicate_column: str | None = None,
        require_replicate_column: bool = False,
        context: str = "batch-correction",
    ) -> ResolvedBatchDesignMetadata:
        required_columns: list[str] = [str(batch_column).strip()]
        required_columns.extend(str(column).strip() for column in condition_columns)
        if replicate_column is not None:
            required_columns.append(str(replicate_column).strip())
        sample_order = self._sample_alignment_validator.run(
            phospho=phospho,
            sample_metadata=sample_metadata,
            required_columns=tuple(required_columns),
            context=context,
        )
        if sample_metadata is None:
            raise PhosPyInputError(
                f"{context} validation requires sample_metadata input data"
            )
        metadata_frame = sample_metadata
        aligned = cast(pd.DataFrame, metadata_frame.copy(deep=False))
        aligned.index = _normalize_sample_index(
            metadata_frame.index,
            field_name=f"{context} sample_metadata.index",
        )
        aligned = cast(pd.DataFrame, aligned.reindex(sample_order))
        batch_by_sample = self._batch_structure_validator.run(
            sample_metadata=aligned,
            sample_order=sample_order,
            batch_column=batch_column,
            context=context,
        )
        condition_by_sample = self._condition_structure_validator.run(
            sample_metadata=aligned,
            sample_order=sample_order,
            condition_columns=condition_columns,
            context=context,
        )
        replicate_by_sample = self._replicate_structure_validator.run(
            sample_metadata=aligned,
            sample_order=sample_order,
            replicate_column=replicate_column,
            required=require_replicate_column,
            context=context,
        )
        return ResolvedBatchDesignMetadata(
            batch_by_sample=batch_by_sample,
            condition_by_sample=condition_by_sample,
            replicate_by_sample=replicate_by_sample,
            sample_order=sample_order,
        )


class DesignRankValidator:
    """Validate full column rank for numeric design matrices."""

    def run(self, design_matrix: pd.DataFrame, *, context: str) -> int:
        if not isinstance(design_matrix, pd.DataFrame):
            raise PhosPyInputError(f"{context} design matrix must be a DataFrame")
        if not design_matrix.index.is_unique:
            raise PhosPyInputError(
                f"{context} design matrix sample labels must be unique"
            )
        if not design_matrix.columns.is_unique:
            raise PhosPyInputError(
                f"{context} design matrix coefficient labels must be unique"
            )
        try:
            values = np.asarray(design_matrix.to_numpy(dtype=float), dtype=float)
        except (TypeError, ValueError) as exc:
            raise PhosPyInputError(
                f"{context} design matrix must contain numeric values"
            ) from exc
        if not np.isfinite(values).all():
            raise PhosPyInputError(
                f"{context} design matrix must contain finite numeric values"
            )
        rank = _matrix_rank(values)
        column_count = int(values.shape[1])
        if rank < column_count:
            coefficients = ", ".join(str(label) for label in design_matrix.columns)
            raise PhosPyInputError(
                f"{context} design matrix is rank-deficient; condition, batch, "
                "or replicate terms are collinear or confounded "
                f"(rank={rank}, columns={column_count}, coefficients={coefficients})"
            )
        return rank


class BatchCorrectionAdequacyValidator:
    """Validate fixed-effect batch-correction design adequacy.

    The validator only inspects sample labels and categorical design rank. It does
    not correct matrices, estimate coefficients, or mutate metadata.
    """

    def __init__(
        self,
        *,
        design_rank_validator: DesignRankValidator | None = None,
    ) -> None:
        self._design_rank_validator = design_rank_validator or DesignRankValidator()

    def run(
        self,
        *,
        batch_by_sample: Mapping[str, object],
        condition_by_sample: Mapping[str, object],
        sample_order: Sequence[str],
        preserve_condition_effects: bool,
    ) -> None:
        samples = _normalize_sample_order(sample_order)
        if preserve_condition_effects is not True:
            raise PhosPyInputError(
                "linear_residualize_batch requires "
                "preprocessing_config.batch_correction.preserve_condition_effects=True; "
                "refusing batch correction because condition effects would not be "
                "explicitly preserved"
            )
        if len(samples) < 2:
            raise PhosPyInputError(
                "linear_residualize_batch requires at least two samples to "
                "estimate batch effects while preserving condition effects"
            )

        batch_labels = _resolve_labels(
            batch_by_sample,
            sample_order=samples,
            label_kind="batch",
        )
        condition_labels = _resolve_labels(
            condition_by_sample,
            sample_order=samples,
            label_kind="condition",
        )
        batch_levels = _levels_in_order(batch_labels)
        condition_levels = _levels_in_order(condition_labels)
        if len(batch_levels) < 2:
            raise PhosPyInputError(
                "linear_residualize_batch requires at least two batch levels; "
                f"observed {len(batch_levels)}"
            )

        singleton_batches = _singleton_levels(batch_labels)
        if singleton_batches:
            raise PhosPyInputError(
                "linear_residualize_batch requires at least two samples in each "
                "batch level to estimate batch effects; singleton batch levels: "
                f"{_format_labels(singleton_batches)}"
            )

        preservation_design = _treatment_coded_design(
            condition_labels,
            include_intercept=True,
        )
        preservation_columns = int(preservation_design.shape[1])
        preservation_rank = _matrix_rank(preservation_design)
        if preservation_rank < preservation_columns:
            raise PhosPyInputError(
                "linear_residualize_batch condition preservation design is "
                "rank-deficient; condition effects cannot be explicitly preserved "
                f"(rank={preservation_rank}, columns={preservation_columns})"
            )
        if len(samples) <= preservation_rank:
            raise PhosPyInputError(
                "linear_residualize_batch condition preservation design is "
                "saturated; batch effects cannot be estimated while preserving "
                "condition effects "
                f"(samples={len(samples)}, condition_design_rank={preservation_rank})"
            )

        batch_terms = _treatment_coded_design(batch_labels, include_intercept=False)
        full_design = np.concatenate((preservation_design, batch_terms), axis=1)
        full_columns = int(full_design.shape[1])
        if len(samples) <= full_columns:
            raise PhosPyInputError(
                "linear_residualize_batch requires more samples than estimable "
                "condition-plus-batch design parameters; "
                f"samples={len(samples)}, design_columns={full_columns}. Add "
                "replicate samples or reduce batch/condition levels."
            )

        if len(condition_levels) > 1 and _condition_is_determined_by_batch(
            batch_labels=batch_labels,
            condition_labels=condition_labels,
        ):
            raise PhosPyInputError(
                "linear_residualize_batch cannot run because batch and condition "
                "are perfectly confounded: each batch level contains only one "
                "condition level, so removing batch would remove biological "
                "condition signal"
            )
        if len(condition_levels) > 1 and _batch_is_determined_by_condition(
            batch_labels=batch_labels,
            condition_labels=condition_labels,
        ):
            raise PhosPyInputError(
                "linear_residualize_batch cannot run because batch and condition "
                "are perfectly confounded: each condition level contains only one "
                "batch level, so batch cannot be estimated while preserving "
                "condition effects"
            )

        full_frame = pd.DataFrame(
            full_design,
            index=pd.Index(samples, name="sample"),
            columns=pd.Index(
                (
                    *(f"condition[{index}]" for index in range(preservation_columns)),
                    *(f"batch[{index}]" for index in range(batch_terms.shape[1])),
                ),
                name="coefficient",
            ),
        )
        try:
            full_rank = self._design_rank_validator.run(
                full_frame,
                context="linear_residualize_batch batch/condition",
            )
        except PhosPyInputError as exc:
            raise PhosPyInputError(
                "linear_residualize_batch batch/condition design is "
                "rank-deficient; batch effects are not estimable while preserving "
                "condition effects "
                f"(rank={_matrix_rank(full_design)}, columns={full_columns})"
            ) from exc
        batch_degrees = int(batch_terms.shape[1])
        batch_rank_after_condition = full_rank - preservation_rank
        if full_rank < full_columns or batch_rank_after_condition < batch_degrees:
            raise PhosPyInputError(
                "linear_residualize_batch batch/condition design is "
                "rank-deficient; batch effects are not estimable while preserving "
                "condition effects "
                f"(rank={full_rank}, columns={full_columns}, "
                f"estimable_batch_degrees={batch_rank_after_condition}, "
                f"batch_degrees={batch_degrees})"
            )


def _normalize_sample_order(sample_order: Sequence[str]) -> tuple[str, ...]:
    samples = tuple(str(sample).strip() for sample in sample_order)
    blank_positions = [
        position for position, sample in enumerate(samples) if sample == ""
    ]
    if blank_positions:
        raise PhosPyInputError(
            "linear_residualize_batch sample_order contains blank sample labels at "
            f"positions {_format_positions(blank_positions)}"
        )
    if len(set(samples)) != len(samples):
        duplicates = list(
            dict.fromkeys(sample for sample in samples if samples.count(sample) > 1)
        )
        raise PhosPyInputError(
            "linear_residualize_batch sample_order contains duplicate sample "
            f"labels: {_format_labels(duplicates)}"
        )
    return samples


def _normalize_sample_index(index: pd.Index, *, field_name: str) -> pd.Index:
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
            f"{field_name} must not contain missing sample labels; positions: "
            f"{_format_positions(missing_positions)}"
        )
    if blank_positions:
        raise PhosPyInputError(
            f"{field_name} must not contain blank sample labels; positions: "
            f"{_format_positions(blank_positions)}"
        )
    return pd.Index(normalized, name=index.name)


def _require_unique_sample_index(index: pd.Index, *, field_name: str) -> None:
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
        f"{field_name} contains duplicate sample labels: {_format_labels(duplicates)}"
    )


def _require_unique_metadata_columns(
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
        f"{_format_labels(duplicate_columns)}"
    )


def _require_metadata_column(
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


def _resolve_column_labels(
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
        if _is_missing_value(value):
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
            f"{label_kind} labels for samples: {_format_labels(missing_samples)}"
        )
    if blank_samples:
        raise PhosPyInputError(
            f"{context} sample_metadata column {column!r} contains blank "
            f"{label_kind} labels for samples: {_format_labels(blank_samples)}"
        )
    return labels


def _resolve_labels(
    labels_by_sample: Mapping[str, object],
    *,
    sample_order: tuple[str, ...],
    label_kind: str,
) -> tuple[str, ...]:
    labels: list[str] = []
    missing_samples: list[str] = []
    blank_samples: list[str] = []
    for sample in sample_order:
        if sample not in labels_by_sample:
            missing_samples.append(sample)
            continue
        value = labels_by_sample[sample]
        if _is_missing_value(value):
            missing_samples.append(sample)
            continue
        label = str(value).strip()
        if label == "":
            blank_samples.append(sample)
            continue
        labels.append(label)

    if missing_samples:
        raise PhosPyInputError(
            f"linear_residualize_batch requires {label_kind} labels for every "
            f"sample; missing {label_kind} labels for samples: "
            f"{_format_labels(missing_samples)}"
        )
    if blank_samples:
        raise PhosPyInputError(
            f"linear_residualize_batch requires {label_kind} labels for every "
            f"sample; blank {label_kind} labels for samples: "
            f"{_format_labels(blank_samples)}"
        )
    return tuple(labels)


def _levels_in_order(labels: Sequence[str]) -> tuple[str, ...]:
    levels: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label in seen:
            continue
        seen.add(label)
        levels.append(label)
    return tuple(levels)


def _singleton_levels(labels: Sequence[str]) -> tuple[str, ...]:
    counts: dict[str, int] = {}
    for label in labels:
        counts[label] = counts.get(label, 0) + 1
    return tuple(label for label in _levels_in_order(labels) if counts[label] == 1)


def _treatment_coded_design(
    labels: Sequence[str],
    *,
    include_intercept: bool,
) -> np.ndarray:
    levels = _levels_in_order(labels)
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


def _matrix_rank(matrix: np.ndarray) -> int:
    if matrix.size == 0:
        return 0
    return int(np.linalg.matrix_rank(matrix))


def _condition_is_determined_by_batch(
    *,
    batch_labels: Sequence[str],
    condition_labels: Sequence[str],
) -> bool:
    conditions_by_batch: dict[str, set[str]] = {}
    for batch, condition in zip(batch_labels, condition_labels, strict=True):
        conditions_by_batch.setdefault(batch, set()).add(condition)
    return all(len(conditions) == 1 for conditions in conditions_by_batch.values())


def _batch_is_determined_by_condition(
    *,
    batch_labels: Sequence[str],
    condition_labels: Sequence[str],
) -> bool:
    batches_by_condition: dict[str, set[str]] = {}
    for batch, condition in zip(batch_labels, condition_labels, strict=True):
        batches_by_condition.setdefault(condition, set()).add(batch)
    return all(len(batches) == 1 for batches in batches_by_condition.values())


def _is_missing_value(value: object) -> bool:
    return bool(pd.Series((value,), dtype="object").isna().iat[0])


def _format_positions(positions: Sequence[int]) -> str:
    preview = ", ".join(str(position) for position in tuple(positions)[:8])
    suffix = "" if len(positions) <= 8 else ", ..."
    return f"[{preview}{suffix}]"


def _format_labels(labels: Sequence[str]) -> str:
    preview = ", ".join(repr(value) for value in tuple(labels)[:5])
    suffix = "" if len(labels) <= 5 else " ..."
    return f"{preview}{suffix}"


__all__ = [
    "BatchCorrectionAdequacyValidator",
    "BatchDesignMetadataValidator",
    "BatchStructureValidator",
    "ConditionStructureValidator",
    "DesignRankValidator",
    "ReplicateStructureValidator",
    "ResolvedBatchDesignMetadata",
    "SampleMetadataAlignmentValidator",
]
