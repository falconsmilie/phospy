"""Batch-correction design adequacy validation."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from phospy.errors.input import PhosPyInputError


class BatchCorrectionAdequacyValidator:
    """Validate fixed-effect batch-correction design adequacy.

    The validator only inspects sample labels and categorical design rank. It does
    not correct matrices, estimate coefficients, or mutate metadata.
    """

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

        full_rank = _matrix_rank(full_design)
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


__all__ = ["BatchCorrectionAdequacyValidator"]
