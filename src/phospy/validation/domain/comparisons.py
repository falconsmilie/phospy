from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

from ...errors import InputCompatibilityError
from ...internal.constants import ComparisonSpec

if TYPE_CHECKING:
    from ...datasets.schema import DatasetSchema


def validate_comparison_specs(
    *,
    comparison_groups: Sequence[str],
    comparisons: Sequence[ComparisonSpec] | None,
    context: str,
) -> tuple[ComparisonSpec, ...] | None:
    """Validate comparison definitions against a set of allowed groups."""

    if comparisons is None:
        return None

    resolved = tuple(comparisons)
    valid_groups = frozenset(comparison_groups)
    seen: set[tuple[str, str]] = set()

    for left_group, right_group in resolved:
        if left_group not in valid_groups:
            msg = f"{context} contains Unknown comparison group: {left_group}"
            raise InputCompatibilityError(msg)
        if right_group not in valid_groups:
            msg = f"{context} contains Unknown comparison group: {right_group}"
            raise InputCompatibilityError(msg)
        if left_group == right_group:
            msg = (
                f"{context} contains Self comparison pair: "
                f"{left_group!r}, {right_group!r}"
            )
            raise InputCompatibilityError(msg)

        canonical_pair = tuple(sorted((left_group, right_group)))
        if canonical_pair in seen:
            msg = (
                f"{context} contains Duplicate comparison pair regardless of "
                f"direction: {left_group!r}, {right_group!r}"
            )
            raise InputCompatibilityError(msg)
        seen.add(canonical_pair)

    return resolved


def validate_dataset_comparisons(
    *,
    schema: DatasetSchema,
    comparisons: Sequence[ComparisonSpec] | None,
    context: str,
) -> tuple[ComparisonSpec, ...] | None:
    """Validate dataset comparison definitions against the dataset schema."""

    try:
        return validate_comparison_specs(
            comparison_groups=schema.comparison_groups,
            comparisons=comparisons,
            context=context,
        )
    except (InputCompatibilityError, TypeError, ValueError):
        raise


__all__ = ["validate_comparison_specs", "validate_dataset_comparisons"]
