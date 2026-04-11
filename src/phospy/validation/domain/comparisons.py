from __future__ import annotations

from collections.abc import Sequence

from ...datasets.schema import DatasetSchema
from ...errors import InputCompatibilityError
from ...internal.constants import ComparisonSpec


def validate_dataset_comparisons(
    *,
    schema: DatasetSchema,
    comparisons: Sequence[ComparisonSpec] | None,
    context: str,
) -> tuple[ComparisonSpec, ...] | None:
    """Validate dataset comparison definitions against the dataset schema."""

    try:
        return schema.validate_comparisons(comparisons, context=context)
    except (InputCompatibilityError, TypeError, ValueError):
        raise


__all__ = ["validate_dataset_comparisons"]
