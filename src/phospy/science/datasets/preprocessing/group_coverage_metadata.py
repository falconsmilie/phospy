"""Coverage-filter validation built on shared sample-group resolution."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.errors.input import PhosPyInputError
from phospy.science.datasets.preprocessing.sample_group_metadata import (
    ResolvedSampleGroups,
    SampleGroupMetadataResolver,
)


@dataclass(frozen=True, slots=True)
class ResolvedGroupCoverageFilterMetadata(ResolvedSampleGroups):
    """Compatibility result type for group coverage consumers."""


class GroupCoverageFilterMetadataValidator:
    """Validate and resolve sample groups for coverage filtering."""

    def __init__(
        self,
        *,
        sample_group_resolver: SampleGroupMetadataResolver | None = None,
    ) -> None:
        self._sample_group_resolver = (
            sample_group_resolver or SampleGroupMetadataResolver()
        )

    def run(
        self,
        *,
        phospho: pd.DataFrame,
        sample_metadata: pd.DataFrame | None,
        group_column: str | None,
        min_groups_passing_threshold: int,
    ) -> ResolvedGroupCoverageFilterMetadata:
        resolved = self._sample_group_resolver.run(
            phospho=phospho,
            sample_metadata=sample_metadata,
            group_column=group_column,
        )
        if min_groups_passing_threshold > len(resolved.sample_order_by_group):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.group_coverage_filter."
                "min_groups_passing_threshold cannot exceed the number of observed "
                f"groups; threshold={int(min_groups_passing_threshold)}, "
                f"observed_groups={int(len(resolved.sample_order_by_group))}"
            )
        return ResolvedGroupCoverageFilterMetadata(
            group_column=resolved.group_column,
            original_sample_order=resolved.original_sample_order,
            sample_order=resolved.sample_order,
            canonical_label_by_original=resolved.canonical_label_by_original,
            original_label_by_canonical=resolved.original_label_by_canonical,
            group_by_sample=resolved.group_by_sample,
            group_by_original_sample=resolved.group_by_original_sample,
            group_by_column_position=resolved.group_by_column_position,
            sample_order_by_group=resolved.sample_order_by_group,
            original_sample_order_by_group=(resolved.original_sample_order_by_group),
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


__all__ = [
    "GroupCoverageFilterMetadataValidator",
    "ResolvedGroupCoverageFilterMetadata",
    "require_numeric_group_coverage_matrix",
]
