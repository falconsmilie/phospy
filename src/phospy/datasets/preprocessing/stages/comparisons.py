"""Comparison-building stage for dataset preprocessing."""

from __future__ import annotations

from dataclasses import replace
from itertools import combinations

import pandas as pd

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS,
    DatasetComparisonPair,
)
from phospy.datasets.preprocessing.models import (
    DATASET_PREPROCESSING_STAGE_COMPARISONS,
    PreprocessingState,
)
from phospy.errors.input import PhosPyInputError

_COMPARISON_OUTPUT_PREFIX = "p_"


class ComparisonsStage:
    """Build dataset-level pairwise comparison columns when requested."""

    stage_key = DATASET_PREPROCESSING_STAGE_COMPARISONS

    def run(self, state: PreprocessingState) -> PreprocessingState:
        policy = state.plan.comparison_building_policy
        if policy == DATASET_COMPARISON_BUILDING_POLICY_NONE:
            return state
        if policy != DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS:
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "comparisons.policy"
            )

        sample_metadata = state.sample_metadata
        if sample_metadata is None:
            raise PhosPyInputError(
                "dataset build request preprocessing comparison-building "
                "policy='sample_metadata_pairs' requires sample_metadata input data"
            )
        if not sample_metadata.index.equals(state.phospho.columns):
            raise PhosPyInputError(
                "dataset build request preprocessing comparison-building requires "
                "sample_metadata.index to exactly match phospho.columns"
            )

        group_column = state.plan.comparison_sample_group_column
        if group_column not in sample_metadata.columns:
            raise PhosPyInputError(
                "dataset build request preprocessing comparison-building requires "
                f"sample_metadata column '{group_column}'"
            )
        group_labels = _resolve_group_labels(
            sample_metadata=sample_metadata,
            group_column=group_column,
        )
        group_to_samples = _resolve_group_to_samples(
            sample_index=sample_metadata.index,
            group_labels=group_labels,
        )
        pairs = _resolve_pairs(
            groups=tuple(group_to_samples.keys()),
            explicit_pairs=state.plan.comparison_pairs,
        )
        group_means = _build_group_means(
            phospho=state.phospho,
            group_to_samples=group_to_samples,
        )
        comparison_matrix = _build_comparison_matrix(
            group_means=group_means,
            pairs=pairs,
        )
        return replace(state, comparisons=comparison_matrix)


def _resolve_group_labels(
    *,
    sample_metadata: pd.DataFrame,
    group_column: str,
) -> pd.Series:
    raw_values = sample_metadata.loc[:, group_column]
    normalized_values = raw_values.astype("string").str.strip()
    invalid_mask = (
        raw_values.isna() | normalized_values.isna() | (normalized_values == "")
    )
    if bool(invalid_mask.any()):
        invalid_samples = sample_metadata.index[invalid_mask].astype(str).tolist()
        preview = ", ".join(invalid_samples[:5])
        suffix = "" if len(invalid_samples) <= 5 else ", ..."
        raise PhosPyInputError(
            "dataset build request preprocessing comparison-building requires "
            f"sample_metadata.{group_column} to contain non-empty group labels "
            f"for all samples; invalid sample IDs: {preview}{suffix}"
        )
    return normalized_values.astype(str)


def _resolve_group_to_samples(
    *,
    sample_index: pd.Index,
    group_labels: pd.Series,
) -> dict[str, list[str]]:
    group_to_samples: dict[str, list[str]] = {}
    for sample_id, group_label in zip(
        sample_index.astype(str), group_labels, strict=True
    ):
        group_to_samples.setdefault(str(group_label), []).append(str(sample_id))
    return group_to_samples


def _resolve_pairs(
    *,
    groups: tuple[str, ...],
    explicit_pairs: tuple[DatasetComparisonPair, ...] | None,
) -> tuple[DatasetComparisonPair, ...]:
    if explicit_pairs is None:
        if len(groups) < 2:
            raise PhosPyInputError(
                "dataset build request preprocessing comparison-building requires "
                "at least two sample groups when inferring comparisons"
            )
        return tuple(combinations(sorted(groups), 2))

    resolved_pairs = tuple(explicit_pairs)
    if not resolved_pairs:
        raise PhosPyInputError(
            "dataset build request preprocessing comparison-building "
            "comparisons.pairs must contain at least one pair when provided"
        )
    valid_groups = set(groups)
    seen_pairs: set[tuple[str, str]] = set()
    for left_group, right_group in resolved_pairs:
        left = str(left_group).strip()
        right = str(right_group).strip()
        if not left or not right:
            raise PhosPyInputError(
                "dataset build request preprocessing comparison-building "
                "comparisons.pairs must contain non-empty group names"
            )
        if left == right:
            raise PhosPyInputError(
                "dataset build request preprocessing comparison-building "
                "comparisons.pairs cannot contain self-comparison pairs"
            )
        canonical_pair = tuple(sorted((left, right)))
        if canonical_pair in seen_pairs:
            raise PhosPyInputError(
                "dataset build request preprocessing comparison-building "
                "comparisons.pairs contains duplicate pairs regardless of direction"
            )
        seen_pairs.add(canonical_pair)
        unknown_groups = [group for group in (left, right) if group not in valid_groups]
        if unknown_groups:
            unknown_list = ", ".join(unknown_groups)
            supported = ", ".join(sorted(valid_groups))
            raise PhosPyInputError(
                "dataset build request preprocessing comparison-building "
                "comparisons.pairs references unknown sample groups: "
                f"{unknown_list}. Supported groups: {supported}"
            )
    return tuple(
        (str(left).strip(), str(right).strip()) for left, right in resolved_pairs
    )


def _build_group_means(
    *,
    phospho: pd.DataFrame,
    group_to_samples: dict[str, list[str]],
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            group_name: phospho.loc[:, sample_columns].mean(axis=1, skipna=True)
            for group_name, sample_columns in group_to_samples.items()
        },
        index=phospho.index.copy(),
    )


def _build_comparison_matrix(
    *,
    group_means: pd.DataFrame,
    pairs: tuple[DatasetComparisonPair, ...],
) -> pd.DataFrame:
    comparison_matrix = pd.DataFrame(index=group_means.index.copy())
    for left_group, right_group in pairs:
        comparison_matrix[f"{_COMPARISON_OUTPUT_PREFIX}{left_group}_{right_group}"] = (
            group_means.loc[:, left_group] - group_means.loc[:, right_group]
        )
    return comparison_matrix


__all__ = ["ComparisonsStage"]
