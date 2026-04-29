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
    ComparisonBuildResult,
    PreprocessingStageResult,
    PreprocessingState,
)
from phospy.datasets.preprocessing.report_rows import (
    report_rows_from_comparison_group_stats_dataframe,
    report_rows_from_comparison_pair_stats_dataframe,
)
from phospy.datasets.preprocessing.report_schema import (
    COMPARISON_GROUP_STATS_COLUMNS,
    COMPARISON_PAIR_STATS_COLUMNS,
)
from phospy.errors.input import PhosPyInputError
from phospy.provenance.hashing import hash_table

_COMPARISON_OUTPUT_PREFIX = "p_"


class ComparisonsStage:
    """Build dataset-level pairwise comparison columns when requested."""

    stage_key = DATASET_PREPROCESSING_STAGE_COMPARISONS

    def run(self, state: PreprocessingState) -> PreprocessingStageResult:
        policy = state.plan.comparison_building_policy
        if policy == DATASET_COMPARISON_BUILDING_POLICY_NONE:
            return PreprocessingStageResult(
                state=state,
                diagnostics={
                    "dropped_row_ids": (),
                    "dropped_row_count": 0,
                    "imputed_cell_count": 0,
                    "imputed_row_ids": (),
                    "notes": "stage executed",
                    "diagnostics": {
                        "policy": policy,
                        "sample_group_column": state.plan.comparison_sample_group_column,
                        "resolved_comparison_pairs": [],
                        "group_labels": [],
                        "output_comparison_hash": None,
                    },
                },
            )
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
        build_result = _build_comparison_build_result(
            phospho=state.phospho,
            group_to_samples=group_to_samples,
            pairs=pairs,
        )
        next_state = replace(
            state,
            comparisons=build_result.comparisons,
            comparison_group_stats=build_result.comparison_group_stats,
            comparison_pair_stats=build_result.comparison_pair_stats,
        )
        return PreprocessingStageResult(
            state=next_state,
            report_rows=(
                report_rows_from_comparison_group_stats_dataframe(
                    build_result.comparison_group_stats
                )
                + report_rows_from_comparison_pair_stats_dataframe(
                    build_result.comparison_pair_stats
                )
            ),
            diagnostics={
                "dropped_row_ids": (),
                "dropped_row_count": 0,
                "imputed_cell_count": 0,
                "imputed_row_ids": (),
                "notes": "stage executed",
                "diagnostics": {
                    "policy": policy,
                    "sample_group_column": state.plan.comparison_sample_group_column,
                    "resolved_comparison_pairs": _resolve_comparison_pairs(next_state),
                    "group_labels": _resolve_group_labels_from_stats(next_state),
                    "output_comparison_hash": hash_table(
                        build_result.comparisons,
                        name="comparisons.output.table",
                    ),
                },
            },
        )


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
        canonical_pair = (left, right) if left <= right else (right, left)
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


def _build_group_statistics(
    *,
    phospho: pd.DataFrame,
    group_to_samples: dict[str, list[str]],
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    group_records: list[pd.DataFrame] = []
    group_summaries: dict[str, pd.DataFrame] = {}
    site_ids = phospho.index.astype(str).tolist()
    for group_name, sample_columns in group_to_samples.items():
        group_values = phospho.loc[:, sample_columns]
        n = group_values.count(axis=1).astype(int)
        mean = group_values.mean(axis=1, skipna=True)
        sd = group_values.std(axis=1, skipna=True)
        sem = sd / n.astype("float64").pow(0.5).where(n > 0)
        median = group_values.median(axis=1, skipna=True)
        minimum = group_values.min(axis=1, skipna=True)
        maximum = group_values.max(axis=1, skipna=True)
        summary = pd.DataFrame(
            {
                "n": n,
                "mean": mean,
                "sd": sd,
                "sem": sem,
                "median": median,
                "min": minimum,
                "max": maximum,
            },
            index=phospho.index.copy(),
        )
        group_summaries[group_name] = summary
        group_records.append(
            pd.DataFrame(
                {
                    "site_id": site_ids,
                    "group": [group_name] * len(site_ids),
                    "n": n.to_numpy(),
                    "mean": mean.to_numpy(),
                    "sd": sd.to_numpy(),
                    "sem": sem.to_numpy(),
                    "median": median.to_numpy(),
                    "min": minimum.to_numpy(),
                    "max": maximum.to_numpy(),
                    "sample_ids": [tuple(sample_columns)] * len(site_ids),
                },
                columns=pd.Index(COMPARISON_GROUP_STATS_COLUMNS),
            )
        )
    if not group_records:
        return (
            pd.DataFrame.from_records([], columns=COMPARISON_GROUP_STATS_COLUMNS),
            group_summaries,
        )
    group_stats = pd.concat(group_records, axis=0, ignore_index=True)
    return group_stats, group_summaries


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


def _build_pair_statistics(
    *,
    phospho: pd.DataFrame,
    pairs: tuple[DatasetComparisonPair, ...],
    comparison_matrix: pd.DataFrame,
    group_summaries: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    pair_records: list[pd.DataFrame] = []
    site_ids = phospho.index.astype(str).tolist()
    for left_group, right_group in pairs:
        left_summary = group_summaries[left_group]
        right_summary = group_summaries[right_group]
        comparison_name = f"{_COMPARISON_OUTPUT_PREFIX}{left_group}_{right_group}"
        effect_size = comparison_matrix.loc[:, comparison_name]
        pair_records.append(
            pd.DataFrame(
                {
                    "site_id": site_ids,
                    "comparison": [comparison_name] * len(site_ids),
                    "left_group": [left_group] * len(site_ids),
                    "right_group": [right_group] * len(site_ids),
                    "left_n": left_summary.loc[:, "n"].to_numpy(),
                    "right_n": right_summary.loc[:, "n"].to_numpy(),
                    "left_mean": left_summary.loc[:, "mean"].to_numpy(),
                    "right_mean": right_summary.loc[:, "mean"].to_numpy(),
                    "left_sd": left_summary.loc[:, "sd"].to_numpy(),
                    "right_sd": right_summary.loc[:, "sd"].to_numpy(),
                    "left_sem": left_summary.loc[:, "sem"].to_numpy(),
                    "right_sem": right_summary.loc[:, "sem"].to_numpy(),
                    "effect_size": effect_size.to_numpy(),
                    "left_median": left_summary.loc[:, "median"].to_numpy(),
                    "right_median": right_summary.loc[:, "median"].to_numpy(),
                    "left_min": left_summary.loc[:, "min"].to_numpy(),
                    "right_min": right_summary.loc[:, "min"].to_numpy(),
                    "left_max": left_summary.loc[:, "max"].to_numpy(),
                    "right_max": right_summary.loc[:, "max"].to_numpy(),
                },
                columns=pd.Index(COMPARISON_PAIR_STATS_COLUMNS),
            )
        )
    if not pair_records:
        return pd.DataFrame.from_records([], columns=COMPARISON_PAIR_STATS_COLUMNS)
    return pd.concat(pair_records, axis=0, ignore_index=True)


def _build_comparison_build_result(
    *,
    phospho: pd.DataFrame,
    group_to_samples: dict[str, list[str]],
    pairs: tuple[DatasetComparisonPair, ...],
) -> ComparisonBuildResult:
    group_stats, group_summaries = _build_group_statistics(
        phospho=phospho,
        group_to_samples=group_to_samples,
    )
    group_means = _build_group_means(
        phospho=phospho,
        group_to_samples=group_to_samples,
    )
    comparisons = _build_comparison_matrix(group_means=group_means, pairs=pairs)
    pair_stats = _build_pair_statistics(
        phospho=phospho,
        pairs=pairs,
        comparison_matrix=comparisons,
        group_summaries=group_summaries,
    )
    return ComparisonBuildResult(
        comparisons=comparisons,
        comparison_group_stats=group_stats,
        comparison_pair_stats=pair_stats,
    )


def _resolve_comparison_pairs(state: PreprocessingState) -> list[tuple[str, str]]:
    pair_stats = state.comparison_pair_stats
    if pair_stats is None or pair_stats.empty:
        return []
    pairs = pair_stats.loc[:, ["left_group", "right_group"]].drop_duplicates()
    return [
        (str(left), str(right))
        for left, right in pairs.itertuples(index=False, name=None)
    ]


def _resolve_group_labels_from_stats(state: PreprocessingState) -> list[str]:
    group_stats = state.comparison_group_stats
    if group_stats is None or group_stats.empty:
        return []
    labels = group_stats.loc[:, "group"].astype(str).drop_duplicates().tolist()
    return [str(label) for label in labels]


__all__ = ["ComparisonsStage"]
