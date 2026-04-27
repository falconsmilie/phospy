"""Public dataset preprocessing configuration models."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from phospy.errors.input import PhosPyInputError

DATASET_MISSING_DATA_POLICY_FORBID = "forbid"
DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN = "impute_row_median"
DatasetMissingDataPolicy = Literal[
    "forbid",
    "impute_row_median",
]
DATASET_MISSING_DATA_POLICIES = frozenset(
    {
        DATASET_MISSING_DATA_POLICY_FORBID,
        DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
    }
)
DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY = "identity"
DATASET_INTENSITY_TRANSFORM_POLICY_LOG2 = "log2"
DatasetIntensityTransformPolicy = Literal["identity", "log2"]
DATASET_INTENSITY_TRANSFORM_POLICIES = frozenset(
    {
        DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY,
        DATASET_INTENSITY_TRANSFORM_POLICY_LOG2,
    }
)
DATASET_NORMALISATION_POLICY_NONE = "none"
DATASET_NORMALISATION_POLICY_MEDIAN_CENTER = "median_center"
DATASET_NORMALISATION_POLICY_QUANTILE = "quantile"
DatasetNormalisationPolicy = Literal["none", "median_center", "quantile"]
DATASET_NORMALISATION_POLICIES = frozenset(
    {
        DATASET_NORMALISATION_POLICY_NONE,
        DATASET_NORMALISATION_POLICY_MEDIAN_CENTER,
        DATASET_NORMALISATION_POLICY_QUANTILE,
    }
)
DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE = "none"
DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL = "subtract_log_total"
DatasetTotalProteinCorrectionPolicy = Literal[
    "none",
    "subtract_log_total",
]
DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES = frozenset(
    {
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
    }
)
DATASET_SITE_MATRIX_POLICY_AS_INPUT = "as_input"
DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA = "build_from_metadata"
DatasetSiteMatrixPolicy = Literal["as_input", "build_from_metadata"]
DATASET_SITE_MATRIX_POLICIES = frozenset(
    {
        DATASET_SITE_MATRIX_POLICY_AS_INPUT,
        DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA,
    }
)
DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL = "max_mean_signal"
DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST = "first"
DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN = "aggregate_mean"
DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN = "aggregate_median"
DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR = "error"
DatasetSiteMatrixDuplicateSitePolicy = Literal[
    "max_mean_signal",
    "first",
    "aggregate_mean",
    "aggregate_median",
    "error",
]
DATASET_SITE_MATRIX_DUPLICATE_POLICIES = frozenset(
    {
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL,
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST,
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN,
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN,
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR,
    }
)
DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING = "drop_any_missing"
DatasetSiteMatrixMissingDataPolicy = Literal["drop_any_missing",]
DATASET_SITE_MATRIX_MISSING_DATA_POLICIES = frozenset(
    {
        DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    }
)
DATASET_COMPARISON_BUILDING_POLICY_NONE = "none"
DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS = "sample_metadata_pairs"
DatasetComparisonBuildingPolicy = Literal["none", "sample_metadata_pairs"]
DatasetComparisonPair = tuple[str, str]
DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN = "comparison_group"
DATASET_COMPARISON_BUILDING_POLICIES = frozenset(
    {
        DATASET_COMPARISON_BUILDING_POLICY_NONE,
        DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS,
    }
)

_INCOMPATIBLE_SITE_MATRIX_MISSING_DATA_POLICIES = frozenset(
    {"retain_missing", "require_min_observed_values"}
)
_SUPPORTED_SITE_MATRIX_MISSING_DATA_POLICY = (
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
)


@dataclass(frozen=True, slots=True)
class DatasetMissingDataConfig:
    """Public missing-data policy options for dataset building.

    - `"forbid"`: do not preprocess missing values (strict default behavior).
    - `"impute_row_median"`: for each site row, drop rows with fewer than
      `min_observed_values` quantified samples, then impute remaining missing
      values with that row's observed-value median.

    `min_observed_values` is required for `"impute_row_median"` and must stay
    unset for `"forbid"`.
    """

    policy: DatasetMissingDataPolicy = DATASET_MISSING_DATA_POLICY_FORBID
    min_observed_values: int | None = None

    def __post_init__(self) -> None:
        policy = self.policy
        if policy not in DATASET_MISSING_DATA_POLICIES:
            supported = ", ".join(sorted(DATASET_MISSING_DATA_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.policy "
                f"must be one of: {supported}"
            )

        min_observed_values = self.min_observed_values
        if policy == DATASET_MISSING_DATA_POLICY_FORBID:
            if min_observed_values is not None:
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.min_observed_values must be "
                    "None when missing_data.policy='forbid'"
                )
            return
        if policy == DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN:
            if isinstance(min_observed_values, bool) or not isinstance(
                min_observed_values, int
            ):
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.min_observed_values must be an "
                    "int when missing_data.policy='impute_row_median'"
                )
            if min_observed_values < 1:
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.min_observed_values "
                    "must be greater than or equal to 1 when "
                    "missing_data.policy='impute_row_median'"
                )
            return
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "missing_data.policy"
        )


@dataclass(frozen=True, slots=True)
class DatasetIntensityTransformConfig:
    """Public intensity transform policy options for dataset building.

    - `"identity"`: no transform (strict default).
    - `"log2"`: apply `log2(value + pseudocount)` to quantitative matrix values.

    `pseudocount` must be non-negative.
    """

    policy: DatasetIntensityTransformPolicy = (
        DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY
    )
    pseudocount: float = 1.0

    def __post_init__(self) -> None:
        policy = self.policy
        if policy not in DATASET_INTENSITY_TRANSFORM_POLICIES:
            supported = ", ".join(sorted(DATASET_INTENSITY_TRANSFORM_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                f"policy must be one of: {supported}"
            )

        pseudocount = self.pseudocount
        if isinstance(pseudocount, bool) or not isinstance(pseudocount, (int, float)):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                "pseudocount must be a float or int"
            )
        if not math.isfinite(float(pseudocount)):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                "pseudocount must be finite"
            )
        if pseudocount < 0:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform."
                "pseudocount must be greater than or equal to 0"
            )


@dataclass(frozen=True, slots=True)
class DatasetNormalisationConfig:
    """Public normalisation policy options for dataset building.

    - `"none"`: no normalisation (strict default).
    - `"median_center"`: subtract sample-wise medians.
    - `"quantile"`: force sample columns to share one empirical distribution.
    """

    policy: DatasetNormalisationPolicy = DATASET_NORMALISATION_POLICY_NONE

    def __post_init__(self) -> None:
        policy = self.policy
        if policy not in DATASET_NORMALISATION_POLICIES:
            supported = ", ".join(sorted(DATASET_NORMALISATION_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.normalisation.policy "
                f"must be one of: {supported}"
            )


@dataclass(frozen=True, slots=True)
class DatasetTotalProteinCorrectionConfig:
    """Public total/protein correction policy options for dataset building.

    - `"none"`: do not apply total/protein correction.
    - `"subtract_log_total"`: subtract matched log-scale total-protein abundance
      from log-scale phosphosite abundance in the builder preprocessing lane:
      `log2_phospho - log2_total`.
    """

    policy: DatasetTotalProteinCorrectionPolicy = (
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )

    def __post_init__(self) -> None:
        policy = self.policy
        if policy not in DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES:
            supported = ", ".join(sorted(DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"policy must be one of: {supported}"
            )


@dataclass(frozen=True, slots=True)
class DatasetSiteMatrixConfig:
    """Public site-matrix policy options for dataset building.

    The supported public builder story is intentionally narrow: PhosPy builds a
    strict, missing-value-free `AnalysisReadyPhosphoDataset`.

    `policy` controls whether the site-matrix stage runs:

    - `"as_input"`: preserve interpreted site-matrix-ready rows as provided.
    - `"build_from_metadata"`: construct site-matrix-ready rows from
      `site_metadata` (`gene_symbol`, `site`) after upstream
      missing-data/total-correction preprocessing, using row-level
      `site_sequence` support from supplied values and/or supported derivation.

    When `policy="build_from_metadata"`, the supported public row-retention mode is
    fixed to `missing_data_policy="drop_any_missing"`, so only complete phospho
    rows enter strict `AnalysisReadyPhosphoDataset` construction.

    `duplicate_site_policy` controls duplicate-site collapse for the retained
    complete-case rows:

    - `"error"`: strict/scientifically cautious mode; fail when duplicate
      constructed site identifiers are present.
    - `"first"`: convenient input-order rule; keep the first encountered row and
      drop later duplicates.
    - `"max_mean_signal"` (default): keep the highest-mean signal row among
      duplicates, which can favour higher-abundance / stronger-signal rows.
    - `"aggregate_mean"`: aggregate duplicate phospho values by column mean,
      preserving all rows numerically but potentially blurring distinct peptide
      contexts.
    - `"aggregate_median"`: aggregate duplicate phospho values by column median,
      preserving all rows numerically with similar context-blurring trade-offs.

    Duplicate-row resolution and metadata conflicts are reported through
    `dataset.preprocessing_report` when this stage executes.

    `minimum_observed_values` remains internal-only state and must
    stay unset in the supported public builder lane.
    """

    policy: DatasetSiteMatrixPolicy = DATASET_SITE_MATRIX_POLICY_AS_INPUT
    duplicate_site_policy: DatasetSiteMatrixDuplicateSitePolicy = (
        DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL
    )
    missing_data_policy: DatasetSiteMatrixMissingDataPolicy = (
        DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
    )
    minimum_observed_values: int | None = None

    def __post_init__(self) -> None:
        policy = self.policy
        if policy not in DATASET_SITE_MATRIX_POLICIES:
            supported = ", ".join(sorted(DATASET_SITE_MATRIX_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix.policy "
                f"must be one of: {supported}"
            )

        duplicate_site_policy = self.duplicate_site_policy
        if duplicate_site_policy not in DATASET_SITE_MATRIX_DUPLICATE_POLICIES:
            supported_duplicates = ", ".join(
                sorted(DATASET_SITE_MATRIX_DUPLICATE_POLICIES)
            )
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "duplicate_site_policy must be one of: "
                f"{supported_duplicates}"
            )

        missing_data_policy = self.missing_data_policy
        if missing_data_policy not in DATASET_SITE_MATRIX_MISSING_DATA_POLICIES:
            if missing_data_policy in _INCOMPATIBLE_SITE_MATRIX_MISSING_DATA_POLICIES:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.site_matrix."
                    f"missing_data_policy='{missing_data_policy}' is not supported "
                    "for strict AnalysisReadyPhosphoDataset construction in the "
                    "public complete-case builder lane. Use "
                    "site_matrix.missing_data_policy="
                    f"'{_SUPPORTED_SITE_MATRIX_MISSING_DATA_POLICY}'."
                )
            supported_missing_policies = ", ".join(
                sorted(DATASET_SITE_MATRIX_MISSING_DATA_POLICIES)
            )
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "missing_data_policy must be one of: "
                f"{supported_missing_policies}"
            )

        minimum_observed_values = self.minimum_observed_values
        if minimum_observed_values is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "minimum_observed_values is not supported for strict "
                "AnalysisReadyPhosphoDataset construction and must be None"
            )

        if policy == DATASET_SITE_MATRIX_POLICY_AS_INPUT and (
            duplicate_site_policy
            != DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "duplicate_site_policy is only valid when "
                "site_matrix.policy='build_from_metadata'"
            )


@dataclass(frozen=True, slots=True)
class DatasetComparisonBuildingConfig:
    """Public comparison-building policy options for dataset building.

    - `"none"`: do not build dataset-level pairwise comparisons.
    - `"sample_metadata_pairs"`: build comparison columns from grouped sample
      metadata.

    For `"sample_metadata_pairs"`:

    - `sample_group_column` must exist in `sample_metadata` and define one
      non-empty group label per sample.
    - `pairs` supports explicit pass-through comparisons as `(left, right)`
      tuples. If omitted, comparisons are inferred from all observed groups.
    """

    policy: DatasetComparisonBuildingPolicy = DATASET_COMPARISON_BUILDING_POLICY_NONE
    sample_group_column: str = DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN
    pairs: tuple[DatasetComparisonPair, ...] | None = None

    def __post_init__(self) -> None:
        policy = self.policy
        if policy not in DATASET_COMPARISON_BUILDING_POLICIES:
            supported = ", ".join(sorted(DATASET_COMPARISON_BUILDING_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.policy "
                f"must be one of: {supported}"
            )
        sample_group_column = self.sample_group_column
        if not isinstance(sample_group_column, str) or not sample_group_column.strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons."
                "sample_group_column must be a non-empty string"
            )
        pairs = self.pairs
        if policy == DATASET_COMPARISON_BUILDING_POLICY_NONE:
            if pairs is not None:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "must be None when comparisons.policy='none'"
                )
            return
        if policy != DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS:
            raise PhosPyInputError(
                "dataset build request preprocessing_config contains an unsupported "
                "comparisons.policy"
            )
        if pairs is None:
            return
        if not isinstance(pairs, (tuple, list)):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs must be "
                "a sequence of (left_group, right_group) pairs when provided"
            )
        resolved_pairs = tuple(pairs)
        if not resolved_pairs:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.pairs must "
                "contain at least one pair when provided"
            )
        seen_pairs: set[tuple[str, str]] = set()
        for pair in resolved_pairs:
            if not isinstance(pair, tuple) or len(pair) != 2:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "must contain only (left_group, right_group) tuples"
                )
            left_group, right_group = pair
            if not isinstance(left_group, str) or not left_group.strip():
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "must contain non-empty left_group strings"
                )
            if not isinstance(right_group, str) or not right_group.strip():
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "must contain non-empty right_group strings"
                )
            left = left_group.strip()
            right = right_group.strip()
            if left == right:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "cannot contain self-comparison pairs"
                )
            canonical_pair = tuple(sorted((left, right)))
            if canonical_pair in seen_pairs:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.comparisons.pairs "
                    "contains duplicate pairs regardless of direction"
                )
            seen_pairs.add(canonical_pair)


__all__ = [
    "DATASET_COMPARISON_BUILDING_POLICIES",
    "DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN",
    "DATASET_COMPARISON_BUILDING_POLICY_NONE",
    "DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS",
    "DATASET_INTENSITY_TRANSFORM_POLICIES",
    "DATASET_INTENSITY_TRANSFORM_POLICY_IDENTITY",
    "DATASET_INTENSITY_TRANSFORM_POLICY_LOG2",
    "DATASET_MISSING_DATA_POLICIES",
    "DATASET_MISSING_DATA_POLICY_FORBID",
    "DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN",
    "DATASET_NORMALISATION_POLICIES",
    "DATASET_NORMALISATION_POLICY_MEDIAN_CENTER",
    "DATASET_NORMALISATION_POLICY_NONE",
    "DATASET_NORMALISATION_POLICY_QUANTILE",
    "DATASET_SITE_MATRIX_POLICIES",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICIES",
    "DATASET_SITE_MATRIX_MISSING_DATA_POLICIES",
    "DATASET_SITE_MATRIX_POLICY_AS_INPUT",
    "DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR",
    "DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL",
    "DatasetComparisonBuildingConfig",
    "DatasetComparisonPair",
    "DatasetComparisonBuildingPolicy",
    "DatasetIntensityTransformConfig",
    "DatasetIntensityTransformPolicy",
    "DatasetMissingDataConfig",
    "DatasetMissingDataPolicy",
    "DatasetNormalisationConfig",
    "DatasetNormalisationPolicy",
    "DatasetSiteMatrixConfig",
    "DatasetSiteMatrixDuplicateSitePolicy",
    "DatasetSiteMatrixMissingDataPolicy",
    "DatasetSiteMatrixPolicy",
    "DatasetTotalProteinCorrectionConfig",
    "DatasetTotalProteinCorrectionPolicy",
]
