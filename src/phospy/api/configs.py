"""Public workflow and dataset-preprocessing configuration models."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

from phospy.errors.input import PhosPyInputError
from phospy.errors.validation import WorkflowValidationError

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
DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_RATIO_TO_TOTAL = "ratio_to_total"
DatasetTotalProteinCorrectionPolicy = Literal[
    "none",
    "subtract_log_total",
    "ratio_to_total",
]
DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES = frozenset(
    {
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL,
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_RATIO_TO_TOTAL,
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
KINASE_SCORING_MIN_SUBSTRATES_FLOOR = 2
KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT = "strict"
KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA = "median_skipna"
KinaseProfileMissingValueStrategy = Literal[
    "strict",
    "median_skipna",
]
KINASE_PROFILE_MISSING_VALUE_STRATEGIES = frozenset(
    {
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT,
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA,
    }
)
KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR = 1
KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR = 1
KINASE_ACTIVITY_DEFAULT_THRESHOLD = 0.6
KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES = 3
KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES = 20
SIGNALOME_MODULE_COUNT_FLOOR = 1
SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT = 0.5
SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT = 0.1
SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR = 1
SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT = 10
SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY = "cutoff_binary"
SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP = "weighted_top"
SignalomeAssignmentPolicy = Literal["cutoff_binary", "weighted_top"]
SIGNALOME_ASSIGNMENT_POLICIES = frozenset(
    {
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
        SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    }
)
SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT = "allow_and_report"
SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP = "error_on_drop"
SignalomeScorePreconditioningPolicy = Literal["allow_and_report", "error_on_drop"]
SIGNALOME_SCORE_PRECONDITIONING_POLICIES = frozenset(
    {
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    }
)
SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY = "positive_only"
SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD = "absolute_threshold"
SIGNALOME_KINASE_NETWORK_POLICY_SIGNED = "signed"
SignalomeKinaseNetworkPolicy = Literal[
    "positive_only",
    "absolute_threshold",
    "signed",
]
SIGNALOME_KINASE_NETWORK_POLICIES = frozenset(
    {
        SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY,
        SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD,
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED,
    }
)
KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING = "deterministic_ranking"
KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE = "adaptive_ensemble"
KinasePredictionMode = Literal[
    "deterministic_ranking",
    "adaptive_ensemble",
]
KINASE_PREDICTION_MODES = frozenset(
    {
        KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
        KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    }
)
KINASE_ADAPTIVE_POLICY_STABLE = "stable"
KINASE_ADAPTIVE_POLICY_R_PARITY = "r_parity"
KinaseAdaptivePolicy = Literal["stable", "r_parity"]
KINASE_ADAPTIVE_POLICIES = frozenset(
    {
        KINASE_ADAPTIVE_POLICY_STABLE,
        KINASE_ADAPTIVE_POLICY_R_PARITY,
    }
)
KINASE_PREDICTION_DEFAULT_ITERATIONS = 5
_INCOMPATIBLE_SITE_MATRIX_MISSING_DATA_POLICIES = frozenset(
    {"retain_missing", "require_min_observed_values"}
)
_SUPPORTED_SITE_MATRIX_MISSING_DATA_POLICY = (
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
)


def resolve_dataset_total_protein_correction_policy(
    policy: DatasetTotalProteinCorrectionPolicy,
) -> DatasetTotalProteinCorrectionPolicy:
    """Resolve deprecated total-protein correction aliases to canonical policy."""

    if policy == DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_RATIO_TO_TOTAL:
        return DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_SUBTRACT_LOG_TOTAL
    return policy


def _require_int_at_least(
    value: object,
    *,
    field_name: str,
    minimum: int,
    error_type: type[Exception],
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise error_type(f"{field_name} must be an int")
    if value < minimum:
        raise error_type(f"{field_name} must be greater than or equal to {minimum}")
    return value


def _require_real_between(
    value: object,
    *,
    field_name: str,
    minimum: float,
    maximum: float,
    error_type: type[Exception],
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise error_type(
            f"{field_name} must be a float between {float(minimum):.1f} and "
            f"{float(maximum):.1f}"
        )
    numeric_value = float(value)
    if not minimum <= numeric_value <= maximum:
        raise error_type(
            f"{field_name} must be between {float(minimum):.1f} and "
            f"{float(maximum):.1f}"
        )
    return numeric_value


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
    - `"ratio_to_total"`: deprecated alias that resolves to
      `"subtract_log_total"`.
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

    `minimum_observed_values` remains internal-only compatibility state and must
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


@dataclass(frozen=True, slots=True)
class DatasetPreprocessingConfig:
    """Public grouped preprocessing policy for dataset building.

    The builder owns this policy surface. Groups are intentionally separated so
    supported preprocessing science remains user-visible:

    - `intensity_transform`: quantitative transform policy.
    - `normalisation`: sample-wise normalisation policy.
    - `missing_data`: missing-value handling policy.
    - `total_protein_correction`: total/protein correction policy.
    - `site_matrix`: site-matrix construction policy.
    - `comparisons`: comparison-building policy.
    """

    intensity_transform: DatasetIntensityTransformConfig = field(
        default_factory=DatasetIntensityTransformConfig
    )
    normalisation: DatasetNormalisationConfig = field(
        default_factory=DatasetNormalisationConfig
    )
    missing_data: DatasetMissingDataConfig = field(
        default_factory=DatasetMissingDataConfig
    )
    total_protein_correction: DatasetTotalProteinCorrectionConfig = field(
        default_factory=DatasetTotalProteinCorrectionConfig
    )
    site_matrix: DatasetSiteMatrixConfig = field(
        default_factory=DatasetSiteMatrixConfig
    )
    comparisons: DatasetComparisonBuildingConfig = field(
        default_factory=DatasetComparisonBuildingConfig
    )

    def __post_init__(self) -> None:
        if not isinstance(self.intensity_transform, DatasetIntensityTransformConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.intensity_transform "
                "must be a DatasetIntensityTransformConfig"
            )
        if not isinstance(self.normalisation, DatasetNormalisationConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.normalisation must be a "
                "DatasetNormalisationConfig"
            )
        if not isinstance(self.missing_data, DatasetMissingDataConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data must be a "
                "DatasetMissingDataConfig"
            )
        if not isinstance(
            self.total_protein_correction, DatasetTotalProteinCorrectionConfig
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction "
                "must be a DatasetTotalProteinCorrectionConfig"
            )
        if not isinstance(self.site_matrix, DatasetSiteMatrixConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix must be a "
                "DatasetSiteMatrixConfig"
            )
        if not isinstance(self.comparisons, DatasetComparisonBuildingConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons must be a "
                "DatasetComparisonBuildingConfig"
            )


@dataclass(frozen=True, slots=True)
class KinaseScoringConfig:
    """Public scoring-stage configuration.

    `min_substrates` is constrained to the public scoring support floor used by
    the supported rewrite contract.

    Supported scoring semantics are stage-pure: score generation is determined
    only by analysis-ready dataset values, resolved reference content, and this
    explicit scoring configuration. Prediction mode and reference input
    provenance (preset vs explicit bundle) do not redefine scoring behavior.

    `include_diagnostic_scoring_tables` controls publication of non-authoritative
    diagnostic scoring outputs (`motif_scores`, `score_fusion_weights`). The
    authoritative downstream lane (`rank_weighted_fusion_scores` with profile
    fallback) is always computed.

    `profile_missing_value_strategy` controls column-wise median behavior when a
    kinase profile is built from multiple quantified substrates:

    - `"strict"` propagates missing values (`median(..., skipna=False)`)
    - `"median_skipna"` ignores missing values (`median(..., skipna=True)`)
    """

    min_substrates: int = KINASE_SCORING_MIN_SUBSTRATES_FLOOR
    include_diagnostic_scoring_tables: bool = False
    profile_missing_value_strategy: KinaseProfileMissingValueStrategy = (
        KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT
    )

    def __post_init__(self) -> None:
        if not isinstance(self.include_diagnostic_scoring_tables, bool):
            raise WorkflowValidationError(
                "scoring_config.include_diagnostic_scoring_tables must be a bool"
            )
        if (
            self.profile_missing_value_strategy
            not in KINASE_PROFILE_MISSING_VALUE_STRATEGIES
        ):
            allowed = ", ".join(sorted(KINASE_PROFILE_MISSING_VALUE_STRATEGIES))
            raise WorkflowValidationError(
                "scoring_config.profile_missing_value_strategy must be one of: "
                f"{allowed}"
            )
        _require_int_at_least(
            self.min_substrates,
            field_name="scoring_config.min_substrates",
            minimum=KINASE_SCORING_MIN_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )


@dataclass(frozen=True, slots=True)
class KinasePredictionConfig:
    """Public prediction-stage configuration.

    `mode` selects the prediction lane:

    - `"deterministic_ranking"`: deterministic top-kinase selection from
      downstream scores.
    - `"adaptive_ensemble"`: real adaptive ensemble execution ported from donor
      science.

    `ensemble_size` is mode-dependent by design and should be interpreted with
    `mode`:

    - deterministic lane: maximum number of selected kinase columns.
    - adaptive lane: number of ensemble executions per kinase.
    """

    top_k: int = 30
    ensemble_size: int = 10
    mode: KinasePredictionMode = KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING
    adaptive_policy: KinaseAdaptivePolicy = KINASE_ADAPTIVE_POLICY_STABLE
    n_iterations: int = KINASE_PREDICTION_DEFAULT_ITERATIONS
    random_state: int | None = None

    def __post_init__(self) -> None:
        if self.mode not in KINASE_PREDICTION_MODES:
            allowed_modes = ", ".join(sorted(KINASE_PREDICTION_MODES))
            raise WorkflowValidationError(
                f"prediction_config.mode must be one of: {allowed_modes}"
            )
        if self.adaptive_policy not in KINASE_ADAPTIVE_POLICIES:
            allowed_policies = ", ".join(sorted(KINASE_ADAPTIVE_POLICIES))
            raise WorkflowValidationError(
                f"prediction_config.adaptive_policy must be one of: {allowed_policies}"
            )
        if self.random_state is not None:
            _require_int_at_least(
                self.random_state,
                field_name="prediction_config.random_state",
                minimum=0,
                error_type=WorkflowValidationError,
            )
        _require_int_at_least(
            self.top_k,
            field_name="prediction_config.top_k",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.ensemble_size,
            field_name="prediction_config.ensemble_size",
            minimum=1,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.n_iterations,
            field_name="prediction_config.n_iterations",
            minimum=1,
            error_type=WorkflowValidationError,
        )


@dataclass(frozen=True, slots=True)
class KinaseActivityConfig:
    """Configuration for the supported kinase activity stage.

    Activity runs inside `KinaseWorkflow` and can be disabled by setting either:

    - `activity_config=None` on `KinaseWorkflowRequest`, or
    - `activity_config.enabled=False`.
    """

    enabled: bool = True
    threshold: float = KINASE_ACTIVITY_DEFAULT_THRESHOLD
    min_substrates: int = KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES
    top_n_substrates: int = KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES

    def __post_init__(self) -> None:
        if not isinstance(self.enabled, bool):
            raise WorkflowValidationError("activity_config.enabled must be a bool")
        _require_real_between(
            self.threshold,
            field_name="activity_config.threshold",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.min_substrates,
            field_name="activity_config.min_substrates",
            minimum=KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.top_n_substrates,
            field_name="activity_config.top_n_substrates",
            minimum=KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR,
            error_type=WorkflowValidationError,
        )


@dataclass(frozen=True, slots=True)
class SignalomeConfig:
    """Public signalome workflow configuration.

    `network_policy` controls how score correlations are thresholded and encoded
    in `kinase_network.edges.correlation`:

    - `"positive_only"`: keep only positive correlations `>= threshold`.
    - `"absolute_threshold"`: keep correlations where `abs(correlation) >=
      threshold` and emit unsigned absolute correlation values.
    - `"signed"`: keep correlations where `abs(correlation) >= threshold` and
      emit signed correlation values.

    `assignment_policy` controls module-support attribution:

    - `"cutoff_binary"`: binary support per kinase from
      `substrate_support_cutoff`.
    - `"weighted_top"`: fractional support propagated from per-site
      `top_kinase_weights` ties.

    `score_preconditioning_policy` controls how all-missing downstream score
    rows are handled before score-driven signalome stages:

    - `"allow_and_report"` (default): drop all-missing rows and continue,
      reporting exact counts in diagnostics.
    - `"error_on_drop"`: fail signalome interpretation when any all-missing
      rows would be dropped.
    """

    substrate_support_cutoff: float = 0.5
    network_correlation_threshold: float = 0.5
    network_policy: SignalomeKinaseNetworkPolicy = (
        SIGNALOME_KINASE_NETWORK_POLICY_SIGNED
    )
    assignment_policy: SignalomeAssignmentPolicy = (
        SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    )
    score_preconditioning_policy: SignalomeScorePreconditioningPolicy = (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )
    module_count: int | None = None
    module_selection_primary_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT
    )
    module_selection_fallback_correlation_threshold: float = (
        SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT
    )
    module_selection_max_clusters: int = SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT

    def __post_init__(self) -> None:
        _require_real_between(
            self.substrate_support_cutoff,
            field_name="signalome workflow request config.substrate_support_cutoff",
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        _require_real_between(
            self.network_correlation_threshold,
            field_name=(
                "signalome workflow request config.network_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        if self.network_policy not in SIGNALOME_KINASE_NETWORK_POLICIES:
            allowed_policies = ", ".join(sorted(SIGNALOME_KINASE_NETWORK_POLICIES))
            raise WorkflowValidationError(
                "signalome workflow request config.network_policy "
                f"must be one of: {allowed_policies}"
            )
        if self.assignment_policy not in SIGNALOME_ASSIGNMENT_POLICIES:
            allowed_policies = ", ".join(sorted(SIGNALOME_ASSIGNMENT_POLICIES))
            raise WorkflowValidationError(
                "signalome workflow request config.assignment_policy "
                f"must be one of: {allowed_policies}"
            )
        if (
            self.score_preconditioning_policy
            not in SIGNALOME_SCORE_PRECONDITIONING_POLICIES
        ):
            allowed_policies = ", ".join(
                sorted(SIGNALOME_SCORE_PRECONDITIONING_POLICIES)
            )
            raise WorkflowValidationError(
                "signalome workflow request config.score_preconditioning_policy "
                f"must be one of: {allowed_policies}"
            )
        if self.module_count is not None:
            _require_int_at_least(
                self.module_count,
                field_name="signalome workflow request config.module_count",
                minimum=SIGNALOME_MODULE_COUNT_FLOOR,
                error_type=WorkflowValidationError,
            )
        _require_real_between(
            self.module_selection_primary_correlation_threshold,
            field_name=(
                "signalome workflow request config."
                "module_selection_primary_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        _require_real_between(
            self.module_selection_fallback_correlation_threshold,
            field_name=(
                "signalome workflow request config."
                "module_selection_fallback_correlation_threshold"
            ),
            minimum=0.0,
            maximum=1.0,
            error_type=WorkflowValidationError,
        )
        _require_int_at_least(
            self.module_selection_max_clusters,
            field_name="signalome workflow request config.module_selection_max_clusters",
            minimum=SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR,
            error_type=WorkflowValidationError,
        )


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
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_RATIO_TO_TOTAL",
    "resolve_dataset_total_protein_correction_policy",
    "DatasetComparisonBuildingConfig",
    "DatasetComparisonPair",
    "DatasetComparisonBuildingPolicy",
    "DatasetIntensityTransformConfig",
    "DatasetIntensityTransformPolicy",
    "DatasetMissingDataConfig",
    "DatasetNormalisationConfig",
    "DatasetNormalisationPolicy",
    "KINASE_ADAPTIVE_POLICIES",
    "KINASE_ADAPTIVE_POLICY_R_PARITY",
    "KINASE_ADAPTIVE_POLICY_STABLE",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGIES",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_MEDIAN_SKIPNA",
    "KINASE_PROFILE_MISSING_VALUE_STRATEGY_STRICT",
    "KINASE_PREDICTION_DEFAULT_ITERATIONS",
    "KINASE_PREDICTION_MODES",
    "KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE",
    "KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING",
    "KINASE_ACTIVITY_DEFAULT_MIN_SUBSTRATES",
    "KINASE_ACTIVITY_DEFAULT_THRESHOLD",
    "KINASE_ACTIVITY_DEFAULT_TOP_N_SUBSTRATES",
    "KINASE_ACTIVITY_MIN_SUBSTRATES_FLOOR",
    "KINASE_ACTIVITY_TOP_N_SUBSTRATES_FLOOR",
    "KINASE_SCORING_MIN_SUBSTRATES_FLOOR",
    "SIGNALOME_ASSIGNMENT_POLICIES",
    "SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY",
    "SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICIES",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT",
    "SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP",
    "SIGNALOME_KINASE_NETWORK_POLICIES",
    "SIGNALOME_KINASE_NETWORK_POLICY_ABSOLUTE_THRESHOLD",
    "SIGNALOME_KINASE_NETWORK_POLICY_POSITIVE_ONLY",
    "SIGNALOME_KINASE_NETWORK_POLICY_SIGNED",
    "SIGNALOME_MODULE_COUNT_FLOOR",
    "SIGNALOME_MODULE_SELECTION_FALLBACK_THRESHOLD_DEFAULT",
    "SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_DEFAULT",
    "SIGNALOME_MODULE_SELECTION_MAX_CLUSTERS_FLOOR",
    "SIGNALOME_MODULE_SELECTION_PRIMARY_THRESHOLD_DEFAULT",
    "DatasetMissingDataPolicy",
    "DatasetPreprocessingConfig",
    "DatasetSiteMatrixConfig",
    "DatasetSiteMatrixDuplicateSitePolicy",
    "DatasetSiteMatrixMissingDataPolicy",
    "DatasetSiteMatrixPolicy",
    "DatasetTotalProteinCorrectionConfig",
    "DatasetTotalProteinCorrectionPolicy",
    "KinaseAdaptivePolicy",
    "KinaseProfileMissingValueStrategy",
    "KinasePredictionMode",
    "KinaseActivityConfig",
    "KinasePredictionConfig",
    "KinaseScoringConfig",
    "SignalomeAssignmentPolicy",
    "SignalomeKinaseNetworkPolicy",
    "SignalomeScorePreconditioningPolicy",
    "SignalomeConfig",
]
