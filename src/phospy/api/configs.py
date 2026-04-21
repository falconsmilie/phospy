"""Public workflow and dataset-preprocessing configuration models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

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
DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE = "none"
DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_RATIO_TO_TOTAL = "ratio_to_total"
DatasetTotalProteinCorrectionPolicy = Literal["none", "ratio_to_total"]
DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES = frozenset(
    {
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE,
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
DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL = "max_mean_signal"
DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST = "first"
DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN = "aggregate_mean"
DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN = "aggregate_median"
DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_ERROR = "error"
DatasetSiteMatrixDuplicateSiteStrategy = Literal[
    "max_mean_signal",
    "first",
    "aggregate_mean",
    "aggregate_median",
    "error",
]
DATASET_SITE_MATRIX_DUPLICATE_STRATEGIES = frozenset(
    {
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL,
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST,
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN,
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN,
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_ERROR,
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


@dataclass(frozen=True, slots=True)
class DatasetTotalProteinCorrectionConfig:
    """Public total/protein correction policy options for dataset building.

    - `"none"`: do not apply total/protein correction.
    - `"ratio_to_total"`: subtract matched total/protein abundance from phosphosite
      abundance in the builder preprocessing lane.
    """

    policy: DatasetTotalProteinCorrectionPolicy = (
        DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE
    )


@dataclass(frozen=True, slots=True)
class DatasetSiteMatrixConfig:
    """Public site-matrix policy options for dataset building.

    `policy` controls whether the stage runs:

    - `"as_input"`: preserve interpreted site matrix as provided.
    - `"build_from_metadata"`: construct site-matrix-ready rows from
      `site_metadata` (`gene_symbol`, `site`) after upstream
      missing-data/total-correction preprocessing, using row-level
      `site_sequence` support from supplied values and/or supported derivation.

    When `policy="build_from_metadata"`:

    - `missing_data_policy` controls row retention before duplicate handling:
      - `"drop_any_missing"`: keep only rows with complete phospho values.
        This complete-case policy is the only supported public mode for
        construction of strict `AnalysisReadyPhosphoDataset` outputs.
        Retained-missingness site-matrix modes
        (`"retain_missing"`, `"require_min_observed_values"`) are internal
        compatibility modes and are rejected in the public builder lane.
    - `duplicate_site_strategy` controls duplicate-site collapse:
      - `"max_mean_signal"` (legacy default): keep row with strongest signal.
      - `"first"`: keep first encountered row for each duplicate site.
      - `"aggregate_mean"`: aggregate duplicate phospho values by column mean.
      - `"aggregate_median"`: aggregate duplicate phospho values by column median.
      - `"error"`: fail when duplicate constructed site identifiers are present.

    `minimum_observed_values` is retained for internal preprocessing compatibility
    and must stay unset in the supported public complete-case builder lane.
    """

    policy: DatasetSiteMatrixPolicy = DATASET_SITE_MATRIX_POLICY_AS_INPUT
    duplicate_site_strategy: DatasetSiteMatrixDuplicateSiteStrategy = (
        DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL
    )
    missing_data_policy: DatasetSiteMatrixMissingDataPolicy = (
        DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
    )
    minimum_observed_values: int | None = None


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


@dataclass(frozen=True, slots=True)
class DatasetPreprocessingConfig:
    """Public grouped preprocessing policy for dataset building.

    The builder owns this policy surface. Groups are intentionally separated so
    supported preprocessing science remains user-visible:

    - `missing_data`: missing-value handling policy.
    - `total_protein_correction`: total/protein correction policy.
    - `site_matrix`: site-matrix construction policy.
    - `comparisons`: comparison-building policy.
    """

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
    diagnostic scoring outputs (`motif_scores`, `weights`). The authoritative
    downstream lane (`combined_scores` with profile fallback) is always computed.

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


@dataclass(frozen=True, slots=True)
class KinasePredictionConfig:
    """Public prediction-stage configuration.

    `mode` selects the prediction lane:

    - `"deterministic_ranking"`: deterministic top-kinase selection from
      downstream scores (legacy rewrite shortcut lane).
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


__all__ = [
    "DATASET_COMPARISON_BUILDING_POLICIES",
    "DATASET_COMPARISON_BUILDING_DEFAULT_SAMPLE_GROUP_COLUMN",
    "DATASET_COMPARISON_BUILDING_POLICY_NONE",
    "DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS",
    "DATASET_MISSING_DATA_POLICIES",
    "DATASET_MISSING_DATA_POLICY_FORBID",
    "DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN",
    "DATASET_SITE_MATRIX_POLICIES",
    "DATASET_SITE_MATRIX_DUPLICATE_STRATEGIES",
    "DATASET_SITE_MATRIX_MISSING_DATA_POLICIES",
    "DATASET_SITE_MATRIX_POLICY_AS_INPUT",
    "DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA",
    "DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL",
    "DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_FIRST",
    "DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEAN",
    "DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_AGGREGATE_MEDIAN",
    "DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_ERROR",
    "DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_NONE",
    "DATASET_TOTAL_PROTEIN_CORRECTION_POLICY_RATIO_TO_TOTAL",
    "DatasetComparisonBuildingConfig",
    "DatasetComparisonPair",
    "DatasetComparisonBuildingPolicy",
    "DatasetMissingDataConfig",
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
    "DatasetSiteMatrixDuplicateSiteStrategy",
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
