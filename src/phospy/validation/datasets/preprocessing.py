"""Dataset-build preprocessing config validation."""

from __future__ import annotations

from phospy.api.configs import (
    DATASET_COMPARISON_BUILDING_POLICIES,
    DATASET_COMPARISON_BUILDING_POLICY_NONE,
    DATASET_COMPARISON_BUILDING_POLICY_SAMPLE_METADATA_PAIRS,
    DATASET_MISSING_DATA_POLICIES,
    DATASET_MISSING_DATA_POLICY_FORBID,
    DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGIES,
    DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL,
    DATASET_SITE_MATRIX_MISSING_DATA_POLICIES,
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    DATASET_SITE_MATRIX_POLICIES,
    DATASET_SITE_MATRIX_POLICY_AS_INPUT,
    DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES,
    DatasetComparisonBuildingConfig,
    DatasetMissingDataConfig,
    DatasetPreprocessingConfig,
    DatasetSiteMatrixConfig,
    DatasetTotalProteinCorrectionConfig,
)
from phospy.errors.input import PhosPyInputError

_INCOMPATIBLE_SITE_MATRIX_MISSING_DATA_POLICIES = frozenset(
    {"retain_missing", "require_min_observed_values"}
)
_SUPPORTED_SITE_MATRIX_MISSING_DATA_POLICY = (
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
)


class DatasetPreprocessingConfigValidator:
    """Validate the public config interpreted into the internal preprocessing plan."""

    def run(self, config: DatasetPreprocessingConfig) -> DatasetPreprocessingConfig:
        if not isinstance(config, DatasetPreprocessingConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config must be a "
                "DatasetPreprocessingConfig"
            )

        self._validate_missing_data(config.missing_data)
        self._validate_total_protein_correction(config.total_protein_correction)
        self._validate_site_matrix(config.site_matrix)
        self._validate_comparisons(config.comparisons)
        return config

    def _validate_missing_data(self, config: DatasetMissingDataConfig) -> None:
        if not isinstance(config, DatasetMissingDataConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data must be a "
                "DatasetMissingDataConfig"
            )

        policy = config.policy
        if policy not in DATASET_MISSING_DATA_POLICIES:
            supported = ", ".join(sorted(DATASET_MISSING_DATA_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.policy "
                f"must be one of: {supported}"
            )

        min_observed_values = config.min_observed_values
        if policy == DATASET_MISSING_DATA_POLICY_FORBID:
            if min_observed_values is not None:
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.min_observed_values must be "
                    "None when missing_data.policy='forbid'"
                )
            return

        if policy == DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN:
            if not isinstance(min_observed_values, int):
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

    def _validate_total_protein_correction(
        self, config: DatasetTotalProteinCorrectionConfig
    ) -> None:
        if not isinstance(config, DatasetTotalProteinCorrectionConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction "
                "must be a DatasetTotalProteinCorrectionConfig"
            )

        policy = config.policy
        if policy not in DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES:
            supported = ", ".join(sorted(DATASET_TOTAL_PROTEIN_CORRECTION_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.total_protein_correction."
                f"policy must be one of: {supported}"
            )

    def _validate_site_matrix(self, config: DatasetSiteMatrixConfig) -> None:
        if not isinstance(config, DatasetSiteMatrixConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix must be a "
                "DatasetSiteMatrixConfig"
            )

        policy = config.policy
        if policy not in DATASET_SITE_MATRIX_POLICIES:
            supported = ", ".join(sorted(DATASET_SITE_MATRIX_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix.policy "
                f"must be one of: {supported}"
            )
        duplicate_site_strategy = config.duplicate_site_strategy
        if duplicate_site_strategy not in DATASET_SITE_MATRIX_DUPLICATE_STRATEGIES:
            supported_duplicates = ", ".join(
                sorted(DATASET_SITE_MATRIX_DUPLICATE_STRATEGIES)
            )
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "duplicate_site_strategy must be one of: "
                f"{supported_duplicates}"
            )

        missing_data_policy = config.missing_data_policy
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

        minimum_observed_values = config.minimum_observed_values
        if minimum_observed_values is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "minimum_observed_values is not supported for strict "
                "AnalysisReadyPhosphoDataset construction and must be None"
            )

        if policy != DATASET_SITE_MATRIX_POLICY_AS_INPUT:
            return

        if (
            duplicate_site_strategy
            != DATASET_SITE_MATRIX_DUPLICATE_STRATEGY_MAX_MEAN_SIGNAL
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "duplicate_site_strategy is only valid when "
                "site_matrix.policy='build_from_metadata'"
            )

        if (
            missing_data_policy
            != DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
        ):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "missing_data_policy is only valid when "
                "site_matrix.policy='build_from_metadata'"
            )

        if minimum_observed_values is not None:
            raise PhosPyInputError(
                "dataset build request preprocessing_config.site_matrix."
                "minimum_observed_values is only valid when "
                "site_matrix.policy='build_from_metadata'"
            )

    def _validate_comparisons(self, config: DatasetComparisonBuildingConfig) -> None:
        if not isinstance(config, DatasetComparisonBuildingConfig):
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons must be a "
                "DatasetComparisonBuildingConfig"
            )

        policy = config.policy
        if policy not in DATASET_COMPARISON_BUILDING_POLICIES:
            supported = ", ".join(sorted(DATASET_COMPARISON_BUILDING_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons.policy "
                f"must be one of: {supported}"
            )
        sample_group_column = config.sample_group_column
        if not isinstance(sample_group_column, str) or not sample_group_column.strip():
            raise PhosPyInputError(
                "dataset build request preprocessing_config.comparisons."
                "sample_group_column must be a non-empty string"
            )
        pairs = config.pairs
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
