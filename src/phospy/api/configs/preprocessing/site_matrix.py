"""Site-matrix preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.errors.input import PhosPyInputError

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
DatasetSiteMatrixMissingDataPolicy = Literal["drop_any_missing"]
DATASET_SITE_MATRIX_MISSING_DATA_POLICIES = frozenset(
    {
        DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING,
    }
)

_INCOMPATIBLE_SITE_MATRIX_MISSING_DATA_POLICIES = frozenset(
    {"retain_missing", "require_min_observed_values"}
)
_SUPPORTED_SITE_MATRIX_MISSING_DATA_POLICY = (
    DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING
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


__all__ = [
    "DATASET_SITE_MATRIX_DUPLICATE_POLICIES",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEAN",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_AGGREGATE_MEDIAN",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_ERROR",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_FIRST",
    "DATASET_SITE_MATRIX_DUPLICATE_POLICY_MAX_MEAN_SIGNAL",
    "DATASET_SITE_MATRIX_MISSING_DATA_POLICIES",
    "DATASET_SITE_MATRIX_MISSING_DATA_POLICY_DROP_ANY_MISSING",
    "DATASET_SITE_MATRIX_POLICIES",
    "DATASET_SITE_MATRIX_POLICY_AS_INPUT",
    "DATASET_SITE_MATRIX_POLICY_BUILD_FROM_METADATA",
    "DatasetSiteMatrixConfig",
    "DatasetSiteMatrixDuplicateSitePolicy",
    "DatasetSiteMatrixMissingDataPolicy",
    "DatasetSiteMatrixPolicy",
]
