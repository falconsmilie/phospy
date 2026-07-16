"""Missing-data preprocessing policy configuration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.contracts.configs.preprocessing._validation import (
    validate_missing_data_config,
)

DATASET_MISSING_DATA_POLICY_FORBID = "forbid"
DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN = "impute_row_median"
DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB = "impute_minprob"
DATASET_MISSING_DATA_POLICY_IMPUTE_KNN = "impute_knn"
DatasetMissingDataPolicy = Literal[
    "forbid",
    "impute_row_median",
    "impute_minprob",
    "impute_knn",
]
DATASET_MISSING_DATA_POLICIES = frozenset(
    {
        DATASET_MISSING_DATA_POLICY_FORBID,
        DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
        DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB,
        DATASET_MISSING_DATA_POLICY_IMPUTE_KNN,
    }
)


@dataclass(frozen=True, slots=True)
class DatasetMissingDataConfig:
    """Public missing-data policy options for dataset building.

    - `"forbid"`: do not preprocess missing values (strict default behavior).
    - `"impute_row_median"`: for each site row, drop rows with fewer than
      `min_observed_values` quantified samples, then impute remaining missing
      values with that row's observed-value median.
    - `"impute_minprob"`: left-censored random imputation on log2-scale data
      using a MinProb-style column-wise normal model.
    - `"impute_knn"`: deterministic PhosPy-owned nearest-neighbour imputation
      with fixed `distance="nan_euclidean"`, explicit row-level missingness
      filtering, deterministic donor tie rules, and retained-column mean
      fallback when a missing cell has no donor with overlapping observed
      values.

    `min_observed_values` is required for `"impute_row_median"` and must stay
    unset for `"forbid"`, `"impute_minprob"`, and `"impute_knn"`.

    For `"impute_minprob"`, required parameters are:

    - `q` with `0 < q < 0.5` (recommended: `0.01`)
    - `width` with `0 < width <= 1.0` (recommended: `0.3`)
    - `seed` with integer `>= 0` (recommended: `12345`)
    - `max_missing_fraction_per_row` with `0 < value <= 1` (recommended: `0.5`)

    For `"impute_knn"`, required parameters are:

    - `k` with integer `>= 1`
    - `distance` fixed to `"nan_euclidean"`
    - `max_missing_fraction_per_row` with `0 < value <= 1`

    The KNN implementation is chunked and guarded for practical preprocessing
    scale. It drops rows above `max_missing_fraction_per_row`; for each
    retained missing cell, donors must be retained rows observed in the target
    column. Distances use the target/donor row's shared observed columns and are
    scaled by `n_columns / shared_observed_column_count`. Donor ties are ordered
    by `(str(row_id), original_position)`. The selected donor values are
    averaged without distance weighting; if no donor has any shared observed
    column, the retained-column mean is used. Requests outside the documented
    KNN performance envelope fail with an actionable `PhosPyInputError`.
    """

    policy: DatasetMissingDataPolicy = DATASET_MISSING_DATA_POLICY_FORBID
    min_observed_values: int | None = None
    q: float | None = None
    width: float | None = None
    seed: int | None = None
    k: int | None = None
    distance: str | None = None
    max_missing_fraction_per_row: float | None = None

    def __post_init__(self) -> None:
        validate_missing_data_config(
            policy=self.policy,
            min_observed_values=self.min_observed_values,
            q=self.q,
            width=self.width,
            seed=self.seed,
            k=self.k,
            distance=self.distance,
            max_missing_fraction_per_row=self.max_missing_fraction_per_row,
            supported_policies=DATASET_MISSING_DATA_POLICIES,
            policy_forbid=DATASET_MISSING_DATA_POLICY_FORBID,
            policy_impute_row_median=DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN,
            policy_impute_minprob=DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB,
            policy_impute_knn=DATASET_MISSING_DATA_POLICY_IMPUTE_KNN,
        )


__all__ = [
    "DATASET_MISSING_DATA_POLICIES",
    "DATASET_MISSING_DATA_POLICY_FORBID",
    "DATASET_MISSING_DATA_POLICY_IMPUTE_KNN",
    "DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB",
    "DATASET_MISSING_DATA_POLICY_IMPUTE_ROW_MEDIAN",
    "DatasetMissingDataConfig",
    "DatasetMissingDataPolicy",
]
