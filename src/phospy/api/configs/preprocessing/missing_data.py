"""Missing-data preprocessing policy configuration."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

from phospy.errors.input import PhosPyInputError

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
    - `"impute_knn"`: nearest-neighbour imputation using scikit-learn
      `KNNImputer` with fixed `metric="nan_euclidean"` and explicit row-level
      missingness filtering.

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
        policy = self.policy
        if policy not in DATASET_MISSING_DATA_POLICIES:
            supported = ", ".join(sorted(DATASET_MISSING_DATA_POLICIES))
            raise PhosPyInputError(
                "dataset build request preprocessing_config.missing_data.policy "
                f"must be one of: {supported}"
            )

        min_observed_values = self.min_observed_values
        q = self.q
        width = self.width
        seed = self.seed
        k = self.k
        distance = self.distance
        max_missing_fraction_per_row = self.max_missing_fraction_per_row
        if policy == DATASET_MISSING_DATA_POLICY_FORBID:
            if min_observed_values is not None:
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.min_observed_values must be "
                    "None when missing_data.policy='forbid'"
                )
            if (
                q is not None
                or width is not None
                or seed is not None
                or k is not None
                or distance is not None
                or max_missing_fraction_per_row is not None
            ):
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.q, .width, .seed, .k, "
                    ".distance, and "
                    ".max_missing_fraction_per_row must be None when "
                    "missing_data.policy='forbid'"
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
            if (
                q is not None
                or width is not None
                or seed is not None
                or k is not None
                or distance is not None
                or max_missing_fraction_per_row is not None
            ):
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.q, .width, .seed, .k, "
                    ".distance, and "
                    ".max_missing_fraction_per_row must be None when "
                    "missing_data.policy='impute_row_median'"
                )
            return
        if policy == DATASET_MISSING_DATA_POLICY_IMPUTE_MINPROB:
            if min_observed_values is not None:
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.min_observed_values must be "
                    "None when missing_data.policy='impute_minprob'"
                )
            if isinstance(q, bool) or not isinstance(q, (int, float)):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.q must "
                    "be a float when missing_data.policy='impute_minprob'"
                )
            q_value = float(q)
            if not math.isfinite(q_value) or not (0.0 < q_value < 0.5):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.q must "
                    "satisfy 0 < q < 0.5 when missing_data.policy='impute_minprob'"
                )
            if isinstance(width, bool) or not isinstance(width, (int, float)):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.width "
                    "must be a float when missing_data.policy='impute_minprob'"
                )
            width_value = float(width)
            if not math.isfinite(width_value) or not (0.0 < width_value <= 1.0):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.width "
                    "must satisfy 0 < width <= 1.0 when "
                    "missing_data.policy='impute_minprob'"
                )
            if isinstance(seed, bool) or not isinstance(seed, int):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.seed "
                    "must be an int when missing_data.policy='impute_minprob'"
                )
            if seed < 0:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.seed "
                    "must be greater than or equal to 0 when "
                    "missing_data.policy='impute_minprob'"
                )
            if isinstance(max_missing_fraction_per_row, bool) or not isinstance(
                max_missing_fraction_per_row, (int, float)
            ):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data."
                    "max_missing_fraction_per_row must be a float when "
                    "missing_data.policy='impute_minprob'"
                )
            max_missing_fraction_value = float(max_missing_fraction_per_row)
            if not math.isfinite(max_missing_fraction_value) or not (
                0.0 < max_missing_fraction_value <= 1.0
            ):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data."
                    "max_missing_fraction_per_row must satisfy 0 < value <= 1 when "
                    "missing_data.policy='impute_minprob'"
                )
            if k is not None or distance is not None:
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.k and .distance "
                    "must be None when missing_data.policy='impute_minprob'"
                )
            return
        if policy == DATASET_MISSING_DATA_POLICY_IMPUTE_KNN:
            if min_observed_values is not None:
                raise PhosPyInputError(
                    "dataset build request "
                    "preprocessing_config.missing_data.min_observed_values must be "
                    "None when missing_data.policy='impute_knn'"
                )
            if q is not None or width is not None or seed is not None:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.q, "
                    ".width, and .seed must be None when "
                    "missing_data.policy='impute_knn'"
                )
            if isinstance(k, bool) or not isinstance(k, int):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.k "
                    "must be an int when missing_data.policy='impute_knn'"
                )
            if k < 1:
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.k "
                    "must be greater than or equal to 1 when "
                    "missing_data.policy='impute_knn'"
                )
            if not isinstance(distance, str) or distance.strip() != "nan_euclidean":
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data.distance "
                    "must be 'nan_euclidean' when missing_data.policy='impute_knn'"
                )
            if isinstance(max_missing_fraction_per_row, bool) or not isinstance(
                max_missing_fraction_per_row, (int, float)
            ):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data."
                    "max_missing_fraction_per_row must be a float when "
                    "missing_data.policy='impute_knn'"
                )
            max_missing_fraction_value = float(max_missing_fraction_per_row)
            if not math.isfinite(max_missing_fraction_value) or not (
                0.0 < max_missing_fraction_value <= 1.0
            ):
                raise PhosPyInputError(
                    "dataset build request preprocessing_config.missing_data."
                    "max_missing_fraction_per_row must satisfy 0 < value <= 1 when "
                    "missing_data.policy='impute_knn'"
                )
            return
        raise PhosPyInputError(
            "dataset build request preprocessing_config contains an unsupported "
            "missing_data.policy"
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
