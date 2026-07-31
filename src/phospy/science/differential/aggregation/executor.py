"""Internal peptide-to-site differential estimate-combination execution.

The public post-hoc lane is withdrawn. This source is retained only for future
scientific design work and is not exported through supported facades.
"""

from __future__ import annotations

import math
from typing import cast

import numpy as np
import pandas as pd
from scipy import stats

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.aggregation.models import (
    PEPTIDE_DIFFERENTIAL_CONSISTENCY_POLICY,
    PEPTIDE_DIFFERENTIAL_CONSISTENCY_TOLERANCE_VERSION,
    PEPTIDE_DIFFERENTIAL_MAPPING_WEIGHT_POLICY_REJECT,
    PEPTIDE_DIFFERENTIAL_P_VALUE_ABS_TOLERANCE,
    PEPTIDE_DIFFERENTIAL_P_VALUE_REL_TOLERANCE,
    PEPTIDE_DIFFERENTIAL_STATISTIC_ABS_TOLERANCE,
    PEPTIDE_DIFFERENTIAL_STATISTIC_REL_TOLERANCE,
    PEPTIDE_TO_SITE_AGGREGATION_LEVEL_POSTHOC,
    PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE,
    PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES,
    PEPTIDE_TO_SITE_FIXED_EFFECT_APPROXIMATION_POLICY,
    PEPTIDE_TO_SITE_FIXED_EFFECT_MIN_ASYMPTOTIC_MODERATED_DF,
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE,
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH,
    PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P,
    PeptideDifferentialEstimateTable,
    PeptideToSiteAggregationConfig,
    PeptideToSiteAggregationResult,
)
from phospy.science.differential.aggregation.scientific_policies import (
    build_peptide_to_site_aggregation_policy,
)
from phospy.science.statistics.multiple_testing import adjust_p_values

_SINGLE_ESTIMATE_DEPENDENCE_ASSUMPTION = (
    "not_applicable_single_estimate_no_cross_peptide_combination"
)
_INDEPENDENT_SOURCE_DEPENDENCE_ASSUMPTION = (
    "independent_source_experiments_or_runs_no_same_sample_peptide_dependence"
)
_SINGLE_ESTIMATE_P_VALUE_METHOD = "original_two_sided_t_distribution_p_value"
_STOUFFER_P_VALUE_METHOD = "stouffer_weighted_z_from_signed_two_sided_input_p_values"
_FIXED_EFFECT_P_VALUE_METHOD = (
    "asymptotic_normal_fixed_effect_inverse_variance_df_ge_1000"
)
_T_DISTRIBUTION = "moderated_t"
_Z_DISTRIBUTION = "standard_normal_z"


class PeptideToSiteAggregationExecutor:
    """Execute internal peptide-level estimate-combination experiments."""

    def run_estimates(
        self,
        *,
        estimates: PeptideDifferentialEstimateTable,
        config: PeptideToSiteAggregationConfig,
        contrast_name: str,
    ) -> PeptideToSiteAggregationResult:
        if not isinstance(estimates, PeptideDifferentialEstimateTable):
            raise PhosPyInputError(
                "peptide-to-site aggregation requires a "
                "PeptideDifferentialEstimateTable"
            )
        if not isinstance(config, PeptideToSiteAggregationConfig):
            raise PhosPyInputError(
                "peptide-to-site aggregation config must be a "
                "PeptideToSiteAggregationConfig"
            )

        estimate_frame = estimates.to_dataframe()
        estimate_identity = _estimate_identity(estimate_frame)
        requested_contrast_name = str(contrast_name).strip()
        resolved_contrast_name = (
            estimate_identity["contrast_id"]
            if requested_contrast_name == "aggregated"
            else contrast_name
        )
        if str(resolved_contrast_name).strip() != estimate_identity["contrast_id"]:
            raise PhosPyInputError(
                "peptide-to-site aggregation contrast_name must match the input "
                "PeptideDifferentialEstimateTable contrast_id; "
                f"contrast_name={contrast_name!r}, "
                f"contrast_id={estimate_identity['contrast_id']!r}"
            )
        rows: list[dict[str, object]] = []
        withheld_below_minimum_count = 0
        for site_id, site_rows in estimate_frame.groupby("site_id", sort=True):
            combined = _combine_site_estimates(
                site_id=str(site_id),
                site_rows=site_rows,
                config=config,
            )
            if cast(int, combined["n_peptides_used"]) < int(
                config.min_estimates_per_site
            ):
                withheld_below_minimum_count += 1
            rows.append({"site_id": str(site_id), **combined})

        site_table = pd.DataFrame(rows).set_index("site_id", drop=True)
        site_table.index = pd.Index(site_table.index.astype(str), name="site_id")
        site_table.loc[:, "adj.P.Val"] = adjust_p_values(
            site_table.loc[:, "P.Value"].to_numpy(dtype=float),
            method=config.multiple_testing_method,
        )
        site_table.loc[:, "correction_method"] = str(config.multiple_testing_method)
        site_table = site_table.loc[
            :,
            [
                "contrast_id",
                "contrast_orientation",
                "effect_scale",
                "effect_unit",
                "model_estimator_id",
                "input_statistic_distribution",
                "input_uncertainty_method_version",
                "logFC",
                "standard_error",
                "uncertainty_statistic",
                "P.Value",
                "adj.P.Val",
                "residual_degrees_of_freedom",
                "moderated_degrees_of_freedom",
                "n_peptide_observations",
                "n_peptides_used",
                "source_experiment_ids",
                "peptide_to_site_mapping_policy",
                "multi_site_estimate_count",
                "aggregation_level",
                "dependence_assumption",
                "uncertainty_method",
                "correction_method",
                "p_value_method",
                "statistic_distribution",
            ],
        ]
        warnings: list[str] = []
        if withheld_below_minimum_count:
            warnings.append(
                "Some peptide-to-site estimate groups did not satisfy "
                "min_estimates_per_site and were emitted with missing statistics."
            )

        mapping_policies = tuple(
            sorted(
                {
                    str(value)
                    for value in estimate_frame.loc[
                        :,
                        "peptide_to_site_mapping_policy",
                    ].tolist()
                }
            )
        )
        policy = build_peptide_to_site_aggregation_policy(
            uncertainty_method=config.uncertainty_method,
            min_estimates_per_site=config.min_estimates_per_site,
            dependence_policy=config.dependence_policy,
            multiple_testing_method=config.multiple_testing_method,
            input_mapping_policies=mapping_policies,
            input_contrast_id=estimate_identity["contrast_id"],
            input_contrast_orientation=estimate_identity["contrast_orientation"],
            input_effect_scale=estimate_identity["effect_scale"],
            input_effect_unit=estimate_identity["effect_unit"],
            input_model_estimator_id=estimate_identity["model_estimator_id"],
            input_statistic_distribution=estimate_identity[
                "input_statistic_distribution"
            ],
            input_uncertainty_method_version=estimate_identity[
                "input_uncertainty_method_version"
            ],
            consistency_policy=PEPTIDE_DIFFERENTIAL_CONSISTENCY_POLICY,
            consistency_tolerance_version=(
                PEPTIDE_DIFFERENTIAL_CONSISTENCY_TOLERANCE_VERSION
            ),
            approximation_policy=_approximation_policy(config.uncertainty_method),
            mapping_weight_policy=PEPTIDE_DIFFERENTIAL_MAPPING_WEIGHT_POLICY_REJECT,
        )
        provenance = {
            "aggregation_level": PEPTIDE_TO_SITE_AGGREGATION_LEVEL_POSTHOC,
            "single_estimate_level": PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE,
            "preferred_phospy_lane": (
                "resolve peptide evidence at sample-intensity level before "
                "differential model fitting"
            ),
            "uncertainty_method": config.uncertainty_method,
            "dependence_policy": config.dependence_policy,
            "dependence_assumptions": (_INDEPENDENT_SOURCE_DEPENDENCE_ASSUMPTION,),
            "multiple_testing_method": config.multiple_testing_method,
            "input_contrast_id": estimate_identity["contrast_id"],
            "input_contrast_orientation": estimate_identity["contrast_orientation"],
            "input_effect_scale": estimate_identity["effect_scale"],
            "input_effect_unit": estimate_identity["effect_unit"],
            "input_model_estimator_id": estimate_identity["model_estimator_id"],
            "input_statistic_distribution": estimate_identity[
                "input_statistic_distribution"
            ],
            "input_uncertainty_method_version": estimate_identity[
                "input_uncertainty_method_version"
            ],
            "consistency_policy": PEPTIDE_DIFFERENTIAL_CONSISTENCY_POLICY,
            "consistency_tolerance_version": (
                PEPTIDE_DIFFERENTIAL_CONSISTENCY_TOLERANCE_VERSION
            ),
            "statistic_consistency_abs_tolerance": (
                PEPTIDE_DIFFERENTIAL_STATISTIC_ABS_TOLERANCE
            ),
            "statistic_consistency_rel_tolerance": (
                PEPTIDE_DIFFERENTIAL_STATISTIC_REL_TOLERANCE
            ),
            "p_value_consistency_abs_tolerance": (
                PEPTIDE_DIFFERENTIAL_P_VALUE_ABS_TOLERANCE
            ),
            "p_value_consistency_rel_tolerance": (
                PEPTIDE_DIFFERENTIAL_P_VALUE_REL_TOLERANCE
            ),
            "approximation_policy": _approximation_policy(config.uncertainty_method),
            "mapping_weight_policy": PEPTIDE_DIFFERENTIAL_MAPPING_WEIGHT_POLICY_REJECT,
            "input_peptide_estimates": int(estimate_frame.shape[0]),
            "output_site_rows": int(site_table.shape[0]),
            "min_estimates_per_site": int(config.min_estimates_per_site),
            "withheld_below_minimum_site_count": int(withheld_below_minimum_count),
            "peptide_to_site_mapping_policies": mapping_policies,
            "multi_site_estimate_count": int(
                _multi_site_estimate_count(estimate_frame)
            ),
            "same_source_duplicate_policy": (
                "reject_same_site_estimates_with_duplicate_source_experiment_id"
            ),
            "t_to_z_policy": (
                "finite_df_t_evidence_is_converted_to_z_only_through_signed_"
                "two_sided_p_values"
            ),
        }
        return PeptideToSiteAggregationResult(
            contrast_name=str(resolved_contrast_name),
            table=site_table,
            warnings=tuple(warnings),
            provenance=provenance,
            scientific_policies=(policy,),
        )

    def run_table(
        self,
        *,
        estimate_table: PeptideDifferentialEstimateTable | pd.DataFrame | None = None,
        estimates: PeptideDifferentialEstimateTable | pd.DataFrame | None = None,
        config: PeptideToSiteAggregationConfig,
        contrast_name: str,
    ) -> PeptideToSiteAggregationResult:
        resolved = estimate_table if estimate_table is not None else estimates
        if resolved is None:
            raise PhosPyInputError(
                "peptide-to-site aggregation requires estimate_table or estimates"
            )
        typed_estimates = (
            resolved
            if isinstance(resolved, PeptideDifferentialEstimateTable)
            else PeptideDifferentialEstimateTable(cast(pd.DataFrame, resolved))
        )
        return self.run_estimates(
            estimates=typed_estimates,
            config=config,
            contrast_name=contrast_name,
        )


def signed_z_from_t_statistic(
    t_statistic: float,
    degrees_of_freedom: float,
) -> float:
    """Convert a finite-df t statistic to signed z through its two-sided p-value."""

    statistic = _finite_float(
        t_statistic,
        field_name="t_statistic",
    )
    df = _positive_finite_float(
        degrees_of_freedom,
        field_name="degrees_of_freedom",
    )
    if statistic == 0.0:
        return 0.0
    p_value = float(2.0 * stats.t.sf(abs(statistic), df=df))
    return signed_z_from_two_sided_p_value(
        p_value,
        sign_source=statistic,
    )


def signed_z_from_two_sided_p_value(
    p_value: float,
    *,
    sign_source: float,
) -> float:
    """Return signed standard-normal z from a two-sided p-value and sign source."""

    p = _p_value(p_value, field_name="p_value")
    sign = _sign(sign_source)
    if sign == 0.0:
        if p == 1.0:
            return 0.0
        raise PhosPyInputError(
            "cannot derive a signed z value from a non-null two-sided p-value "
            "when both effect and statistic signs are zero"
        )
    if p == 1.0:
        return 0.0
    return float(sign * stats.norm.isf(p / 2.0))


def _estimate_identity(frame: pd.DataFrame) -> dict[str, str]:
    return {
        "contrast_id": _first_column_value(frame, "contrast_id"),
        "contrast_orientation": _first_column_value(frame, "contrast_orientation"),
        "effect_scale": _first_column_value(frame, "effect_scale"),
        "effect_unit": _first_column_value(frame, "effect_unit"),
        "model_estimator_id": _first_column_value(frame, "model_estimator_id"),
        "input_statistic_distribution": _first_column_value(
            frame,
            "statistic_distribution",
        ),
        "input_uncertainty_method_version": _first_column_value(
            frame,
            "uncertainty_method_version",
        ),
    }


def _identity_result_columns(site_rows: pd.DataFrame) -> dict[str, object]:
    identity = _estimate_identity(site_rows)
    return {
        "contrast_id": identity["contrast_id"],
        "contrast_orientation": identity["contrast_orientation"],
        "effect_scale": identity["effect_scale"],
        "effect_unit": identity["effect_unit"],
        "model_estimator_id": identity["model_estimator_id"],
        "input_statistic_distribution": identity["input_statistic_distribution"],
        "input_uncertainty_method_version": identity[
            "input_uncertainty_method_version"
        ],
    }


def _first_column_value(frame: pd.DataFrame, column_name: str) -> str:
    column: pd.Series = frame.loc[:, column_name]  # pyright: ignore[reportUnknownMemberType, reportArgumentType, reportCallIssue] - pandas-stubs cannot express scalar-column DataFrame.loc selection.
    return str(column.iloc[0])


def _approximation_policy(uncertainty_method: str) -> str:
    if (
        uncertainty_method
        == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE
    ):
        return PEPTIDE_TO_SITE_FIXED_EFFECT_APPROXIMATION_POLICY
    if uncertainty_method == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P:
        return "finite_df_t_preserved_by_signed_two_sided_p_to_z_stouffer_v1"
    if (
        uncertainty_method
        == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH
    ):
        return "none_single_estimate_passthrough_original_finite_df_t_v1"
    return "unsupported_uncertainty_method"


def _combine_site_estimates(
    *,
    site_id: str,
    site_rows: pd.DataFrame,
    config: PeptideToSiteAggregationConfig,
) -> dict[str, object]:
    n_observed = int(site_rows.shape[0])
    if n_observed > 1:
        _enforce_independent_source_estimates(site_id=site_id, site_rows=site_rows)
    if n_observed < int(config.min_estimates_per_site):
        return _missing_site_result(
            site_rows=site_rows,
            n_observed=n_observed,
            n_used=n_observed,
            reason_uncertainty_method=config.uncertainty_method,
            correction_method=config.multiple_testing_method,
        )
    if n_observed == 1:
        return _single_estimate_result(
            site_rows=site_rows,
            correction_method=config.multiple_testing_method,
        )
    if (
        config.uncertainty_method
        == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH
    ):
        raise PhosPyInputError(
            "single_estimate_passthrough cannot combine multiple peptide "
            f"estimates for site_id={site_id!r}; choose a supported "
            "multi-estimate method with independent source estimates"
        )

    if (
        config.uncertainty_method
        == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P
    ):
        return _stouffer_signed_p_result(
            site_rows=site_rows,
            correction_method=config.multiple_testing_method,
        )
    if (
        config.uncertainty_method
        == PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE
    ):
        return _fixed_effect_inverse_variance_result(
            site_rows=site_rows,
            correction_method=config.multiple_testing_method,
        )
    raise PhosPyInputError(
        f"unsupported peptide-to-site uncertainty_method: {config.uncertainty_method!r}"
    )


def _single_estimate_result(
    *,
    site_rows: pd.DataFrame,
    correction_method: str,
) -> dict[str, object]:
    row = site_rows.iloc[0]
    return {
        **_identity_result_columns(site_rows),
        "logFC": float(row["effect"]),
        "standard_error": float(row["standard_error"]),
        "uncertainty_statistic": float(row["statistic"]),
        "P.Value": float(row["p_value"]),
        "residual_degrees_of_freedom": float(row["residual_degrees_of_freedom"]),
        "moderated_degrees_of_freedom": float(row["moderated_degrees_of_freedom"]),
        "n_peptide_observations": int(site_rows.shape[0]),
        "n_peptides_used": 1,
        "source_experiment_ids": str(row["source_experiment_id"]),
        "peptide_to_site_mapping_policy": str(row["peptide_to_site_mapping_policy"]),
        "multi_site_estimate_count": _multi_site_estimate_count(site_rows),
        "aggregation_level": PEPTIDE_TO_SITE_AGGREGATION_LEVEL_SINGLE_ESTIMATE,
        "dependence_assumption": _SINGLE_ESTIMATE_DEPENDENCE_ASSUMPTION,
        "uncertainty_method": PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_SINGLE_ESTIMATE_PASSTHROUGH,
        "correction_method": str(correction_method),
        "p_value_method": _SINGLE_ESTIMATE_P_VALUE_METHOD,
        "statistic_distribution": _T_DISTRIBUTION,
    }


def _stouffer_signed_p_result(
    *,
    site_rows: pd.DataFrame,
    correction_method: str,
) -> dict[str, object]:
    p_values = site_rows.loc[:, "p_value"].to_numpy(dtype=float)
    statistics = site_rows.loc[:, "statistic"].to_numpy(dtype=float)
    effects = site_rows.loc[:, "effect"].to_numpy(dtype=float)
    standard_errors = site_rows.loc[:, "standard_error"].to_numpy(dtype=float)
    sign_sources = np.asarray(
        [
            _sign_source(effect=float(effect), statistic=float(statistic))
            for effect, statistic in zip(effects, statistics, strict=True)
        ],
        dtype=float,
    )
    z_values = np.asarray(
        [
            signed_z_from_two_sided_p_value(
                float(p_value),
                sign_source=float(sign_source),
            )
            for p_value, sign_source in zip(p_values, sign_sources, strict=True)
        ],
        dtype=float,
    )
    weights_for_z = 1.0 / standard_errors
    denominator = math.sqrt(float(np.sum(weights_for_z**2)))
    if denominator <= 0.0 or not math.isfinite(denominator):
        raise PhosPyInputError(
            "peptide-to-site Stouffer combination requires finite positive "
            "standard-error weights"
        )
    combined_z = float(np.sum(weights_for_z * z_values) / denominator)
    p_value = float(2.0 * stats.norm.sf(abs(combined_z)))
    effect, standard_error = _inverse_variance_effect_summary(site_rows)
    return {
        **_identity_result_columns(site_rows),
        "logFC": effect,
        "standard_error": standard_error,
        "uncertainty_statistic": combined_z,
        "P.Value": p_value,
        "residual_degrees_of_freedom": _minimum_column(
            site_rows, "residual_degrees_of_freedom"
        ),
        "moderated_degrees_of_freedom": _minimum_column(
            site_rows,
            "moderated_degrees_of_freedom",
        ),
        "n_peptide_observations": int(site_rows.shape[0]),
        "n_peptides_used": int(site_rows.shape[0]),
        "source_experiment_ids": _joined_unique_sorted(
            site_rows, "source_experiment_id"
        ),
        "peptide_to_site_mapping_policy": _joined_unique_sorted(
            site_rows,
            "peptide_to_site_mapping_policy",
        ),
        "multi_site_estimate_count": _multi_site_estimate_count(site_rows),
        "aggregation_level": PEPTIDE_TO_SITE_AGGREGATION_LEVEL_POSTHOC,
        "dependence_assumption": _INDEPENDENT_SOURCE_DEPENDENCE_ASSUMPTION,
        "uncertainty_method": PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_STOUFFER_SIGNED_P,
        "correction_method": str(correction_method),
        "p_value_method": _STOUFFER_P_VALUE_METHOD,
        "statistic_distribution": _Z_DISTRIBUTION,
    }


def _fixed_effect_inverse_variance_result(
    *,
    site_rows: pd.DataFrame,
    correction_method: str,
) -> dict[str, object]:
    _enforce_fixed_effect_asymptotic_eligibility(site_rows)
    effect, standard_error = _inverse_variance_effect_summary(site_rows)
    z_value = effect / standard_error
    p_value = float(2.0 * stats.norm.sf(abs(z_value)))
    return {
        **_identity_result_columns(site_rows),
        "logFC": effect,
        "standard_error": standard_error,
        "uncertainty_statistic": z_value,
        "P.Value": p_value,
        "residual_degrees_of_freedom": _minimum_column(
            site_rows, "residual_degrees_of_freedom"
        ),
        "moderated_degrees_of_freedom": _minimum_column(
            site_rows,
            "moderated_degrees_of_freedom",
        ),
        "n_peptide_observations": int(site_rows.shape[0]),
        "n_peptides_used": int(site_rows.shape[0]),
        "source_experiment_ids": _joined_unique_sorted(
            site_rows, "source_experiment_id"
        ),
        "peptide_to_site_mapping_policy": _joined_unique_sorted(
            site_rows,
            "peptide_to_site_mapping_policy",
        ),
        "multi_site_estimate_count": _multi_site_estimate_count(site_rows),
        "aggregation_level": PEPTIDE_TO_SITE_AGGREGATION_LEVEL_POSTHOC,
        "dependence_assumption": _INDEPENDENT_SOURCE_DEPENDENCE_ASSUMPTION,
        "uncertainty_method": (
            PEPTIDE_TO_SITE_UNCERTAINTY_METHOD_FIXED_EFFECT_INVERSE_VARIANCE
        ),
        "correction_method": str(correction_method),
        "p_value_method": _FIXED_EFFECT_P_VALUE_METHOD,
        "statistic_distribution": _Z_DISTRIBUTION,
    }


def _inverse_variance_effect_summary(site_rows: pd.DataFrame) -> tuple[float, float]:
    effect = site_rows.loc[:, "effect"].to_numpy(dtype=float)
    standard_error = site_rows.loc[:, "standard_error"].to_numpy(dtype=float)
    weights = 1.0 / (standard_error**2)
    denominator = float(np.sum(weights))
    if denominator <= 0.0 or not math.isfinite(denominator):
        raise PhosPyInputError(
            "peptide-to-site estimate combination requires finite positive "
            "inverse-variance weights"
        )
    weighted_effect = float(np.sum(weights * effect) / denominator)
    combined_standard_error = float(math.sqrt(1.0 / denominator))
    return weighted_effect, combined_standard_error


def _enforce_fixed_effect_asymptotic_eligibility(site_rows: pd.DataFrame) -> None:
    minimum_moderated_df = _minimum_column(site_rows, "moderated_degrees_of_freedom")
    if minimum_moderated_df < PEPTIDE_TO_SITE_FIXED_EFFECT_MIN_ASYMPTOTIC_MODERATED_DF:
        raise PhosPyInputError(
            "fixed_effect_inverse_variance_independent requires the documented "
            "asymptotic-normal input envelope; all moderated_degrees_of_freedom "
            "values must be >= "
            f"{PEPTIDE_TO_SITE_FIXED_EFFECT_MIN_ASYMPTOTIC_MODERATED_DF:g}. "
            f"Observed minimum moderated_degrees_of_freedom={minimum_moderated_df!r}."
        )


def _missing_site_result(
    *,
    site_rows: pd.DataFrame,
    n_observed: int,
    n_used: int,
    reason_uncertainty_method: str,
    correction_method: str,
) -> dict[str, object]:
    return {
        **_identity_result_columns(site_rows),
        "logFC": float("nan"),
        "standard_error": float("nan"),
        "uncertainty_statistic": float("nan"),
        "P.Value": float("nan"),
        "residual_degrees_of_freedom": float("nan"),
        "moderated_degrees_of_freedom": float("nan"),
        "n_peptide_observations": int(n_observed),
        "n_peptides_used": int(n_used),
        "source_experiment_ids": _joined_unique_sorted(
            site_rows, "source_experiment_id"
        ),
        "peptide_to_site_mapping_policy": _joined_unique_sorted(
            site_rows,
            "peptide_to_site_mapping_policy",
        ),
        "multi_site_estimate_count": _multi_site_estimate_count(site_rows),
        "aggregation_level": PEPTIDE_TO_SITE_AGGREGATION_LEVEL_POSTHOC,
        "dependence_assumption": _INDEPENDENT_SOURCE_DEPENDENCE_ASSUMPTION,
        "uncertainty_method": str(reason_uncertainty_method),
        "correction_method": str(correction_method),
        "p_value_method": "not_computed_minimum_evidence_not_met",
        "statistic_distribution": "not_computed",
    }


def _enforce_independent_source_estimates(
    *,
    site_id: str,
    site_rows: pd.DataFrame,
) -> None:
    dependence_values = tuple(
        str(value) for value in site_rows.loc[:, "dependence_policy"].tolist()
    )
    unsupported = tuple(
        value
        for value in dependence_values
        if value != PEPTIDE_TO_SITE_DEPENDENCE_POLICY_INDEPENDENT_SOURCES
    )
    if unsupported:
        raise PhosPyInputError(
            "same-sample or otherwise correlated peptide estimates cannot be "
            "combined by the withdrawn post-hoc peptide-to-site lane; "
            f"site_id={site_id!r}, dependence_policy_values={sorted(set(unsupported))}. "
            "Resolve peptide evidence at sample-intensity level before fitting "
            "the differential model, or add a supported dependence-aware method."
        )
    source_ids = site_rows.loc[:, "source_experiment_id"].astype(str)
    duplicated_sources = sorted(
        set(source_ids.loc[source_ids.duplicated(keep=False)].tolist())
    )
    if duplicated_sources:
        raise PhosPyInputError(
            "same-experiment peptide estimates for one site are rejected by the "
            "withdrawn post-hoc peptide-to-site lane because same-sample peptide "
            "dependence is not modelled; "
            f"site_id={site_id!r}, duplicated_source_experiment_id_values="
            f"{duplicated_sources}"
        )


def _sign_source(*, effect: float, statistic: float) -> float:
    if statistic != 0.0:
        return statistic
    if effect != 0.0:
        return effect
    return 0.0


def _joined_unique_sorted(frame: pd.DataFrame, column_name: str) -> str:
    values = sorted({str(value) for value in frame.loc[:, column_name].tolist()})
    return "|".join(values)


def _minimum_column(frame: pd.DataFrame, column_name: str) -> float:
    return float(np.min(frame.loc[:, column_name].to_numpy(dtype=float)))


def _multi_site_estimate_count(frame: pd.DataFrame) -> int:
    if "mapping_uncertainty" in frame.columns:
        return int(frame.loc[:, "mapping_uncertainty"].astype(bool).sum())
    return int(
        sum(
            1
            for site_id in frame.loc[:, "site_id"].astype(str).tolist()
            if "," in site_id
        )
    )


def _finite_float(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise PhosPyInputError(f"{field_name} must be numeric")
    resolved = float(value)
    if not math.isfinite(resolved):
        raise PhosPyInputError(f"{field_name} must be finite")
    return resolved


def _positive_finite_float(value: object, *, field_name: str) -> float:
    resolved = _finite_float(value, field_name=field_name)
    if resolved <= 0.0:
        raise PhosPyInputError(f"{field_name} must be > 0.0")
    return resolved


def _p_value(value: object, *, field_name: str) -> float:
    resolved = _finite_float(value, field_name=field_name)
    if resolved < 0.0 or resolved > 1.0:
        raise PhosPyInputError(f"{field_name} must be within [0.0, 1.0]")
    return resolved


def _sign(value: float) -> float:
    if value > 0.0:
        return 1.0
    if value < 0.0:
        return -1.0
    return 0.0


__all__: list[str] = []
