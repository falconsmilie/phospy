"""Experimental/internal peptide-to-site differential aggregation execution."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from phospy.errors.input import PhosPyInputError
from phospy.science.differential.aggregation.models import (
    PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE,
    PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META,
    PEPTIDE_TO_SITE_STRATEGY_INVERSE_VARIANCE_WEIGHTED,
    PEPTIDE_TO_SITE_STRATEGY_RANDOM_EFFECT_META,
    PEPTIDE_TO_SITE_STRATEGY_STOUFFER_Z,
    STOUFFER_WEIGHTING_INVERSE_VARIANCE,
    PeptideToSiteAggregationConfig,
    PeptideToSiteAggregationResult,
)
from phospy.science.differential.aggregation.scientific_policies import (
    build_peptide_to_site_aggregation_policy,
)
from phospy.science.differential.models import DifferentialAnalysisResult
from phospy.science.differential.multiple_testing import benjamini_hochberg
from phospy.science.evidence.models import PeptideEvidenceTable

_COMPATIBILITY_WARNING = (
    "compat_best_p_value selects the minimum peptide p-value per site and is "
    "provided for compatibility only; this strategy is statistically biased. "
    "All retained peptide-to-site differential aggregation strategies are "
    "experimental/internal and are not supported for production site-level "
    "inference."
)


class PeptideToSiteAggregationExecutor:
    """Execute unsupported experimental post-hoc site summaries."""

    def run_differential_result(
        self,
        *,
        differential_result: DifferentialAnalysisResult,
        evidence: PeptideEvidenceTable,
        config: PeptideToSiteAggregationConfig,
    ) -> dict[str, PeptideToSiteAggregationResult]:
        if not isinstance(differential_result, DifferentialAnalysisResult):
            raise PhosPyInputError(
                "peptide-to-site aggregation requires a DifferentialAnalysisResult"
            )
        results: dict[str, PeptideToSiteAggregationResult] = {}
        for contrast_name, table in differential_result.contrast_tables.items():
            results[contrast_name] = self.run_table(
                peptide_differential_table=table,
                evidence=evidence,
                config=config,
                contrast_name=contrast_name,
            )
        return results

    def run_table(
        self,
        *,
        peptide_differential_table: pd.DataFrame,
        evidence: PeptideEvidenceTable,
        config: PeptideToSiteAggregationConfig,
        contrast_name: str,
    ) -> PeptideToSiteAggregationResult:
        if not isinstance(evidence, PeptideEvidenceTable):
            raise PhosPyInputError(
                "peptide-to-site aggregation requires a PeptideEvidenceTable"
            )
        if not isinstance(config, PeptideToSiteAggregationConfig):
            raise PhosPyInputError(
                "peptide-to-site aggregation config must be a "
                "PeptideToSiteAggregationConfig"
            )
        table = _validate_peptide_differential_table(peptide_differential_table)
        feature_mapping = _build_feature_site_mapping(evidence=evidence)
        merged = _merge_peptide_table_with_site_mapping(
            peptide_differential_table=table,
            feature_site_mapping=feature_mapping,
        )
        if merged.empty:
            raise PhosPyInputError(
                "peptide-to-site aggregation found no overlapping peptide features "
                "between differential table index and evidence mapping"
            )

        rows: list[dict[str, float | int | str]] = []
        dropped_missing_variance = 0
        for site_id, site_rows in merged.groupby("site_id", sort=False):
            aggregated, dropped = _aggregate_site_rows(
                site_rows=site_rows,
                strategy=config.strategy,
                config=config,
            )
            dropped_missing_variance += dropped
            rows.append(
                {
                    "site_id": str(site_id),
                    "logFC": aggregated["logFC"],
                    "uncertainty_statistic": aggregated["uncertainty_statistic"],
                    "P.Value": aggregated["P.Value"],
                    "n_peptide_observations": int(site_rows.shape[0]),
                    "n_peptides_used": int(aggregated["n_peptides_used"]),
                }
            )
        site_table = pd.DataFrame(rows).set_index("site_id", drop=True)
        site_table.loc[:, "adj.P.Val"] = benjamini_hochberg(
            site_table.loc[:, "P.Value"].to_numpy(dtype=float)
        )
        site_table = site_table.loc[
            :,
            [
                "logFC",
                "uncertainty_statistic",
                "P.Value",
                "adj.P.Val",
                "n_peptide_observations",
                "n_peptides_used",
            ],
        ]
        warnings: list[str] = []
        compatibility_warning = (
            config.strategy == PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE
        )
        if compatibility_warning:
            warnings.append(_COMPATIBILITY_WARNING)
        policy = build_peptide_to_site_aggregation_policy(
            strategy=config.strategy,
            min_peptides_per_site=config.min_peptides_per_site,
            missing_variance_policy=config.missing_variance_policy,
            stouffer_weighting=config.stouffer_weighting,
            random_effect_tau2_floor=config.random_effect_tau2_floor,
            compatibility_mode_warning=compatibility_warning,
        )
        provenance = {
            "aggregation_strategy": config.strategy,
            "min_peptides_per_site": int(config.min_peptides_per_site),
            "missing_variance_policy": config.missing_variance_policy,
            "stouffer_weighting": config.stouffer_weighting,
            "random_effect_tau2_floor": float(config.random_effect_tau2_floor),
            "multi_site_handling": evidence.multi_site_policy_provenance(),
            "input_peptide_rows": int(merged.shape[0]),
            "output_site_rows": int(site_table.shape[0]),
            "dropped_missing_variance_rows": int(dropped_missing_variance),
        }
        return PeptideToSiteAggregationResult(
            contrast_name=contrast_name,
            table=site_table,
            warnings=tuple(warnings),
            provenance=provenance,
            scientific_policies=(policy,),
            _assume_owned=True,
        )


def _validate_peptide_differential_table(table: object) -> pd.DataFrame:
    if not isinstance(table, pd.DataFrame):
        raise PhosPyInputError("peptide differential input must be a pandas DataFrame")
    if table.empty:
        raise PhosPyInputError("peptide differential input must be non-empty")
    required_columns = ("logFC", "t", "P.Value")
    missing = [column for column in required_columns if column not in table.columns]
    if missing:
        joined = ", ".join(missing)
        raise PhosPyInputError(
            f"peptide differential input is missing required columns: {joined}"
        )
    if not table.index.is_unique:
        raise PhosPyInputError(
            "peptide differential input index must be unique peptide feature IDs"
        )
    for column_name in required_columns:
        if not pd.api.types.is_numeric_dtype(table[column_name]):
            raise PhosPyInputError(
                f"peptide differential column {column_name!r} must be numeric"
            )
    return table.copy(deep=True)


def _build_feature_site_mapping(*, evidence: PeptideEvidenceTable) -> pd.DataFrame:
    evidence_frame = evidence.to_dataframe()
    if evidence_frame.loc[:, "unique_feature_id"].duplicated().any():
        raise PhosPyInputError(
            "peptide evidence unique_feature_id values must be unique for "
            "peptide-to-site aggregation"
        )
    feature_rows = evidence_frame.loc[:, ["peptide_row_id", "unique_feature_id"]].copy(
        deep=True
    )
    mapping = evidence.site_mapping.to_dataframe()
    if mapping.empty:
        raise PhosPyInputError(
            "peptide evidence does not provide any peptide-to-site mappings"
        )
    merged = mapping.merge(
        feature_rows,
        how="left",
        on="peptide_row_id",
        indicator=True,
    )
    unresolved = merged.loc[
        merged.loc[:, "_merge"] != "both", "peptide_row_id"
    ].tolist()
    if unresolved:
        preview = ", ".join(repr(value) for value in unresolved[:5])
        suffix = "" if len(unresolved) <= 5 else " ..."
        raise PhosPyInputError(
            "site mapping references peptide_row_id values absent from evidence rows: "
            f"{preview}{suffix}"
        )
    resolved = merged.loc[:, ["unique_feature_id", "site_id"]].drop_duplicates()
    return resolved


def _merge_peptide_table_with_site_mapping(
    *,
    peptide_differential_table: pd.DataFrame,
    feature_site_mapping: pd.DataFrame,
) -> pd.DataFrame:
    feature_series = pd.Series(
        peptide_differential_table.index.astype("object"),
        index=peptide_differential_table.index.copy(),
        name="unique_feature_id",
    )
    joined = peptide_differential_table.copy(deep=True)
    joined.loc[:, "unique_feature_id"] = feature_series
    merged = feature_site_mapping.merge(
        joined,
        how="inner",
        on="unique_feature_id",
    )
    return merged


def _aggregate_site_rows(
    *,
    site_rows: pd.DataFrame,
    strategy: str,
    config: PeptideToSiteAggregationConfig,
) -> tuple[dict[str, float | int], int]:
    if strategy == PEPTIDE_TO_SITE_STRATEGY_COMPAT_BEST_P_VALUE:
        result = _aggregate_compat_best_p(site_rows=site_rows)
        return result, 0
    if strategy in (
        PEPTIDE_TO_SITE_STRATEGY_INVERSE_VARIANCE_WEIGHTED,
        PEPTIDE_TO_SITE_STRATEGY_FIXED_EFFECT_META,
    ):
        result, dropped = _aggregate_fixed_effect(
            site_rows=site_rows,
            min_peptides_per_site=config.min_peptides_per_site,
        )
        return result, dropped
    if strategy == PEPTIDE_TO_SITE_STRATEGY_RANDOM_EFFECT_META:
        result, dropped = _aggregate_random_effect(
            site_rows=site_rows,
            min_peptides_per_site=config.min_peptides_per_site,
            tau2_floor=config.random_effect_tau2_floor,
        )
        return result, dropped
    if strategy == PEPTIDE_TO_SITE_STRATEGY_STOUFFER_Z:
        result, dropped = _aggregate_stouffer(
            site_rows=site_rows,
            min_peptides_per_site=config.min_peptides_per_site,
            weighting=config.stouffer_weighting,
        )
        return result, dropped
    raise PhosPyInputError(
        f"unsupported peptide-to-site aggregation strategy: {strategy!r}"
    )


def _aggregate_compat_best_p(*, site_rows: pd.DataFrame) -> dict[str, float | int]:
    p_values = site_rows.loc[:, "P.Value"].to_numpy(dtype=float)
    finite = np.isfinite(p_values)
    if np.any(finite):
        finite_indices = np.flatnonzero(finite)
        min_pos = finite_indices[int(np.argmin(p_values[finite]))]
    else:
        min_pos = 0
    row = site_rows.iloc[int(min_pos)]
    return {
        "logFC": float(row["logFC"]),
        "uncertainty_statistic": float(row["t"]),
        "P.Value": float(row["P.Value"]),
        "n_peptides_used": 1,
    }


def _aggregate_fixed_effect(
    *,
    site_rows: pd.DataFrame,
    min_peptides_per_site: int,
) -> tuple[dict[str, float | int], int]:
    effect = site_rows.loc[:, "logFC"].to_numpy(dtype=float)
    variance, valid_mask = _derive_variance(site_rows=site_rows)
    dropped = int((~valid_mask).sum())
    if int(valid_mask.sum()) < min_peptides_per_site:
        return _nan_aggregate_result(n_used=int(valid_mask.sum())), dropped
    effect_used = effect[valid_mask]
    variance_used = variance[valid_mask]
    weights = 1.0 / variance_used
    weighted_effect = float(np.sum(weights * effect_used) / np.sum(weights))
    standard_error = float(np.sqrt(1.0 / np.sum(weights)))
    z_value = weighted_effect / standard_error
    p_value = float(2.0 * stats.norm.sf(abs(z_value)))
    return {
        "logFC": weighted_effect,
        "uncertainty_statistic": float(z_value),
        "P.Value": p_value,
        "n_peptides_used": int(valid_mask.sum()),
    }, dropped


def _aggregate_random_effect(
    *,
    site_rows: pd.DataFrame,
    min_peptides_per_site: int,
    tau2_floor: float,
) -> tuple[dict[str, float | int], int]:
    effect = site_rows.loc[:, "logFC"].to_numpy(dtype=float)
    variance, valid_mask = _derive_variance(site_rows=site_rows)
    dropped = int((~valid_mask).sum())
    n_used = int(valid_mask.sum())
    if n_used < min_peptides_per_site:
        return _nan_aggregate_result(n_used=n_used), dropped
    effect_used = effect[valid_mask]
    variance_used = variance[valid_mask]
    if n_used == 1:
        return _aggregate_fixed_effect(
            site_rows=site_rows.loc[valid_mask, :],
            min_peptides_per_site=1,
        )

    weights_fixed = 1.0 / variance_used
    mu_fixed = float(np.sum(weights_fixed * effect_used) / np.sum(weights_fixed))
    q_stat = float(np.sum(weights_fixed * (effect_used - mu_fixed) ** 2))
    c_value = float(
        np.sum(weights_fixed) - (np.sum(weights_fixed**2) / np.sum(weights_fixed))
    )
    if c_value <= 0.0:
        tau2 = float(tau2_floor)
    else:
        tau2 = max(float((q_stat - (n_used - 1)) / c_value), float(tau2_floor))
    weights_random = 1.0 / (variance_used + tau2)
    mu_random = float(np.sum(weights_random * effect_used) / np.sum(weights_random))
    standard_error = float(np.sqrt(1.0 / np.sum(weights_random)))
    z_value = mu_random / standard_error
    p_value = float(2.0 * stats.norm.sf(abs(z_value)))
    return {
        "logFC": mu_random,
        "uncertainty_statistic": float(z_value),
        "P.Value": p_value,
        "n_peptides_used": n_used,
    }, dropped


def _aggregate_stouffer(
    *,
    site_rows: pd.DataFrame,
    min_peptides_per_site: int,
    weighting: str,
) -> tuple[dict[str, float | int], int]:
    z_values = _derive_z_values(site_rows=site_rows)
    finite_mask = np.isfinite(z_values)
    if weighting == STOUFFER_WEIGHTING_INVERSE_VARIANCE:
        variance, valid_variance = _derive_variance(site_rows=site_rows)
        finite_mask = finite_mask & valid_variance
        weights = (
            np.sqrt(1.0 / variance[finite_mask])
            if np.any(finite_mask)
            else np.array([])
        )
        dropped = int((~(np.isfinite(z_values) & valid_variance)).sum())
    else:
        weights = np.ones(int(finite_mask.sum()), dtype=float)
        dropped = int((~finite_mask).sum())

    n_used = int(finite_mask.sum())
    if n_used < min_peptides_per_site or n_used == 0:
        return _nan_aggregate_result(n_used=n_used), dropped
    z_used = z_values[finite_mask]
    if weights.size != z_used.size:
        raise PhosPyInputError(
            "internal aggregation error: stouffer weight vector size mismatch"
        )
    denominator = float(np.sqrt(np.sum(weights**2)))
    if not np.isfinite(denominator) or denominator <= 0.0:
        return _nan_aggregate_result(n_used=n_used), dropped
    combined_z = float(np.sum(weights * z_used) / denominator)
    p_value = float(2.0 * stats.norm.sf(abs(combined_z)))
    log_fc_values = site_rows.loc[:, "logFC"].to_numpy(dtype=float)
    log_fc_used = log_fc_values[finite_mask]
    log_fc = float(np.mean(log_fc_used))
    return {
        "logFC": log_fc,
        "uncertainty_statistic": combined_z,
        "P.Value": p_value,
        "n_peptides_used": n_used,
    }, dropped


def _derive_variance(*, site_rows: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    log_fc = site_rows.loc[:, "logFC"].to_numpy(dtype=float)
    t_stat = site_rows.loc[:, "t"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        standard_error = np.abs(log_fc / t_stat)
        variance = standard_error**2
    valid_mask = (
        np.isfinite(variance)
        & (variance > 0.0)
        & np.isfinite(log_fc)
        & np.isfinite(t_stat)
    )
    return variance, valid_mask


def _derive_z_values(*, site_rows: pd.DataFrame) -> np.ndarray:
    t_stat = site_rows.loc[:, "t"].to_numpy(dtype=float)
    p_values = site_rows.loc[:, "P.Value"].to_numpy(dtype=float)
    log_fc = site_rows.loc[:, "logFC"].to_numpy(dtype=float)
    z = np.full(t_stat.shape, np.nan, dtype=float)
    finite_t = np.isfinite(t_stat)
    z[finite_t] = t_stat[finite_t]
    unresolved = ~finite_t
    if np.any(unresolved):
        p_unresolved = p_values[unresolved]
        sign_unresolved = np.sign(log_fc[unresolved])
        finite_p = (
            np.isfinite(p_unresolved) & (p_unresolved > 0.0) & (p_unresolved <= 1.0)
        )
        if np.any(finite_p):
            z_unresolved = np.full(p_unresolved.shape, np.nan, dtype=float)
            z_unresolved[finite_p] = sign_unresolved[finite_p] * stats.norm.isf(
                p_unresolved[finite_p] / 2.0
            )
            z[unresolved] = z_unresolved
    return z


def _nan_aggregate_result(*, n_used: int) -> dict[str, float | int]:
    return {
        "logFC": float("nan"),
        "uncertainty_statistic": float("nan"),
        "P.Value": float("nan"),
        "n_peptides_used": int(n_used),
    }
