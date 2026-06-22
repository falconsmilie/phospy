"""KSEA-style substrate-set enrichment activity score method."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from phospy.science.activities.models import (
    KSEA_ZSCORE_ACTIVITY_METHOD,
    ActivityMethodSummary,
    KinaseActivityInputs,
    KinaseActivityResult,
    KseaZScoreActivityDiagnostics,
)
from phospy.science.activities.scientific_policies import (
    build_ksea_zscore_activity_policy,
)
from phospy.science.activities.statistics import (
    benjamini_hochberg_q_values,
    two_sided_normal_p_value,
)
from phospy.science.activities.threshold_membership import (
    build_activity_threshold_membership_diagnostics,
    resolve_activity_threshold_membership_policy,
    threshold_membership_filtered_frame,
    threshold_membership_mask_array,
)

KSEA_STATUS_COMPUTED = "computed"
KSEA_STATUS_INSUFFICIENT_SUBSTRATES = "insufficient_substrates"
KSEA_STATUS_ZERO_BACKGROUND_VARIANCE = "zero_background_variance"
KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES = "no_finite_background_values"
KSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES = "no_finite_substrate_values"

KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION = "normal_approximation"
KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG = "benjamini_hochberg"


@dataclass(frozen=True, slots=True)
class KseaZScoreActivityMethod:
    """KSEA v1: unweighted substrate-set enrichment kinase activity score."""

    evidence_threshold: float
    min_substrates: int
    p_value_method: str = KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION
    adjust_p_values: bool = True
    q_value_method: str = KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG

    def run(self, inputs: KinaseActivityInputs) -> KinaseActivityResult:
        if self.p_value_method != KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION:
            raise ValueError("ksea p_value_method must be 'normal_approximation' in v1")
        aligned_pred_mat, aligned_matrix = _align_activity_inputs(
            pred_mat=inputs.pred_mat,
            phospho_matrix=inputs.phospho_matrix,
        )
        kinases = aligned_pred_mat.columns.astype(str)
        conditions = aligned_matrix.columns.astype(str)
        kinase_index = pd.Index(kinases, name="kinase")
        condition_index = pd.Index(conditions, name=aligned_matrix.columns.name)

        evidence_values = aligned_pred_mat.to_numpy(dtype=float, copy=False)
        phospho_values = aligned_matrix.to_numpy(dtype=float, copy=False)
        threshold_diagnostics = build_activity_threshold_membership_diagnostics(
            threshold_parameter="evidence_threshold",
            threshold_value=float(self.evidence_threshold),
        )
        threshold_policy = resolve_activity_threshold_membership_policy()
        membership_mask = threshold_membership_mask_array(
            evidence_values,
            threshold=float(self.evidence_threshold),
        )
        finite_phospho_mask = np.isfinite(phospho_values)

        z_scores = pd.DataFrame(
            np.nan,
            index=kinase_index,
            columns=condition_index,
            dtype=float,
        )
        p_value_matrix = pd.DataFrame(
            np.nan,
            index=kinase_index,
            columns=condition_index,
            dtype=float,
        )
        q_value_matrix = pd.DataFrame(
            np.nan,
            index=kinase_index,
            columns=condition_index,
            dtype=float,
        )
        substrate_means = pd.DataFrame(
            np.nan,
            index=kinase_index,
            columns=condition_index,
            dtype=float,
        )
        substrate_count_table = pd.DataFrame(
            0,
            index=kinase_index,
            columns=condition_index,
            dtype=int,
        )
        rows: list[dict[str, object]] = []
        counts = {
            KSEA_STATUS_COMPUTED: 0,
            KSEA_STATUS_INSUFFICIENT_SUBSTRATES: 0,
            KSEA_STATUS_ZERO_BACKGROUND_VARIANCE: 0,
            KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES: 0,
            KSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES: 0,
        }

        n_kinases = int(len(kinase_index))
        n_conditions = int(len(condition_index))
        for condition_position, condition_name in enumerate(condition_index):
            background_mask = finite_phospho_mask[:, condition_position]
            background_values = phospho_values[background_mask, condition_position]
            n_background = int(background_values.size)
            mean_background = (
                float(background_values.mean()) if n_background > 0 else np.nan
            )
            sd_background = (
                float(background_values.std(ddof=1)) if n_background >= 2 else np.nan
            )
            background_status: str | None = None
            background_reason = ""
            if n_background < 2:
                background_status = KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES
                background_reason = (
                    "background requires at least 2 finite values for sample variance"
                )
            elif not np.isfinite(mean_background) or not np.isfinite(sd_background):
                background_status = KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES
                background_reason = (
                    "background mean or standard deviation is not finite"
                )
            elif float(sd_background) == 0.0:
                background_status = KSEA_STATUS_ZERO_BACKGROUND_VARIANCE
                background_reason = "background sample standard deviation is zero"

            for kinase_position, kinase_name in enumerate(kinase_index):
                substrate_mask = membership_mask[:, kinase_position] & background_mask
                n_substrates = int(substrate_mask.sum())
                substrate_count_table.iat[kinase_position, condition_position] = (
                    n_substrates
                )
                substrate_values = phospho_values[substrate_mask, condition_position]
                if n_substrates > 0:
                    substrate_means.iat[kinase_position, condition_position] = float(
                        substrate_values.mean()
                    )

                status, reason = _resolve_status(
                    n_substrates=n_substrates,
                    min_substrates=self.min_substrates,
                    background_status=background_status,
                    background_reason=background_reason,
                )
                z_score = np.nan
                p_value = np.nan
                if status == KSEA_STATUS_COMPUTED:
                    mean_substrate = float(substrate_values.mean())
                    if not np.isfinite(mean_substrate):
                        status = KSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES
                        reason = "substrate mean is not finite"
                    else:
                        z_score = (
                            (mean_substrate - mean_background)
                            * np.sqrt(float(n_substrates))
                            / float(sd_background)
                        )
                        p_value = two_sided_normal_p_value(float(z_score))
                        z_scores.iat[kinase_position, condition_position] = z_score
                        p_value_matrix.iat[kinase_position, condition_position] = (
                            p_value
                        )

                counts[status] += 1
                rows.append(
                    {
                        "kinase": str(kinase_name),
                        "condition": str(condition_name),
                        "z_score": float(z_score) if np.isfinite(z_score) else np.nan,
                        "p_value": float(p_value) if np.isfinite(p_value) else np.nan,
                        "q_value": np.nan,
                        "n_substrates": int(n_substrates),
                        "n_background_sites": int(n_background),
                        "evidence_threshold": float(self.evidence_threshold),
                        "evidence_threshold_operator": threshold_policy.operator,
                        "evidence_threshold_description": (
                            threshold_policy.description
                        ),
                        "min_substrates": int(self.min_substrates),
                        "computability_status": status,
                        "reason": reason,
                    }
                )

        statistics_table = pd.DataFrame.from_records(
            rows,
            columns=[
                "kinase",
                "condition",
                "z_score",
                "p_value",
                "q_value",
                "n_substrates",
                "n_background_sites",
                "evidence_threshold",
                "evidence_threshold_operator",
                "evidence_threshold_description",
                "min_substrates",
                "computability_status",
                "reason",
            ],
        )
        if self.adjust_p_values:
            for condition_name in condition_index:
                condition_mask = statistics_table.loc[:, "condition"].astype(
                    str
                ) == str(condition_name)
                computed_mask = (
                    statistics_table.loc[:, "computability_status"]
                    == KSEA_STATUS_COMPUTED
                )
                selected = condition_mask & computed_mask
                if not bool(selected.any()):
                    continue
                condition_p_values = statistics_table.loc[selected, "p_value"].astype(
                    float
                )
                q_values = benjamini_hochberg_q_values(condition_p_values)
                statistics_table.loc[selected, "q_value"] = q_values.to_numpy(
                    dtype=float,
                    copy=False,
                )
                condition_rows = statistics_table.loc[selected, "kinase"].astype(str)
                for kinase_name, q_value in zip(
                    condition_rows.tolist(),
                    q_values.to_numpy(dtype=float, copy=False).tolist(),
                    strict=True,
                ):
                    q_value_matrix.at[str(kinase_name), str(condition_name)] = float(
                        q_value
                    )

        target_counts = pd.Series(
            membership_mask.sum(axis=0).astype("int64"),
            index=kinase_index.copy(),
            name="n_targets",
            dtype="int64",
        ).sort_values(ascending=False)
        # Compatibility sidecar: global evidence-membership counts after thresholding.
        thresholded_substrate_counts = target_counts.rename("n_substrates")
        target_table = _build_target_table(
            pred_mat=aligned_pred_mat,
            evidence_threshold=float(self.evidence_threshold),
        )

        summary = ActivityMethodSummary(
            kinases_evaluated=n_kinases,
            kinase_condition_pairs_evaluated=n_kinases * n_conditions,
            kinase_condition_pairs_computed=counts[KSEA_STATUS_COMPUTED],
            kinase_condition_pairs_insufficient_substrates=counts[
                KSEA_STATUS_INSUFFICIENT_SUBSTRATES
            ],
            kinase_condition_pairs_invalid_background_variance=(
                counts[KSEA_STATUS_ZERO_BACKGROUND_VARIANCE]
                + counts[KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES]
            ),
            kinase_condition_pairs_no_finite_background_values=counts[
                KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES
            ],
            kinase_condition_pairs_no_finite_substrate_values=counts[
                KSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES
            ],
        )
        diagnostics = KseaZScoreActivityDiagnostics(
            method_summary=summary,
            threshold_membership_diagnostics=threshold_diagnostics,
            statistics_table=statistics_table,
        )
        policy = build_ksea_zscore_activity_policy(
            evidence_threshold=float(self.evidence_threshold),
            min_substrates=int(self.min_substrates),
            p_value_method=str(self.p_value_method),
            adjust_p_values=bool(self.adjust_p_values),
            q_value_method=(str(self.q_value_method) if self.adjust_p_values else None),
        )

        return KinaseActivityResult._from_owned(
            weighted_activity=z_scores,
            p_value_matrix=p_value_matrix,
            q_value_matrix=q_value_matrix if self.adjust_p_values else None,
            thresholded_substrate_mean_activity=substrate_means,
            thresholded_substrate_counts=thresholded_substrate_counts,
            activity_substrate_counts=substrate_count_table,
            substrate_count_matrix=substrate_count_table,
            target_counts=target_counts,
            target_table=target_table,
            threshold_membership_diagnostics=threshold_diagnostics,
            statistics_table=statistics_table,
            method_summary=summary,
            method_diagnostics=diagnostics,
            policy_provenance=(policy,),
            activity_method=KSEA_ZSCORE_ACTIVITY_METHOD,
        )


def _resolve_status(
    *,
    n_substrates: int,
    min_substrates: int,
    background_status: str | None,
    background_reason: str,
) -> tuple[str, str]:
    if n_substrates < int(min_substrates):
        return (
            KSEA_STATUS_INSUFFICIENT_SUBSTRATES,
            f"n_substrates={n_substrates} is below min_substrates={min_substrates}",
        )
    if background_status is not None:
        return background_status, background_reason
    return KSEA_STATUS_COMPUTED, ""


def _align_activity_inputs(
    *,
    pred_mat: pd.DataFrame,
    phospho_matrix: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    common_sites = pred_mat.index.intersection(phospho_matrix.index)
    return pred_mat.loc[common_sites], phospho_matrix.loc[common_sites]


def _build_target_table(
    *,
    pred_mat: pd.DataFrame,
    evidence_threshold: float,
) -> pd.DataFrame:
    filtered = threshold_membership_filtered_frame(
        pred_mat,
        threshold=evidence_threshold,
    )
    try:
        edges = filtered.stack(future_stack=True).rename("score").reset_index()
    except TypeError:
        edges = filtered.stack().rename("score").reset_index()
    edges = edges.loc[edges["score"].notna()]
    edges.columns = ["site_id", "kinase", "score"]
    return edges.sort_values(["kinase", "score"], ascending=[True, False])


__all__ = [
    "KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION",
    "KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG",
    "KSEA_STATUS_COMPUTED",
    "KSEA_STATUS_INSUFFICIENT_SUBSTRATES",
    "KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES",
    "KSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES",
    "KSEA_STATUS_ZERO_BACKGROUND_VARIANCE",
    "KseaZScoreActivityMethod",
]
