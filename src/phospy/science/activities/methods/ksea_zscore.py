"""KSEA-style substrate-set enrichment activity score method."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd

from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.activities.membership import (
    ACTIVITY_MEMBERSHIP_SOURCE_INCOMPLETE,
    ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN,
    ActivityMembershipSelection,
    KseaMembershipInferentialDecision,
    fingerprint_ksea_tested_quantitative_matrix,
    selected_substrate_universe_from_prediction_matrix,
)
from phospy.science.activities.method_contracts import (
    ksea_zscore_activity_input_contract,
)
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
    ActivityThresholdMembershipDiagnostics,
    ActivityThresholdMembershipPolicy,
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
KSEA_INFERENTIAL_STATUS_ORDINARY_P_Q_AVAILABLE = "ordinary_p_q_available"
KSEA_INFERENTIAL_STATUS_P_Q_UNAVAILABLE = "ordinary_p_q_unavailable"


@dataclass(frozen=True, slots=True)
class _KseaPreparedInputs:
    aligned_pred_mat: pd.DataFrame
    aligned_matrix: pd.DataFrame
    kinase_index: pd.Index
    profile_index: pd.Index
    evidence_threshold: float
    min_substrates: int
    phospho_values: npt.NDArray[np.float64]
    membership_mask: npt.NDArray[np.bool_]
    finite_phospho_mask: npt.NDArray[np.bool_]
    threshold_diagnostics: ActivityThresholdMembershipDiagnostics
    threshold_policy: ActivityThresholdMembershipPolicy
    membership_selection: ActivityMembershipSelection
    inferential_decision: KseaMembershipInferentialDecision
    condition_ids_by_profile: dict[str, str]


@dataclass(frozen=True, slots=True)
class _KseaProfileBackground:
    mask: npt.NDArray[np.bool_]
    values: npt.NDArray[np.float64]
    count: int
    mean: float
    sd: float
    status: str | None
    reason: str


@dataclass(frozen=True, slots=True)
class _KseaPairScore:
    substrate_values: npt.NDArray[np.float64]
    n_substrates: int
    status: str
    reason: str
    z_score: float
    p_value: float


@dataclass(frozen=True, slots=True)
class _KseaScoredOutputs:
    z_scores: pd.DataFrame
    p_value_matrix: pd.DataFrame | None
    q_value_matrix: pd.DataFrame | None
    substrate_means: pd.DataFrame
    substrate_count_table: pd.DataFrame
    statistics_table: pd.DataFrame
    status_counts: dict[str, int]


@dataclass(frozen=True, slots=True)
class _KseaMutableScoringTables:
    z_scores: pd.DataFrame
    p_value_matrix: pd.DataFrame | None
    q_value_matrix: pd.DataFrame | None
    substrate_means: pd.DataFrame
    substrate_count_table: pd.DataFrame
    rows: list[dict[str, object]]
    status_counts: dict[str, int]


@dataclass(frozen=True, slots=True)
class KseaZScoreActivityMethod:
    """KSEA v1: unweighted substrate-set enrichment kinase activity score."""

    evidence_threshold: float
    min_substrates: int
    p_value_method: str = KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION
    adjust_p_values: bool = True
    q_value_method: str = KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG

    def run(self, inputs: KinaseActivityInputs) -> KinaseActivityResult:
        _validate_ksea_method_configuration(self)
        prepared = _prepare_ksea_inputs(method=self, inputs=inputs)
        scored = _score_ksea_profiles(method=self, prepared=prepared)
        corrected = _apply_ksea_multiple_testing(
            method=self,
            prepared=prepared,
            scored=scored,
        )
        return _assemble_ksea_result(
            method=self,
            inputs=inputs,
            prepared=prepared,
            scored=corrected,
        )


def _validate_ksea_method_configuration(method: KseaZScoreActivityMethod) -> None:
    if method.p_value_method != KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION:
        raise ValueError("ksea p_value_method must be 'normal_approximation' in v1")


def _prepare_ksea_inputs(
    *,
    method: KseaZScoreActivityMethod,
    inputs: KinaseActivityInputs,
) -> _KseaPreparedInputs:
    from phospy.science.quantitative_method_contracts import (
        resolve_activity_input_contract,
    )

    if inputs.activity_input is None:
        raise WorkflowBoundaryError(
            "KSEA-style activity requires typed ActivityInputMatrix semantics"
        )
    resolve_activity_input_contract(
        activity_input=inputs.activity_input,
        contract=ksea_zscore_activity_input_contract(),
        context="KSEA-style activity input",
    )
    aligned_pred_mat, aligned_matrix = _align_activity_inputs(
        pred_mat=inputs.pred_mat,
        phospho_matrix=inputs.phospho_matrix,
    )
    kinase_index = pd.Index(aligned_pred_mat.columns.astype(str), name="kinase")
    profile_index = pd.Index(
        aligned_matrix.columns.astype(str),
        name=aligned_matrix.columns.name,
    )
    evidence_values: npt.NDArray[np.float64] = aligned_pred_mat.to_numpy(
        dtype=float,
        copy=False,
    )
    phospho_values: npt.NDArray[np.float64] = aligned_matrix.to_numpy(
        dtype=float,
        copy=False,
    )
    threshold_diagnostics = build_activity_threshold_membership_diagnostics(
        threshold_parameter="evidence_threshold",
        threshold_value=float(method.evidence_threshold),
    )
    threshold_policy = resolve_activity_threshold_membership_policy()
    membership_mask: npt.NDArray[np.bool_] = threshold_membership_mask_array(
        evidence_values,
        threshold=float(method.evidence_threshold),
    )
    membership_selection = _require_membership_selection(inputs)
    _validate_ksea_membership_boundary(
        membership_selection=membership_selection,
        aligned_pred_mat=aligned_pred_mat,
        aligned_matrix=aligned_matrix,
        evidence_threshold=float(method.evidence_threshold),
    )
    finite_phospho_mask: npt.NDArray[np.bool_] = np.isfinite(phospho_values)
    return _KseaPreparedInputs(
        aligned_pred_mat=aligned_pred_mat,
        aligned_matrix=aligned_matrix,
        kinase_index=kinase_index,
        profile_index=profile_index,
        evidence_threshold=float(method.evidence_threshold),
        min_substrates=int(method.min_substrates),
        phospho_values=phospho_values,
        membership_mask=membership_mask,
        finite_phospho_mask=finite_phospho_mask,
        threshold_diagnostics=threshold_diagnostics,
        threshold_policy=threshold_policy,
        membership_selection=membership_selection,
        inferential_decision=membership_selection.inferential_decision,
        condition_ids_by_profile=_condition_ids_by_profile(
            profile_ids=profile_index,
            condition_ids=inputs.profile_metadata.condition_ids,
            include_condition=inputs.input_semantics.has_real_condition_contract,
        ),
    )


def _score_ksea_profiles(
    *,
    method: KseaZScoreActivityMethod,
    prepared: _KseaPreparedInputs,
) -> _KseaScoredOutputs:
    tables = _initialise_ksea_scoring_tables(method=method, prepared=prepared)
    for profile_position, profile_id in enumerate(prepared.profile_index):
        background = _profile_background(
            prepared=prepared,
            profile_position=profile_position,
        )
        for kinase_position, kinase_name in enumerate(prepared.kinase_index):
            score = _score_ksea_pair(
                method=method,
                prepared=prepared,
                background=background,
                kinase_position=kinase_position,
                profile_position=profile_position,
            )
            _record_ksea_pair_score(
                tables=tables,
                prepared=prepared,
                background=background,
                score=score,
                kinase_position=kinase_position,
                kinase_name=str(kinase_name),
                profile_position=profile_position,
                profile_id=str(profile_id),
            )
    return _KseaScoredOutputs(
        z_scores=tables.z_scores,
        p_value_matrix=tables.p_value_matrix,
        q_value_matrix=tables.q_value_matrix,
        substrate_means=tables.substrate_means,
        substrate_count_table=tables.substrate_count_table,
        statistics_table=pd.DataFrame.from_records(
            tables.rows,
            columns=_ksea_statistics_columns(prepared),
        ),
        status_counts=tables.status_counts,
    )


def _initialise_ksea_scoring_tables(
    *,
    method: KseaZScoreActivityMethod,
    prepared: _KseaPreparedInputs,
) -> _KseaMutableScoringTables:
    z_scores = pd.DataFrame(
        np.nan,
        index=prepared.kinase_index,
        columns=prepared.profile_index,
        dtype=float,
    )
    p_value_matrix = (
        pd.DataFrame(
            np.nan,
            index=prepared.kinase_index,
            columns=prepared.profile_index,
            dtype=float,
        )
        if prepared.inferential_decision.ordinary_p_q_available
        and method.adjust_p_values
        else None
    )
    q_value_matrix = (
        pd.DataFrame(
            np.nan,
            index=prepared.kinase_index,
            columns=prepared.profile_index,
            dtype=float,
        )
        if prepared.inferential_decision.ordinary_p_q_available
        else None
    )
    return _KseaMutableScoringTables(
        z_scores=z_scores,
        p_value_matrix=p_value_matrix,
        q_value_matrix=q_value_matrix,
        substrate_means=pd.DataFrame(
            np.nan,
            index=prepared.kinase_index,
            columns=prepared.profile_index,
            dtype=float,
        ),
        substrate_count_table=pd.DataFrame(
            0,
            index=prepared.kinase_index,
            columns=prepared.profile_index,
            dtype=int,
        ),
        rows=[],
        status_counts={
            KSEA_STATUS_COMPUTED: 0,
            KSEA_STATUS_INSUFFICIENT_SUBSTRATES: 0,
            KSEA_STATUS_ZERO_BACKGROUND_VARIANCE: 0,
            KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES: 0,
            KSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES: 0,
        },
    )


def _profile_background(
    *,
    prepared: _KseaPreparedInputs,
    profile_position: int,
) -> _KseaProfileBackground:
    background_mask = prepared.finite_phospho_mask[:, profile_position]
    background_values = prepared.phospho_values[background_mask, profile_position]
    n_background = int(background_values.size)
    mean_background = float(background_values.mean()) if n_background > 0 else np.nan
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
        background_reason = "background mean or standard deviation is not finite"
    elif float(sd_background) == 0.0:
        background_status = KSEA_STATUS_ZERO_BACKGROUND_VARIANCE
        background_reason = "background sample standard deviation is zero"
    return _KseaProfileBackground(
        mask=background_mask,
        values=background_values,
        count=n_background,
        mean=mean_background,
        sd=sd_background,
        status=background_status,
        reason=background_reason,
    )


def _score_ksea_pair(
    *,
    method: KseaZScoreActivityMethod,
    prepared: _KseaPreparedInputs,
    background: _KseaProfileBackground,
    kinase_position: int,
    profile_position: int,
) -> _KseaPairScore:
    substrate_mask = prepared.membership_mask[:, kinase_position] & background.mask
    substrate_values = prepared.phospho_values[substrate_mask, profile_position]
    n_substrates = int(substrate_mask.sum())
    status, reason = _resolve_status(
        n_substrates=n_substrates,
        min_substrates=method.min_substrates,
        background_status=background.status,
        background_reason=background.reason,
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
                (mean_substrate - background.mean)
                * np.sqrt(float(n_substrates))
                / float(background.sd)
            )
            if prepared.inferential_decision.ordinary_p_q_available:
                p_value = two_sided_normal_p_value(float(z_score))
    return _KseaPairScore(
        substrate_values=substrate_values,
        n_substrates=n_substrates,
        status=status,
        reason=reason,
        z_score=float(z_score) if np.isfinite(z_score) else np.nan,
        p_value=float(p_value) if np.isfinite(p_value) else np.nan,
    )


def _record_ksea_pair_score(
    *,
    tables: _KseaMutableScoringTables,
    prepared: _KseaPreparedInputs,
    background: _KseaProfileBackground,
    score: _KseaPairScore,
    kinase_position: int,
    kinase_name: str,
    profile_position: int,
    profile_id: str,
) -> None:
    tables.substrate_count_table.iat[kinase_position, profile_position] = (
        score.n_substrates
    )
    if score.n_substrates > 0:
        tables.substrate_means.iat[kinase_position, profile_position] = float(
            score.substrate_values.mean()
        )
    if score.status == KSEA_STATUS_COMPUTED:
        tables.z_scores.iat[kinase_position, profile_position] = score.z_score
        if tables.p_value_matrix is not None:
            tables.p_value_matrix.iat[kinase_position, profile_position] = score.p_value
    tables.status_counts[score.status] += 1
    tables.rows.append(
        _ksea_statistics_row(
            prepared=prepared,
            kinase_name=kinase_name,
            profile_id=profile_id,
            score=score,
            n_background=background.count,
        )
    )


def _ksea_statistics_row(
    *,
    prepared: _KseaPreparedInputs,
    kinase_name: str,
    profile_id: str,
    score: _KseaPairScore,
    n_background: int,
) -> dict[str, object]:
    row: dict[str, object] = {
        "kinase": kinase_name,
        "profile_id": profile_id,
        "z_score": score.z_score,
        "p_value": score.p_value,
        "q_value": np.nan,
        "n_substrates": int(score.n_substrates),
        "n_background_sites": int(n_background),
        "evidence_threshold": prepared.evidence_threshold,
        "evidence_threshold_operator": prepared.threshold_policy.operator,
        "evidence_threshold_description": prepared.threshold_policy.description,
        "min_substrates": prepared.min_substrates,
        "computability_status": score.status,
        "reason": score.reason,
        "inferential_eligible": bool(
            prepared.inferential_decision.ordinary_p_q_available
        ),
        "inferential_status": prepared.inferential_decision.status,
        "inferential_reason": prepared.inferential_decision.reason,
        "membership_source_category": prepared.membership_selection.source_category,
        "membership_selection_method": prepared.membership_selection.selection_method,
    }
    condition_id = prepared.condition_ids_by_profile.get(profile_id)
    if condition_id is not None:
        row["condition"] = condition_id
    return row


def _ksea_statistics_columns(prepared: _KseaPreparedInputs) -> list[str]:
    statistics_columns = [
        "kinase",
        "profile_id",
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
        "inferential_eligible",
        "inferential_status",
        "inferential_reason",
        "membership_source_category",
        "membership_selection_method",
    ]
    if prepared.condition_ids_by_profile:
        statistics_columns.insert(2, "condition")
    return statistics_columns


def _apply_ksea_multiple_testing(
    *,
    method: KseaZScoreActivityMethod,
    prepared: _KseaPreparedInputs,
    scored: _KseaScoredOutputs,
) -> _KseaScoredOutputs:
    if not (
        method.adjust_p_values and prepared.inferential_decision.ordinary_p_q_available
    ):
        return scored
    for profile_id in prepared.profile_index:
        profile_mask = scored.statistics_table.loc[:, "profile_id"].astype(str) == str(
            profile_id
        )
        computed_mask = (
            scored.statistics_table.loc[:, "computability_status"]
            == KSEA_STATUS_COMPUTED
        )
        selected = profile_mask & computed_mask
        if not bool(selected.any()):
            continue
        profile_p_values = scored.statistics_table.loc[selected, "p_value"].astype(
            float
        )
        finite_p_values = np.isfinite(
            profile_p_values.to_numpy(dtype=float, copy=False)
        )
        if not bool(finite_p_values.all()):
            selected = selected & scored.statistics_table.loc[:, "p_value"].notna()
            profile_p_values = scored.statistics_table.loc[selected, "p_value"].astype(
                float
            )
        if profile_p_values.empty:
            continue
        q_values = benjamini_hochberg_q_values(profile_p_values)
        scored.statistics_table.loc[selected, "q_value"] = q_values.to_numpy(
            dtype=float,
            copy=False,
        )
        profile_rows = scored.statistics_table.loc[selected, "kinase"].astype(str)
        for kinase_name, q_value in zip(
            profile_rows.tolist(),
            q_values.to_numpy(dtype=float, copy=False).tolist(),
            strict=True,
        ):
            if scored.q_value_matrix is not None:
                scored.q_value_matrix.at[str(kinase_name), str(profile_id)] = float(
                    q_value
                )
    return scored


def _assemble_ksea_result(
    *,
    method: KseaZScoreActivityMethod,
    inputs: KinaseActivityInputs,
    prepared: _KseaPreparedInputs,
    scored: _KseaScoredOutputs,
) -> KinaseActivityResult:
    target_counts = pd.Series(
        prepared.membership_mask.sum(axis=0).astype("int64"),
        index=prepared.kinase_index.copy(),
        name="n_targets",
        dtype="int64",
    ).sort_values(ascending=False)
    # Compatibility sidecar: global evidence-membership counts after thresholding.
    thresholded_substrate_counts = target_counts.rename("n_substrates")
    target_table = _build_target_table(
        pred_mat=prepared.aligned_pred_mat,
        evidence_threshold=float(method.evidence_threshold),
    )
    n_kinases = int(len(prepared.kinase_index))
    n_profiles = int(len(prepared.profile_index))
    summary = ActivityMethodSummary(
        kinases_evaluated=n_kinases,
        kinase_profile_pairs_evaluated=n_kinases * n_profiles,
        kinase_profile_pairs_computed=scored.status_counts[KSEA_STATUS_COMPUTED],
        kinase_profile_pairs_insufficient_substrates=scored.status_counts[
            KSEA_STATUS_INSUFFICIENT_SUBSTRATES
        ],
        kinase_profile_pairs_invalid_background_variance=(
            scored.status_counts[KSEA_STATUS_ZERO_BACKGROUND_VARIANCE]
            + scored.status_counts[KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES]
        ),
        kinase_profile_pairs_no_finite_background_values=scored.status_counts[
            KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES
        ],
        kinase_profile_pairs_no_finite_substrate_values=scored.status_counts[
            KSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES
        ],
    )
    diagnostics = KseaZScoreActivityDiagnostics(
        method_summary=summary,
        threshold_membership_diagnostics=prepared.threshold_diagnostics,
        statistics_table=scored.statistics_table,
    )
    policy = build_ksea_zscore_activity_policy(
        evidence_threshold=float(method.evidence_threshold),
        min_substrates=int(method.min_substrates),
        p_value_method=str(method.p_value_method),
        adjust_p_values=bool(method.adjust_p_values),
        q_value_method=(str(method.q_value_method) if method.adjust_p_values else None),
        membership_inferential_eligible=bool(
            prepared.inferential_decision.ordinary_p_q_available
        ),
    )

    return KinaseActivityResult.from_trusted_owned(
        weighted_activity=scored.z_scores,
        p_value_matrix=scored.p_value_matrix,
        q_value_matrix=scored.q_value_matrix,
        thresholded_substrate_mean_activity=scored.substrate_means,
        thresholded_substrate_counts=thresholded_substrate_counts,
        activity_substrate_counts=scored.substrate_count_table,
        substrate_count_matrix=scored.substrate_count_table,
        target_counts=target_counts,
        target_table=target_table,
        threshold_membership_diagnostics=prepared.threshold_diagnostics,
        statistics_table=scored.statistics_table,
        method_summary=summary,
        method_diagnostics=diagnostics,
        policy_provenance=(policy,),
        activity_method=KSEA_ZSCORE_ACTIVITY_METHOD,
        input_semantics=inputs.input_semantics,
        profile_metadata=inputs.profile_metadata,
        membership_selection=prepared.membership_selection,
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
    return pred_mat.reindex(index=phospho_matrix.index.copy()), phospho_matrix


def _require_membership_selection(
    inputs: KinaseActivityInputs,
) -> ActivityMembershipSelection:
    membership_selection = inputs.membership_selection
    if membership_selection is None:
        raise WorkflowBoundaryError(
            "KSEA-style activity requires typed membership_selection provenance"
        )
    return membership_selection


def _validate_ksea_membership_boundary(
    *,
    membership_selection: ActivityMembershipSelection,
    aligned_pred_mat: pd.DataFrame,
    aligned_matrix: pd.DataFrame,
    evidence_threshold: float,
) -> None:
    tested_fingerprint = membership_selection.tested_quantitative_matrix_fingerprint
    if tested_fingerprint is not None:
        actual_fingerprint = fingerprint_ksea_tested_quantitative_matrix(aligned_matrix)
        if tested_fingerprint != actual_fingerprint:
            raise WorkflowBoundaryError(
                "KSEA membership provenance tested_quantitative_matrix_fingerprint "
                "does not match the actual KSEA background phospho matrix; "
                "next_action=construct ActivityMembershipSelection with a "
                "fingerprint of the exact phospho_matrix passed to "
                "KseaZScoreActivityMethod.run"
            )
    expected_kinases = tuple(str(value) for value in aligned_pred_mat.columns.tolist())
    if membership_selection.selected_kinase_universe != expected_kinases:
        raise WorkflowBoundaryError(
            "KSEA membership provenance selected_kinase_universe does not match "
            "the effective prediction-matrix kinase columns; "
            f"expected={expected_kinases!r}, got="
            f"{membership_selection.selected_kinase_universe!r}; "
            "next_action=rebuild membership provenance from the prediction matrix "
            "used for KSEA"
        )
    if membership_selection.source_category in {
        ACTIVITY_MEMBERSHIP_SOURCE_UNKNOWN,
        ACTIVITY_MEMBERSHIP_SOURCE_INCOMPLETE,
    }:
        return
    expected_substrates = selected_substrate_universe_from_prediction_matrix(
        aligned_pred_mat,
        threshold=float(evidence_threshold),
    )
    if membership_selection.selected_substrate_universe != expected_substrates:
        raise WorkflowBoundaryError(
            "KSEA membership provenance selected_substrate_universe does not "
            "match thresholded membership after alignment to the KSEA background; "
            f"expected={expected_substrates!r}, got="
            f"{membership_selection.selected_substrate_universe!r}; "
            "next_action=rebuild membership provenance from the aligned KSEA "
            "prediction membership matrix"
        )


def _condition_ids_by_profile(
    *,
    profile_ids: pd.Index,
    condition_ids: tuple[str, ...],
    include_condition: bool,
) -> dict[str, str]:
    if not include_condition:
        return {}
    return {
        str(profile_id): str(condition_id)
        for profile_id, condition_id in zip(profile_ids, condition_ids, strict=True)
    }


def _build_target_table(
    *,
    pred_mat: pd.DataFrame,
    evidence_threshold: float,
) -> pd.DataFrame:
    filtered = threshold_membership_filtered_frame(
        pred_mat,
        threshold=evidence_threshold,
    )
    stacked = _stack_membership_scores(filtered)
    stacked.name = "score"
    edges = stacked.reset_index()
    edges = edges.loc[edges["score"].notna()]
    edges.columns = ["site_id", "kinase", "score"]
    return edges.sort_values(["kinase", "score"], ascending=[True, False])


def _stack_membership_scores(filtered: pd.DataFrame) -> pd.Series:
    try:
        stacked = filtered.stack(future_stack=True)
    except TypeError:
        stacked = filtered.stack()
    if not isinstance(stacked, pd.Series):
        raise WorkflowBoundaryError(
            "KSEA target-table assembly expected pandas stack() to produce a "
            "Series for a single-level column membership frame"
        )
    return stacked


__all__ = [
    "KSEA_P_VALUE_METHOD_NORMAL_APPROXIMATION",
    "KSEA_Q_VALUE_METHOD_BENJAMINI_HOCHBERG",
    "KSEA_INFERENTIAL_STATUS_ORDINARY_P_Q_AVAILABLE",
    "KSEA_INFERENTIAL_STATUS_P_Q_UNAVAILABLE",
    "KSEA_STATUS_COMPUTED",
    "KSEA_STATUS_INSUFFICIENT_SUBSTRATES",
    "KSEA_STATUS_NO_FINITE_BACKGROUND_VALUES",
    "KSEA_STATUS_NO_FINITE_SUBSTRATE_VALUES",
    "KSEA_STATUS_ZERO_BACKGROUND_VARIANCE",
    "KseaZScoreActivityMethod",
]
