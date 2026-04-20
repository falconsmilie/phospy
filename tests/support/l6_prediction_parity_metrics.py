from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import pandas as pd

from phospy import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    ReferencePreset,
)
from tests.support.rewrite_fixture_data import (
    build_rat_l6_dataset,
    load_l6_prediction_reference_candidate_substrates,
    load_l6_prediction_reference_combined_scores,
    load_l6_prediction_reference_predmat,
    load_l6_prediction_reference_profile_scores,
    load_l6_prediction_reference_top30,
    load_l6_prediction_reference_weights,
)


@dataclass(frozen=True, slots=True)
class TableParityMetrics:
    observed_shape: tuple[int, int]
    expected_shape: tuple[int, int]
    shared_site_count: int
    shared_kinase_count: int
    mean_abs_diff: float
    max_abs_diff: float
    mean_pearson_corr: float
    mean_spearman_corr: float


@dataclass(frozen=True, slots=True)
class CandidateParityMetrics:
    observed_rows: int
    expected_rows: int
    overlap_rows: int
    overlap_precision: float
    overlap_recall: float
    overlap_f1: float
    observed_kinase_count: int
    expected_kinase_count: int
    shared_kinase_count: int


@dataclass(frozen=True, slots=True)
class RankingParityMetrics:
    kinases_compared: int
    mean_spearman_rank_corr: float
    mean_top10_overlap: float
    mean_top20_overlap: float
    mean_top30_overlap: float
    top_rank_matches: int
    top_rank_total: int
    good_top10_count: int


@dataclass(frozen=True, slots=True)
class L6PredictionParityMetrics:
    profile: TableParityMetrics
    combined: TableParityMetrics
    weights: TableParityMetrics
    candidates: CandidateParityMetrics
    stable_ranking: RankingParityMetrics
    r_parity_ranking: RankingParityMetrics
    cross_policy_prediction_corr: float
    cross_policy_prediction_mae: float
    cross_policy_mean_top10_overlap: float
    cross_policy_mean_top20_overlap: float
    cross_policy_mean_top30_overlap: float


def _mean_column_correlation(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    method: str,
) -> float:
    common_columns = observed.columns.intersection(expected.columns)
    common_index = observed.index.intersection(expected.index)
    correlations: list[float] = []
    for column in common_columns:
        corr = observed.loc[common_index, column].corr(
            expected.loc[common_index, column],
            method=method,
        )
        if pd.notna(corr):
            correlations.append(float(corr))
    if not correlations:
        return 0.0
    return float(pd.Series(correlations).mean())


def _collect_table_parity_metrics(
    *,
    observed: pd.DataFrame,
    expected: pd.DataFrame,
) -> TableParityMetrics:
    shared_sites = observed.index.intersection(expected.index)
    shared_kinases = observed.columns.intersection(expected.columns)
    aligned_observed = observed.loc[shared_sites, shared_kinases].astype(float)
    aligned_expected = expected.loc[shared_sites, shared_kinases].astype(float)
    absolute_delta = (aligned_observed - aligned_expected).abs()
    flattened = absolute_delta.to_numpy(dtype=float).ravel()
    return TableParityMetrics(
        observed_shape=(int(observed.shape[0]), int(observed.shape[1])),
        expected_shape=(int(expected.shape[0]), int(expected.shape[1])),
        shared_site_count=int(shared_sites.size),
        shared_kinase_count=int(shared_kinases.size),
        mean_abs_diff=float(pd.Series(flattened).mean()) if flattened.size else 0.0,
        max_abs_diff=float(pd.Series(flattened).max()) if flattened.size else 0.0,
        mean_pearson_corr=_mean_column_correlation(
            aligned_observed,
            aligned_expected,
            method="pearson",
        ),
        mean_spearman_corr=_mean_column_correlation(
            aligned_observed,
            aligned_expected,
            method="spearman",
        ),
    )


def _collect_candidate_metrics(
    *,
    observed: pd.DataFrame,
    expected: pd.DataFrame,
) -> CandidateParityMetrics:
    observed_set = set(
        map(
            tuple,
            observed.loc[:, ["kinase", "site_id"]].itertuples(index=False, name=None),
        )
    )
    expected_set = set(
        map(
            tuple,
            expected.loc[:, ["kinase", "site_id"]].itertuples(index=False, name=None),
        )
    )
    overlap_rows = int(len(observed_set & expected_set))
    observed_rows = int(len(observed_set))
    expected_rows = int(len(expected_set))
    precision = float(overlap_rows / observed_rows) if observed_rows else 0.0
    recall = float(overlap_rows / expected_rows) if expected_rows else 0.0
    f1 = (
        (2.0 * precision * recall / (precision + recall))
        if (precision + recall) > 0.0
        else 0.0
    )
    observed_kinases = set(observed.loc[:, "kinase"].astype(str).tolist())
    expected_kinases = set(expected.loc[:, "kinase"].astype(str).tolist())
    return CandidateParityMetrics(
        observed_rows=observed_rows,
        expected_rows=expected_rows,
        overlap_rows=overlap_rows,
        overlap_precision=precision,
        overlap_recall=recall,
        overlap_f1=f1,
        observed_kinase_count=len(observed_kinases),
        expected_kinase_count=len(expected_kinases),
        shared_kinase_count=len(observed_kinases & expected_kinases),
    )


def _top_n_overlap(
    expected_ranked: list[str], observed_ranked: list[str], n: int
) -> float:
    expected_top = expected_ranked[:n]
    observed_top = observed_ranked[:n]
    if not expected_top:
        return 0.0
    return len(set(expected_top) & set(observed_top)) / float(len(expected_top))


def _collect_ranking_metrics(
    *,
    observed_pred_mat: pd.DataFrame,
    expected_pred_mat: pd.DataFrame,
    expected_top30: pd.DataFrame,
) -> RankingParityMetrics:
    expected_ranked_map = {
        str(kinase): group.sort_values("rank").loc[:, "site_id"].astype(str).tolist()
        for kinase, group in expected_top30.groupby("kinase", sort=False)
    }
    kinases = sorted(
        set(observed_pred_mat.columns.astype(str))
        & set(expected_pred_mat.columns.astype(str))
        & set(expected_ranked_map)
    )
    top10: list[float] = []
    top20: list[float] = []
    top30: list[float] = []
    spearman: list[float] = []
    top_rank_matches = 0

    for kinase in kinases:
        observed_series = observed_pred_mat.loc[:, kinase].astype(float).dropna()
        if observed_series.empty:
            continue
        observed_ranked = (
            observed_pred_mat.loc[:, kinase]
            .astype(float)
            .sort_values(ascending=False)
            .index.astype(str)
            .tolist()
        )
        expected_ranked = expected_ranked_map[kinase]
        top10.append(_top_n_overlap(expected_ranked, observed_ranked, 10))
        top20.append(_top_n_overlap(expected_ranked, observed_ranked, 20))
        top30.append(_top_n_overlap(expected_ranked, observed_ranked, 30))
        if (
            observed_ranked
            and expected_ranked
            and observed_ranked[0] == expected_ranked[0]
        ):
            top_rank_matches += 1
        expected_series = expected_pred_mat.loc[observed_series.index, kinase].astype(
            float
        )
        corr = expected_series.rank(ascending=False, method="average").corr(
            observed_series.rank(ascending=False, method="average"),
            method="spearman",
        )
        if pd.notna(corr):
            spearman.append(float(corr))

    kinases_compared = len(spearman)
    return RankingParityMetrics(
        kinases_compared=kinases_compared,
        mean_spearman_rank_corr=float(pd.Series(spearman).mean()) if spearman else 0.0,
        mean_top10_overlap=float(pd.Series(top10).mean()) if top10 else 0.0,
        mean_top20_overlap=float(pd.Series(top20).mean()) if top20 else 0.0,
        mean_top30_overlap=float(pd.Series(top30).mean()) if top30 else 0.0,
        top_rank_matches=int(top_rank_matches),
        top_rank_total=int(kinases_compared),
        good_top10_count=int(sum(overlap >= 0.70 for overlap in top10)),
    )


def _ranked_sites_by_kinase(pred_mat: pd.DataFrame) -> dict[str, list[str]]:
    ranked: dict[str, list[str]] = {}
    for kinase in pred_mat.columns.astype(str):
        ordered = (
            pred_mat.loc[:, kinase]
            .astype(float)
            .sort_values(ascending=False, kind="mergesort")
            .index.astype(str)
            .tolist()
        )
        ranked[kinase] = ordered
    return ranked


def _mean_top_n_overlap(
    ranked_left: dict[str, list[str]],
    ranked_right: dict[str, list[str]],
    *,
    top_n: int,
) -> float:
    shared_kinases = sorted(set(ranked_left) & set(ranked_right))
    overlaps: list[float] = []
    for kinase in shared_kinases:
        left_top = ranked_left[kinase][:top_n]
        right_top = ranked_right[kinase][:top_n]
        if not left_top:
            continue
        overlaps.append(len(set(left_top) & set(right_top)) / float(len(left_top)))
    if not overlaps:
        return 0.0
    return float(pd.Series(overlaps).mean())


@lru_cache(maxsize=1)
def _run_l6_workflow(adaptive_policy: str):
    dataset = build_rat_l6_dataset(n_sites=None)
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
            ),
            prediction_config=KinasePredictionConfig(
                top_k=30,
                ensemble_size=10,
                mode="adaptive_ensemble",
                adaptive_policy=adaptive_policy,
                n_iterations=5,
                random_state=1,
            ),
            activity_config=None,
        )
    )


@lru_cache(maxsize=1)
def collect_l6_prediction_parity_metrics() -> L6PredictionParityMetrics:
    stable_result = _run_l6_workflow("stable")
    r_parity_result = _run_l6_workflow("r_parity")

    observed_profile = stable_result.scoring_result.profile_scores
    expected_profile = load_l6_prediction_reference_profile_scores()
    observed_combined = stable_result.scoring_result.combined_scores
    if observed_combined is None:
        raise AssertionError("expected combined scores to be present for L6 parity")
    expected_combined = load_l6_prediction_reference_combined_scores()
    observed_weights = stable_result.scoring_result.weights
    if observed_weights is None:
        raise AssertionError("expected weight table to be present for L6 parity")
    expected_weights = load_l6_prediction_reference_weights()
    expected_candidates = load_l6_prediction_reference_candidate_substrates()
    expected_pred_mat = load_l6_prediction_reference_predmat()
    expected_top30 = load_l6_prediction_reference_top30()

    observed_candidates = stable_result.prediction_result.substrate_list
    if observed_candidates is None:
        raise AssertionError("expected substrate_list to be present for L6 parity")
    observed_candidates_frame = observed_candidates.loc[
        :, ["kinase", "substrate_site"]
    ].rename(columns={"substrate_site": "site_id"})

    stable_pred = stable_result.prediction_result.pred_mat
    r_parity_pred = r_parity_result.prediction_result.pred_mat
    stable_long = (
        stable_pred.stack(dropna=False, future_stack=False)
        .rename("score_stable")
        .reset_index()
        .rename(
            columns={
                "level_1": "kinase",
                stable_pred.index.name or "level_0": "site_id",
            }
        )
    )
    r_parity_long = (
        r_parity_pred.stack(dropna=False, future_stack=False)
        .rename("score_r_parity")
        .reset_index()
        .rename(
            columns={
                "level_1": "kinase",
                r_parity_pred.index.name or "level_0": "site_id",
            }
        )
    )
    merged = stable_long.merge(
        r_parity_long,
        on=["site_id", "kinase"],
        how="inner",
        validate="one_to_one",
    )
    merged_non_null = merged.dropna(subset=["score_stable", "score_r_parity"])
    cross_policy_delta = (
        merged_non_null.loc[:, "score_stable"]
        - merged_non_null.loc[:, "score_r_parity"]
    ).abs()
    stable_ranked = _ranked_sites_by_kinase(stable_pred)
    r_parity_ranked = _ranked_sites_by_kinase(r_parity_pred)

    return L6PredictionParityMetrics(
        profile=_collect_table_parity_metrics(
            observed=observed_profile,
            expected=expected_profile,
        ),
        combined=_collect_table_parity_metrics(
            observed=observed_combined,
            expected=expected_combined,
        ),
        weights=_collect_table_parity_metrics(
            observed=observed_weights,
            expected=expected_weights,
        ),
        candidates=_collect_candidate_metrics(
            observed=observed_candidates_frame,
            expected=expected_candidates,
        ),
        stable_ranking=_collect_ranking_metrics(
            observed_pred_mat=stable_pred,
            expected_pred_mat=expected_pred_mat,
            expected_top30=expected_top30,
        ),
        r_parity_ranking=_collect_ranking_metrics(
            observed_pred_mat=r_parity_pred,
            expected_pred_mat=expected_pred_mat,
            expected_top30=expected_top30,
        ),
        cross_policy_prediction_corr=float(
            merged_non_null.loc[:, "score_stable"].corr(
                merged_non_null.loc[:, "score_r_parity"]
            )
        )
        if not merged_non_null.empty
        else 0.0,
        cross_policy_prediction_mae=float(cross_policy_delta.mean())
        if not cross_policy_delta.empty
        else 0.0,
        cross_policy_mean_top10_overlap=_mean_top_n_overlap(
            stable_ranked,
            r_parity_ranked,
            top_n=10,
        ),
        cross_policy_mean_top20_overlap=_mean_top_n_overlap(
            stable_ranked,
            r_parity_ranked,
            top_n=20,
        ),
        cross_policy_mean_top30_overlap=_mean_top_n_overlap(
            stable_ranked,
            r_parity_ranked,
            top_n=30,
        ),
    )
