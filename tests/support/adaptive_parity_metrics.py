from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import pandas as pd

from phospy.advanced import KinasePredictionConfig
from phospy.science.prediction.candidates import build_candidate_substrate_list
from phospy.science.prediction.execution import run_adaptive_ensemble_prediction
from tests.support.rewrite_fixture_data import (
    load_adaptive_sampling_edge_rank_weighted_fusion_scores,
    load_adaptive_sampling_edge_trace_candidates,
    load_adaptive_sampling_edge_trace_predictions,
    load_adaptive_sampling_edge_trace_top,
)

ADAPTIVE_PARITY_TOP_K = 4
ADAPTIVE_PARITY_SCORE_THRESHOLD = 0.8
ADAPTIVE_PARITY_INCLUSION = 1
ADAPTIVE_PARITY_ENSEMBLE_SIZE = 1
ADAPTIVE_PARITY_N_ITERATIONS = 2
ADAPTIVE_PARITY_RANDOM_STATE = 18
ADAPTIVE_POLICIES = ("stable", "r_parity")
TOP_OVERLAP_LEVELS = (10, 20, 30)
POLICY_DISPLAY_LABELS = {
    "stable": "stable (default)",
    "r_parity": "r_parity",
}


@dataclass(frozen=True, slots=True)
class AdaptivePolicyLaneMetrics:
    adaptive_policy: str
    policy_label: str
    prediction_shape: tuple[int, int]
    kinases_compared: int
    candidate_count: int
    candidate_kinase_count: int
    selected_trace_candidate_count: int
    donor_prediction_rows: int
    donor_prediction_corr: float
    donor_prediction_mae: float
    donor_prediction_max_abs_diff: float
    donor_top_rank_matches: int
    donor_top_rank_total: int
    donor_top_prob_mae: float
    donor_top_set_overlap_matches: int
    donor_top_set_overlap_total: int
    donor_mean_top10_overlap: float
    donor_mean_top20_overlap: float
    donor_mean_top30_overlap: float
    _observed_prediction_frame: pd.DataFrame
    _observed_top_frame: pd.DataFrame


@dataclass(frozen=True, slots=True)
class AdaptivePolicyComparisonMetrics:
    stable: AdaptivePolicyLaneMetrics
    r_parity: AdaptivePolicyLaneMetrics
    cross_policy_prediction_corr: float
    cross_policy_prediction_mae: float
    cross_policy_prediction_max_abs_diff: float
    cross_policy_top_rank_matches: int
    cross_policy_top_rank_total: int
    cross_policy_top_set_overlap_matches: int
    cross_policy_top_set_overlap_total: int
    cross_policy_mean_top10_overlap: float
    cross_policy_mean_top20_overlap: float
    cross_policy_mean_top30_overlap: float


@lru_cache(maxsize=1)
def collect_adaptive_policy_comparison_metrics() -> AdaptivePolicyComparisonMetrics:
    stable = _collect_lane_metrics(adaptive_policy="stable")
    r_parity = _collect_lane_metrics(adaptive_policy="r_parity")

    cross_prediction = stable._observed_prediction_frame.merge(
        r_parity._observed_prediction_frame,
        on=["kinase", "site"],
        suffixes=("_stable", "_r_parity"),
        validate="one_to_one",
    )
    prediction_delta = (
        cross_prediction.loc[:, "prob_class_1_r_parity"]
        - cross_prediction.loc[:, "prob_class_1_stable"]
    ).abs()

    cross_top_rank = stable._observed_top_frame.merge(
        r_parity._observed_top_frame,
        on=["kinase", "rank"],
        suffixes=("_stable", "_r_parity"),
        validate="one_to_one",
    )
    cross_top_rank_matches = int(
        (
            cross_top_rank.loc[:, "site_stable"]
            == cross_top_rank.loc[:, "site_r_parity"]
        ).sum()
    )

    cross_top_set_overlap_matches, cross_top_set_overlap_total = _top_set_overlap(
        stable_top=stable._observed_top_frame,
        r_parity_top=r_parity._observed_top_frame,
    )
    stable_ranked = _ranked_sites_by_kinase(stable._observed_prediction_frame)
    r_parity_ranked = _ranked_sites_by_kinase(r_parity._observed_prediction_frame)

    return AdaptivePolicyComparisonMetrics(
        stable=stable,
        r_parity=r_parity,
        cross_policy_prediction_corr=float(
            cross_prediction.loc[:, "prob_class_1_stable"].corr(
                cross_prediction.loc[:, "prob_class_1_r_parity"]
            )
        ),
        cross_policy_prediction_mae=float(prediction_delta.mean()),
        cross_policy_prediction_max_abs_diff=float(prediction_delta.max()),
        cross_policy_top_rank_matches=cross_top_rank_matches,
        cross_policy_top_rank_total=int(cross_top_rank.shape[0]),
        cross_policy_top_set_overlap_matches=cross_top_set_overlap_matches,
        cross_policy_top_set_overlap_total=cross_top_set_overlap_total,
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


def _collect_lane_metrics(*, adaptive_policy: str) -> AdaptivePolicyLaneMetrics:
    rank_weighted_fusion_scores = (
        load_adaptive_sampling_edge_rank_weighted_fusion_scores()
    )
    candidate_substrates = build_candidate_substrate_list(
        scores=rank_weighted_fusion_scores,
        top=ADAPTIVE_PARITY_TOP_K,
        score_threshold=ADAPTIVE_PARITY_SCORE_THRESHOLD,
        inclusion=ADAPTIVE_PARITY_INCLUSION,
    )
    observed = run_adaptive_ensemble_prediction(
        prediction_score_matrix=rank_weighted_fusion_scores,
        candidate_substrates=candidate_substrates,
        prediction_config=KinasePredictionConfig(
            top_k=ADAPTIVE_PARITY_TOP_K,
            deterministic_max_selected_kinases=ADAPTIVE_PARITY_ENSEMBLE_SIZE,
            adaptive_ensemble_runs=ADAPTIVE_PARITY_ENSEMBLE_SIZE,
            mode="adaptive_ensemble",
            adaptive_policy=adaptive_policy,
            n_iterations=ADAPTIVE_PARITY_N_ITERATIONS,
            random_state=ADAPTIVE_PARITY_RANDOM_STATE,
        ),
        random_state=ADAPTIVE_PARITY_RANDOM_STATE,
    )
    observed_prediction_frame = _prediction_long_frame(observed)
    expected_prediction_frame = (
        load_adaptive_sampling_edge_trace_predictions()
        .drop_duplicates()
        .reset_index(drop=True)
    )
    donor_merged_prediction = observed_prediction_frame.merge(
        expected_prediction_frame,
        on=["kinase", "site"],
        suffixes=("_py", "_donor"),
        validate="one_to_one",
    )
    prediction_delta = (
        donor_merged_prediction.loc[:, "prob_class_1_py"]
        - donor_merged_prediction.loc[:, "prob_class_1_donor"]
    ).abs()

    observed_top_frame = _top_rank_frame(observed)
    expected_top_frame = load_adaptive_sampling_edge_trace_top().drop_duplicates(
        subset=["kinase", "rank", "site"]
    )
    donor_merged_top = observed_top_frame.merge(
        expected_top_frame,
        on=["kinase", "rank"],
        suffixes=("_py", "_donor"),
        validate="one_to_one",
    )
    donor_top_rank_matches = int(
        (
            donor_merged_top.loc[:, "site_py"] == donor_merged_top.loc[:, "site_donor"]
        ).sum()
    )
    donor_top_set_overlap_matches, donor_top_set_overlap_total = _top_set_overlap(
        stable_top=observed_top_frame,
        r_parity_top=expected_top_frame,
    )

    trace_candidates = load_adaptive_sampling_edge_trace_candidates()
    selected_trace_candidate_count = int(
        trace_candidates.loc[:, "selected_candidate"].sum()
    )
    candidate_count = int(sum(len(sites) for sites in candidate_substrates.values()))
    observed_ranked = _ranked_sites_by_kinase(observed_prediction_frame)
    donor_ranked = _ranked_sites_by_kinase(expected_prediction_frame)

    return AdaptivePolicyLaneMetrics(
        adaptive_policy=adaptive_policy,
        policy_label=POLICY_DISPLAY_LABELS[adaptive_policy],
        prediction_shape=(int(observed.shape[0]), int(observed.shape[1])),
        kinases_compared=int(donor_merged_prediction.loc[:, "kinase"].nunique()),
        candidate_count=candidate_count,
        candidate_kinase_count=int(len(candidate_substrates)),
        selected_trace_candidate_count=selected_trace_candidate_count,
        donor_prediction_rows=int(donor_merged_prediction.shape[0]),
        donor_prediction_corr=float(
            donor_merged_prediction.loc[:, "prob_class_1_py"].corr(
                donor_merged_prediction.loc[:, "prob_class_1_donor"]
            )
        ),
        donor_prediction_mae=float(prediction_delta.mean()),
        donor_prediction_max_abs_diff=float(prediction_delta.max()),
        donor_top_rank_matches=donor_top_rank_matches,
        donor_top_rank_total=int(donor_merged_top.shape[0]),
        donor_top_prob_mae=float(
            (
                donor_merged_top.loc[:, "prob_class_1_py"]
                - donor_merged_top.loc[:, "prob_class_1_donor"]
            )
            .abs()
            .mean()
        ),
        donor_top_set_overlap_matches=donor_top_set_overlap_matches,
        donor_top_set_overlap_total=donor_top_set_overlap_total,
        donor_mean_top10_overlap=_mean_top_n_overlap(
            observed_ranked,
            donor_ranked,
            top_n=10,
        ),
        donor_mean_top20_overlap=_mean_top_n_overlap(
            observed_ranked,
            donor_ranked,
            top_n=20,
        ),
        donor_mean_top30_overlap=_mean_top_n_overlap(
            observed_ranked,
            donor_ranked,
            top_n=30,
        ),
        _observed_prediction_frame=observed_prediction_frame,
        _observed_top_frame=observed_top_frame,
    )


def _prediction_long_frame(observed: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "kinase": str(kinase),
            "site": str(site),
            "prob_class_1": float(score),
        }
        for kinase in observed.columns.astype(str)
        for site, score in observed.loc[:, kinase].items()
    ]
    return pd.DataFrame(rows)


def _top_rank_frame(observed: pd.DataFrame) -> pd.DataFrame:
    top_rows: list[dict[str, object]] = []
    for kinase in observed.columns.astype(str):
        ranked = (
            observed.loc[:, kinase]
            .sort_values(ascending=False)
            .head(ADAPTIVE_PARITY_TOP_K)
        )
        for rank, (site, score) in enumerate(ranked.items(), start=1):
            top_rows.append(
                {
                    "kinase": str(kinase),
                    "rank": int(rank),
                    "site": str(site),
                    "prob_class_1": float(score),
                }
            )
    return pd.DataFrame(top_rows)


def _top_set_overlap(
    *,
    stable_top: pd.DataFrame,
    r_parity_top: pd.DataFrame,
) -> tuple[int, int]:
    stable_by_kinase = {
        str(kinase): set(group.loc[:, "site"].astype(str).tolist())
        for kinase, group in stable_top.groupby("kinase", sort=False)
    }
    r_parity_by_kinase = {
        str(kinase): set(group.loc[:, "site"].astype(str).tolist())
        for kinase, group in r_parity_top.groupby("kinase", sort=False)
    }
    shared_kinases = sorted(set(stable_by_kinase) & set(r_parity_by_kinase))
    overlap_matches = 0
    overlap_total = 0
    for kinase in shared_kinases:
        left = stable_by_kinase[kinase]
        right = r_parity_by_kinase[kinase]
        overlap_matches += len(left & right)
        overlap_total += len(left)
    return overlap_matches, overlap_total


def _ranked_sites_by_kinase(frame: pd.DataFrame) -> dict[str, list[str]]:
    ranked: dict[str, list[str]] = {}
    for kinase, group in frame.groupby("kinase", sort=False):
        ordered = group.sort_values(
            ["prob_class_1", "site"],
            ascending=[False, True],
            kind="mergesort",
        )
        ranked[str(kinase)] = ordered.loc[:, "site"].astype(str).tolist()
    return ranked


def _mean_top_n_overlap(
    ranked_reference: dict[str, list[str]],
    ranked_candidate: dict[str, list[str]],
    *,
    top_n: int,
) -> float:
    shared_kinases = sorted(set(ranked_reference) & set(ranked_candidate))
    overlaps: list[float] = []
    for kinase in shared_kinases:
        reference = ranked_reference[kinase]
        candidate = ranked_candidate[kinase]
        denominator = min(top_n, len(reference), len(candidate))
        if denominator <= 0:
            continue
        reference_top = reference[:denominator]
        candidate_top = candidate[:denominator]
        overlap = len(set(reference_top) & set(candidate_top)) / float(denominator)
        overlaps.append(overlap)
    if not overlaps:
        return 0.0
    return float(sum(overlaps) / len(overlaps))
