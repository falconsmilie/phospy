from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from phospy.prediction.candidates import build_candidate_substrate_list
from phospy.prediction.policies import (
    PredictionSamplingRandomSource,
    resolve_prediction_sampling_policy,
)
from phospy.prediction.sampling_runtime import (
    normalize_probabilities,
    transform_resampling_probabilities,
)
from phospy.prediction.svm import (
    aligned_binary_decision_vector,
    make_svm,
    require_sklearn,
)
from tests.support.rewrite_fixture_data import (
    load_adaptive_sampling_replay_rank_weighted_fusion_scores,
    load_adaptive_sampling_replay_trace_candidates,
    load_adaptive_sampling_replay_trace_final_predictions,
    load_adaptive_sampling_replay_trace_final_top,
    load_adaptive_sampling_replay_trace_initial_negatives,
    load_adaptive_sampling_replay_trace_iteration_samples,
)

TRACE_TOP_K = 30
TRACE_SCORE_THRESHOLD = 0.6
TRACE_INCLUSION = 5
TRACE_ENSEMBLE_SIZE = 10
TRACE_N_ITERATIONS = 5
TRACE_RANDOM_STATE = 1


@dataclass(frozen=True, slots=True)
class AdaptiveTraceLaneMetrics:
    adaptive_policy: str
    candidate_count: int
    candidate_kinase_count: int
    selected_trace_kinase_count: int
    initial_overlap_matches: int
    initial_overlap_total: int
    sample_overlap_matches: int
    sample_overlap_total: int
    donor_final_prediction_corr: float
    donor_final_prediction_mae: float
    donor_final_prediction_max_abs_diff: float
    donor_top_rank_matches: int
    donor_top_rank_total: int
    donor_top_set_overlap_matches: int
    donor_top_set_overlap_total: int
    donor_mean_top10_overlap: float
    donor_mean_top20_overlap: float
    donor_mean_top30_overlap: float
    deterministic_under_seed: bool
    _aggregated_predictions: pd.DataFrame
    _top_frame: pd.DataFrame


@dataclass(frozen=True, slots=True)
class AdaptiveTraceComparisonMetrics:
    stable: AdaptiveTraceLaneMetrics
    r_parity: AdaptiveTraceLaneMetrics
    cross_policy_prediction_corr: float
    cross_policy_prediction_mae: float
    cross_policy_top_set_overlap_matches: int
    cross_policy_top_set_overlap_total: int
    cross_policy_mean_top10_overlap: float
    cross_policy_mean_top20_overlap: float
    cross_policy_mean_top30_overlap: float


@dataclass(frozen=True, slots=True)
class _TracedLaneResult:
    initial_negatives: pd.DataFrame
    iteration_samples: pd.DataFrame
    final_predictions: pd.DataFrame
    final_top: pd.DataFrame
    aggregated_predictions: pd.DataFrame


def _eligible_trace_kinases(
    *,
    candidate_substrates: dict[str, list[str]],
) -> list[str]:
    expected_candidates = load_adaptive_sampling_replay_trace_candidates()
    expected_initial = load_adaptive_sampling_replay_trace_initial_negatives()
    expected_samples = load_adaptive_sampling_replay_trace_iteration_samples()
    trace_selected_counts = (
        expected_candidates.loc[expected_candidates.loc[:, "selected_candidate"]]
        .groupby("kinase")
        .size()
        .to_dict()
    )

    eligible_kinases: list[str] = []
    for kinase in sorted(trace_selected_counts):
        trace_count = int(trace_selected_counts[kinase])
        current_count = len(candidate_substrates.get(str(kinase), []))
        initial_counts = (
            expected_initial.loc[
                expected_initial.loc[:, "kinase"].astype(str) == kinase
            ]
            .groupby("ensemble")
            .size()
        )
        sample_counts = (
            expected_samples.loc[
                expected_samples.loc[:, "kinase"].astype(str) == kinase
            ]
            .groupby(["ensemble", "iteration", "class_label"])
            .size()
        )
        initial_ok = (not initial_counts.empty) and bool(
            (initial_counts == trace_count).all()
        )
        sample_ok = (not sample_counts.empty) and bool(
            (sample_counts == trace_count).all()
        )
        if current_count == trace_count and initial_ok and sample_ok:
            eligible_kinases.append(str(kinase))
    return eligible_kinases


def _run_traced_lane(
    *,
    adaptive_policy: str,
    kinases: list[str],
    initial_overrides: pd.DataFrame | None = None,
    sample_overrides: pd.DataFrame | None = None,
) -> _TracedLaneResult:
    rank_weighted_fusion_scores = (
        load_adaptive_sampling_replay_rank_weighted_fusion_scores()
    )
    candidate_substrates = build_candidate_substrate_list(
        scores=rank_weighted_fusion_scores,
        top=TRACE_TOP_K,
        score_threshold=TRACE_SCORE_THRESHOLD,
        inclusion=TRACE_INCLUSION,
    )
    policy = resolve_prediction_sampling_policy(adaptive_policy)
    random_source = PredictionSamplingRandomSource(
        policy=policy,
        random_state=TRACE_RANDOM_STATE,
    )
    standard_scaler, svc = require_sklearn()

    feature_values = rank_weighted_fusion_scores.to_numpy(dtype=float, copy=False)
    feature_index = rank_weighted_fusion_scores.index
    all_positions = np.arange(feature_values.shape[0], dtype=int)
    site_position = {
        str(site_id): int(position) for position, site_id in enumerate(feature_index)
    }

    initial_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []
    final_prediction_rows: list[dict[str, object]] = []
    final_top_rows: list[dict[str, object]] = []

    for kinase in kinases:
        positive_sites = candidate_substrates.get(kinase, [])
        positive_positions = np.asarray(
            [
                site_position[site_id]
                for site_id in positive_sites
                if site_id in site_position
            ],
            dtype=int,
        )
        if positive_positions.size == 0:
            continue
        negative_mask = np.ones(len(all_positions), dtype=bool)
        negative_mask[positive_positions] = False
        negative_positions = all_positions[negative_mask]
        if negative_positions.size == 0:
            continue
        positive_values = feature_values[positive_positions, :]
        base_labels = np.concatenate(
            [
                np.repeat(1, len(positive_positions)),
                np.repeat(2, len(positive_positions)),
            ]
        )
        negative_rng, resampling_rng = random_source.generators_for_kinase(
            kinase=kinase
        )

        for ensemble in range(1, TRACE_ENSEMBLE_SIZE + 1):
            sampled_negative_positions: np.ndarray
            if initial_overrides is not None:
                initial_site_overrides = (
                    initial_overrides.loc[
                        (initial_overrides.loc[:, "kinase"].astype(str) == kinase)
                        & (initial_overrides.loc[:, "ensemble"].astype(int) == ensemble)
                    ]
                    .sort_values("draw")
                    .loc[:, "site"]
                    .astype(str)
                    .tolist()
                )
                if len(initial_site_overrides) == len(positive_positions):
                    sampled_negative_positions = np.asarray(
                        [
                            site_position[site_id]
                            for site_id in initial_site_overrides
                            if site_id in site_position
                        ],
                        dtype=int,
                    )
                    if len(sampled_negative_positions) != len(positive_positions):
                        sampled_negative_positions = negative_rng.choice(
                            negative_positions,
                            size=len(positive_positions),
                            replace=len(negative_positions) < len(positive_positions),
                        )
                else:
                    sampled_negative_positions = negative_rng.choice(
                        negative_positions,
                        size=len(positive_positions),
                        replace=len(negative_positions) < len(positive_positions),
                    )
            else:
                sampled_negative_positions = negative_rng.choice(
                    negative_positions,
                    size=len(positive_positions),
                    replace=len(negative_positions) < len(positive_positions),
                )
            for draw, position in enumerate(sampled_negative_positions, start=1):
                initial_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": int(ensemble),
                        "draw": int(draw),
                        "site": str(feature_index[int(position)]),
                    }
                )

            base_positions = np.concatenate(
                [positive_positions, sampled_negative_positions]
            )
            base_values = np.concatenate(
                [
                    positive_values,
                    feature_values[sampled_negative_positions, :],
                ],
                axis=0,
            )
            current_values = base_values
            current_labels = base_labels
            model = None

            for iteration in range(1, TRACE_N_ITERATIONS + 1):
                model = make_svm(
                    StandardScaler=standard_scaler,
                    SVC=svc,
                    kernel="rbf",
                    use_r_parity_scaler=policy.adaptive_policy == "r_parity",
                )
                model.fit(current_values, current_labels)
                probability_matrix = model.predict_proba(base_values)
                model_classes = np.asarray(model.classes_, dtype=int)
                positive_idx = np.flatnonzero(model_classes == 1)
                positive_probabilities = (
                    probability_matrix[:, int(positive_idx[0])]
                    if len(positive_idx) == 1
                    else None
                )
                _ = aligned_binary_decision_vector(
                    model=model,
                    values=base_values,
                    positive_probabilities=positive_probabilities,
                )
                class_weights: dict[int, np.ndarray | None] = {}
                for class_index, class_label in enumerate(model_classes):
                    class_mask = base_labels == class_label
                    transformed = transform_resampling_probabilities(
                        probability_matrix[class_mask, class_index],
                        sampling_policy=policy,
                    )
                    class_weights[int(class_label)] = normalize_probabilities(
                        transformed
                    )

                resampled_values: list[np.ndarray] = []
                resampled_labels: list[np.ndarray] = []
                for class_label in model_classes:
                    class_mask = base_labels == class_label
                    class_values = base_values[class_mask]
                    class_positions = base_positions[class_mask]
                    sampled_indices: np.ndarray
                    if sample_overrides is not None:
                        sample_rows_expected = (
                            sample_overrides.loc[
                                (
                                    sample_overrides.loc[:, "kinase"].astype(str)
                                    == kinase
                                )
                                & (
                                    sample_overrides.loc[:, "ensemble"].astype(int)
                                    == ensemble
                                )
                                & (
                                    sample_overrides.loc[:, "iteration"].astype(int)
                                    == iteration
                                )
                                & (
                                    sample_overrides.loc[:, "class_label"].astype(str)
                                    == str(int(class_label))
                                )
                            ]
                            .sort_values("draw")
                            .loc[:, "site"]
                            .astype(str)
                            .tolist()
                        )
                        if len(sample_rows_expected) == class_values.shape[0]:
                            class_site_lookup: dict[str, int] = {
                                str(feature_index[int(position)]): idx
                                for idx, position in enumerate(class_positions)
                            }
                            sampled_indices = np.asarray(
                                [
                                    class_site_lookup[site_id]
                                    for site_id in sample_rows_expected
                                    if site_id in class_site_lookup
                                ],
                                dtype=int,
                            )
                            if len(sampled_indices) != class_values.shape[0]:
                                sampled_indices = resampling_rng.choice(
                                    class_values.shape[0],
                                    size=class_values.shape[0],
                                    replace=True,
                                    p=class_weights[int(class_label)],
                                )
                        else:
                            sampled_indices = resampling_rng.choice(
                                class_values.shape[0],
                                size=class_values.shape[0],
                                replace=True,
                                p=class_weights[int(class_label)],
                            )
                    else:
                        sampled_indices = resampling_rng.choice(
                            class_values.shape[0],
                            size=class_values.shape[0],
                            replace=True,
                            p=class_weights[int(class_label)],
                        )
                    for draw, sampled_idx in enumerate(sampled_indices, start=1):
                        sample_rows.append(
                            {
                                "kinase": kinase,
                                "ensemble": int(ensemble),
                                "iteration": int(iteration),
                                "class_label": str(int(class_label)),
                                "draw": int(draw),
                                "site": str(
                                    feature_index[
                                        int(class_positions[int(sampled_idx)])
                                    ]
                                ),
                            }
                        )
                    resampled_values.append(class_values[sampled_indices])
                    resampled_labels.append(
                        np.repeat(class_label, class_values.shape[0])
                    )
                current_values = np.vstack(resampled_values)
                current_labels = np.concatenate(resampled_labels)

            if model is None:
                continue
            final_probabilities = np.asarray(
                model.predict_proba(feature_values), dtype=float
            )
            model_classes = np.asarray(model.classes_, dtype=int)
            class_to_column = {
                int(class_label): int(column)
                for column, class_label in enumerate(model_classes)
            }
            prob_class_1 = final_probabilities[:, class_to_column[1]]
            prob_class_2 = final_probabilities[:, class_to_column[2]]
            for position, site_id in enumerate(feature_index):
                final_prediction_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": int(ensemble),
                        "site": str(site_id),
                        "prob_class_1": float(prob_class_1[position]),
                        "prob_class_2": float(prob_class_2[position]),
                    }
                )
            ranked_positions = np.argsort(-prob_class_1, kind="stable")[:TRACE_TOP_K]
            for rank, position in enumerate(ranked_positions, start=1):
                final_top_rows.append(
                    {
                        "kinase": kinase,
                        "ensemble": int(ensemble),
                        "rank": int(rank),
                        "site": str(feature_index[int(position)]),
                        "prob_class_1": float(prob_class_1[int(position)]),
                    }
                )

    initial_negatives = pd.DataFrame(initial_rows)
    iteration_samples = pd.DataFrame(sample_rows)
    final_predictions = pd.DataFrame(final_prediction_rows)
    final_top = pd.DataFrame(final_top_rows)
    aggregated_predictions = final_predictions.groupby(
        ["kinase", "site"], as_index=False
    )["prob_class_1"].mean()
    return _TracedLaneResult(
        initial_negatives=initial_negatives,
        iteration_samples=iteration_samples,
        final_predictions=final_predictions,
        final_top=final_top,
        aggregated_predictions=aggregated_predictions,
    )


def _set_overlap(
    *,
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    columns: list[str],
) -> tuple[int, int]:
    observed_set = set(
        map(tuple, observed.loc[:, columns].itertuples(index=False, name=None))
    )
    expected_set = set(
        map(tuple, expected.loc[:, columns].itertuples(index=False, name=None))
    )
    return int(len(observed_set & expected_set)), int(len(expected_set))


def _top_set_overlap(
    *,
    left: pd.DataFrame,
    right: pd.DataFrame,
) -> tuple[int, int]:
    left_by_lane = {
        (str(kinase), int(ensemble)): set(group.loc[:, "site"].astype(str).tolist())
        for (kinase, ensemble), group in left.groupby(
            ["kinase", "ensemble"], sort=False
        )
    }
    right_by_lane = {
        (str(kinase), int(ensemble)): set(group.loc[:, "site"].astype(str).tolist())
        for (kinase, ensemble), group in right.groupby(
            ["kinase", "ensemble"], sort=False
        )
    }
    shared = sorted(set(left_by_lane) & set(right_by_lane))
    overlap_matches = 0
    overlap_total = 0
    for key in shared:
        overlap_matches += len(left_by_lane[key] & right_by_lane[key])
        overlap_total += len(right_by_lane[key])
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
        overlap = len(
            set(reference[:denominator]) & set(candidate[:denominator])
        ) / float(denominator)
        overlaps.append(overlap)
    if not overlaps:
        return 0.0
    return float(pd.Series(overlaps).mean())


def _collect_lane_metrics(*, adaptive_policy: str) -> AdaptiveTraceLaneMetrics:
    rank_weighted_fusion_scores = (
        load_adaptive_sampling_replay_rank_weighted_fusion_scores()
    )
    candidate_substrates = build_candidate_substrate_list(
        scores=rank_weighted_fusion_scores,
        top=TRACE_TOP_K,
        score_threshold=TRACE_SCORE_THRESHOLD,
        inclusion=TRACE_INCLUSION,
    )
    eligible_kinases = _eligible_trace_kinases(
        candidate_substrates=candidate_substrates
    )
    if not eligible_kinases:
        raise AssertionError(
            "no replay-aligned kinases available for adaptive trace parity"
        )

    expected_initial = load_adaptive_sampling_replay_trace_initial_negatives()
    expected_samples = load_adaptive_sampling_replay_trace_iteration_samples().astype(
        {"class_label": str}
    )
    expected_final = load_adaptive_sampling_replay_trace_final_predictions()
    expected_top = load_adaptive_sampling_replay_trace_final_top()
    expected_initial = expected_initial.loc[
        expected_initial.loc[:, "kinase"].isin(eligible_kinases)
    ]
    expected_samples = expected_samples.loc[
        expected_samples.loc[:, "kinase"].isin(eligible_kinases)
    ]
    expected_final = expected_final.loc[
        expected_final.loc[:, "kinase"].isin(eligible_kinases)
    ]
    expected_top = expected_top.loc[
        expected_top.loc[:, "kinase"].isin(eligible_kinases)
    ]

    observed = _run_traced_lane(
        adaptive_policy=adaptive_policy,
        kinases=eligible_kinases,
        initial_overrides=expected_initial,
        sample_overrides=expected_samples,
    )
    repeated = _run_traced_lane(
        adaptive_policy=adaptive_policy,
        kinases=eligible_kinases,
        initial_overrides=expected_initial,
        sample_overrides=expected_samples,
    )
    deterministic = (
        observed.initial_negatives.equals(repeated.initial_negatives)
        and observed.iteration_samples.equals(repeated.iteration_samples)
        and observed.final_predictions.equals(repeated.final_predictions)
        and observed.final_top.equals(repeated.final_top)
    )

    initial_overlap_matches, initial_overlap_total = _set_overlap(
        observed=observed.initial_negatives,
        expected=expected_initial,
        columns=["kinase", "ensemble", "site"],
    )
    sample_overlap_matches, sample_overlap_total = _set_overlap(
        observed=observed.iteration_samples,
        expected=expected_samples,
        columns=["kinase", "ensemble", "iteration", "class_label", "site"],
    )

    merged_predictions = observed.final_predictions.merge(
        expected_final,
        on=["kinase", "ensemble", "site"],
        suffixes=("_py", "_donor"),
        validate="one_to_one",
    )
    prediction_delta = (
        merged_predictions.loc[:, "prob_class_1_py"]
        - merged_predictions.loc[:, "prob_class_1_donor"]
    ).abs()

    merged_top = observed.final_top.merge(
        expected_top,
        on=["kinase", "ensemble", "rank"],
        suffixes=("_py", "_donor"),
        validate="one_to_one",
    )
    top_rank_matches = int(
        (merged_top.loc[:, "site_py"] == merged_top.loc[:, "site_donor"]).sum()
    )
    top_set_overlap_matches, top_set_overlap_total = _top_set_overlap(
        left=observed.final_top,
        right=expected_top,
    )

    expected_aggregated = expected_final.groupby(["kinase", "site"], as_index=False)[
        "prob_class_1"
    ].mean()
    observed_ranked = _ranked_sites_by_kinase(observed.aggregated_predictions)
    expected_ranked = _ranked_sites_by_kinase(expected_aggregated)

    return AdaptiveTraceLaneMetrics(
        adaptive_policy=adaptive_policy,
        candidate_count=int(sum(len(sites) for sites in candidate_substrates.values())),
        candidate_kinase_count=int(len(candidate_substrates)),
        selected_trace_kinase_count=int(len(eligible_kinases)),
        initial_overlap_matches=initial_overlap_matches,
        initial_overlap_total=initial_overlap_total,
        sample_overlap_matches=sample_overlap_matches,
        sample_overlap_total=sample_overlap_total,
        donor_final_prediction_corr=float(
            merged_predictions.loc[:, "prob_class_1_py"].corr(
                merged_predictions.loc[:, "prob_class_1_donor"]
            )
        ),
        donor_final_prediction_mae=float(prediction_delta.mean()),
        donor_final_prediction_max_abs_diff=float(prediction_delta.max()),
        donor_top_rank_matches=top_rank_matches,
        donor_top_rank_total=int(merged_top.shape[0]),
        donor_top_set_overlap_matches=top_set_overlap_matches,
        donor_top_set_overlap_total=top_set_overlap_total,
        donor_mean_top10_overlap=_mean_top_n_overlap(
            observed_ranked,
            expected_ranked,
            top_n=10,
        ),
        donor_mean_top20_overlap=_mean_top_n_overlap(
            observed_ranked,
            expected_ranked,
            top_n=20,
        ),
        donor_mean_top30_overlap=_mean_top_n_overlap(
            observed_ranked,
            expected_ranked,
            top_n=30,
        ),
        deterministic_under_seed=deterministic,
        _aggregated_predictions=observed.aggregated_predictions,
        _top_frame=observed.final_top,
    )


@lru_cache(maxsize=1)
def collect_adaptive_trace_replay_metrics() -> AdaptiveTraceComparisonMetrics:
    stable = _collect_lane_metrics(adaptive_policy="stable")
    r_parity = _collect_lane_metrics(adaptive_policy="r_parity")

    merged = stable._aggregated_predictions.merge(
        r_parity._aggregated_predictions,
        on=["kinase", "site"],
        suffixes=("_stable", "_r_parity"),
        validate="one_to_one",
    )
    cross_delta = (
        merged.loc[:, "prob_class_1_stable"] - merged.loc[:, "prob_class_1_r_parity"]
    ).abs()
    cross_top_set_overlap_matches, cross_top_set_overlap_total = _top_set_overlap(
        left=stable._top_frame,
        right=r_parity._top_frame,
    )
    stable_ranked = _ranked_sites_by_kinase(stable._aggregated_predictions)
    r_parity_ranked = _ranked_sites_by_kinase(r_parity._aggregated_predictions)

    return AdaptiveTraceComparisonMetrics(
        stable=stable,
        r_parity=r_parity,
        cross_policy_prediction_corr=float(
            merged.loc[:, "prob_class_1_stable"].corr(
                merged.loc[:, "prob_class_1_r_parity"]
            )
        ),
        cross_policy_prediction_mae=float(cross_delta.mean()),
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
