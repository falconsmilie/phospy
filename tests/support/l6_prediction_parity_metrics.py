from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import pandas as pd

from phospy import (
    KinaseWorkflow,
)
from phospy.advanced import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api import (
    KinaseWorkflowRequest,
    ReferencePreset,
)
from tests.support.rewrite_fixture_data import (
    build_rat_l6_dataset,
    load_l6_prediction_reference_candidate_substrates,
    load_l6_prediction_reference_predmat,
    load_l6_prediction_reference_profile_scores,
    load_l6_prediction_reference_rank_weighted_fusion_scores,
    load_l6_prediction_reference_score_fusion_weights,
    load_l6_prediction_reference_top30,
)


def _canonical_kinase_label(value: object) -> str:
    return str(value).strip().upper()


def _canonicalise_kinase_columns(frame: pd.DataFrame) -> pd.DataFrame:
    canonical = frame.copy(deep=True)
    canonical.columns = pd.Index(
        [_canonical_kinase_label(value) for value in canonical.columns],
        name=canonical.columns.name,
    )
    return canonical.T.groupby(level=0, sort=False).mean().T


def _canonicalise_kinase_index(frame: pd.DataFrame) -> pd.DataFrame:
    canonical = frame.copy(deep=True)
    canonical.index = pd.Index(
        [_canonical_kinase_label(value) for value in canonical.index],
        name=canonical.index.name,
    )
    return canonical.groupby(level=0, sort=False).mean()


def _canonicalise_kinase_column(
    frame: pd.DataFrame,
    *,
    column_name: str,
) -> pd.DataFrame:
    canonical = frame.copy(deep=True)
    canonical.loc[:, column_name] = canonical.loc[:, column_name].map(
        _canonical_kinase_label
    )
    return canonical


def _with_dataset_display_index(frame: pd.DataFrame, result: object) -> pd.DataFrame:
    dataset = result.dataset
    site_metadata = dataset.site_metadata.reindex(frame.index)
    if "display_id" not in site_metadata.columns:
        return frame
    display_ids = site_metadata.loc[:, "display_id"].astype(str).tolist()
    with_display = frame.copy(deep=True)
    with_display.index = pd.Index(display_ids, name="site_id")
    return with_display


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
class RankingSurfaceContract:
    surface_name: str
    observed_source: str
    expected_source: str
    observed_policy: str
    expected_policy: str


@dataclass(frozen=True, slots=True)
class PredictionMatrixSurface:
    frame: pd.DataFrame


@dataclass(frozen=True, slots=True)
class RankedTopKExportSurface:
    frame: pd.DataFrame


@dataclass(frozen=True, slots=True)
class RankedSiteOrderSurface:
    surface_name: str
    ranked_by_kinase: dict[str, list[str]]


@dataclass(frozen=True, slots=True)
class L6PredictionParityMetrics:
    profile: TableParityMetrics
    combined: TableParityMetrics
    weights: TableParityMetrics
    prediction_matrix: TableParityMetrics
    candidates: CandidateParityMetrics
    prediction_matrix_ranking: RankingParityMetrics
    prediction_matrix_ranking_contract: RankingSurfaceContract
    ranked_topk_export: RankingParityMetrics
    ranked_topk_export_contract: RankingSurfaceContract
    policy_divergence: PolicyDivergenceMetrics


@dataclass(frozen=True, slots=True)
class PolicyDivergenceMetrics:
    prediction_matrix_score_corr: float
    prediction_matrix_score_mae: float
    prediction_matrix_ranking: RankingParityMetrics
    prediction_matrix_ranking_contract: RankingSurfaceContract
    ranked_topk_export: RankingParityMetrics
    ranked_topk_export_contract: RankingSurfaceContract


def _stack_frame(frame: pd.DataFrame) -> pd.Series:
    try:
        return frame.stack(future_stack=True)
    except TypeError:
        # pandas<2.1 compatibility path (no future_stack argument)
        return frame.stack(dropna=False)


def _safe_corr(
    left: pd.Series,
    right: pd.Series,
    *,
    method: str,
) -> float | None:
    aligned = (
        pd.concat(
            [left.astype(float), right.astype(float)],
            axis=1,
            join="inner",
        )
        .dropna()
        .reset_index(drop=True)
    )
    if aligned.shape[0] < 2:
        return None
    left_aligned = aligned.iloc[:, 0]
    right_aligned = aligned.iloc[:, 1]
    if left_aligned.nunique(dropna=True) <= 1:
        return None
    if right_aligned.nunique(dropna=True) <= 1:
        return None
    corr = left_aligned.corr(right_aligned, method=method)
    if pd.isna(corr):
        return None
    return float(corr)


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
        corr = _safe_corr(
            observed.loc[common_index, column],
            expected.loc[common_index, column],
            method=method,
        )
        if corr is not None:
            correlations.append(corr)
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


def _spearman_rank_corr_from_ranked_sites(
    expected_ranked: list[str],
    observed_ranked: list[str],
) -> float | None:
    expected_rank_by_site = {
        site_id: rank for rank, site_id in enumerate(expected_ranked, start=1)
    }
    observed_rank_by_site = {
        site_id: rank for rank, site_id in enumerate(observed_ranked, start=1)
    }
    shared_sites = [
        site_id for site_id in expected_ranked if site_id in observed_rank_by_site
    ]
    if len(shared_sites) < 2:
        return None
    expected_ranks = pd.Series(
        [expected_rank_by_site[site_id] for site_id in shared_sites],
        dtype=float,
    )
    observed_ranks = pd.Series(
        [observed_rank_by_site[site_id] for site_id in shared_sites],
        dtype=float,
    )
    return _safe_corr(expected_ranks, observed_ranks, method="spearman")


def _collect_ranking_metrics(
    *,
    observed: RankedSiteOrderSurface,
    expected: RankedSiteOrderSurface,
) -> RankingParityMetrics:
    if observed.surface_name != expected.surface_name:
        raise AssertionError(
            "ranking parity expects aligned surfaces "
            f"(observed={observed.surface_name}, expected={expected.surface_name})"
        )

    kinases = sorted(set(observed.ranked_by_kinase) & set(expected.ranked_by_kinase))
    top10: list[float] = []
    top20: list[float] = []
    top30: list[float] = []
    spearman: list[float] = []
    top_rank_matches = 0

    for kinase in kinases:
        observed_ranked = observed.ranked_by_kinase[kinase]
        expected_ranked = expected.ranked_by_kinase[kinase]
        if not expected_ranked:
            continue
        top10.append(_top_n_overlap(expected_ranked, observed_ranked, 10))
        top20.append(_top_n_overlap(expected_ranked, observed_ranked, 20))
        top30.append(_top_n_overlap(expected_ranked, observed_ranked, 30))
        if (
            observed_ranked
            and expected_ranked
            and observed_ranked[0] == expected_ranked[0]
        ):
            top_rank_matches += 1
        corr = _spearman_rank_corr_from_ranked_sites(expected_ranked, observed_ranked)
        if corr is not None:
            spearman.append(corr)

    kinases_compared = len(top10)
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


def _rank_prediction_matrix_surface(
    surface: PredictionMatrixSurface,
) -> RankedSiteOrderSurface:
    canonical_frame = _canonicalise_kinase_columns(surface.frame)
    ranked_by_kinase: dict[str, list[str]] = {}
    for kinase in canonical_frame.columns.astype(str):
        ordered = (
            canonical_frame.loc[:, kinase]
            .astype(float)
            .sort_values(ascending=False, kind="mergesort")
            .index.astype(str)
            .tolist()
        )
        ranked_by_kinase[kinase] = ordered
    return RankedSiteOrderSurface(
        surface_name="prediction_matrix",
        ranked_by_kinase=ranked_by_kinase,
    )


def _normalize_prediction_topk_frame(frame: pd.DataFrame) -> pd.DataFrame:
    normalized = frame.rename(columns={"substrate_site": "site_id"}).copy(deep=True)
    required_columns = {"kinase", "site_id"}
    if not required_columns.issubset(normalized.columns):
        raise AssertionError(
            "expected top-k prediction frame to contain kinase/site_id columns"
        )
    normalized = normalized.astype({"kinase": str, "site_id": str})
    normalized.loc[:, "kinase"] = normalized.loc[:, "kinase"].map(
        _canonical_kinase_label
    )
    if "rank" in normalized.columns:
        normalized.loc[:, "rank"] = normalized.loc[:, "rank"].astype("int64")
    else:
        if "pred_score" in normalized.columns:
            normalized = normalized.sort_values(
                ["kinase", "pred_score"],
                ascending=[True, False],
                kind="mergesort",
            )
        else:
            normalized = normalized.sort_values(["kinase"], kind="mergesort")
        normalized.loc[:, "rank"] = (
            normalized.groupby("kinase", sort=False).cumcount() + 1
        )
    return (
        normalized.loc[:, ["kinase", "site_id", "rank"]]
        .sort_values(["kinase", "rank"], kind="mergesort")
        .reset_index(drop=True)
    )


def _rank_topk_export_surface(
    surface: RankedTopKExportSurface,
) -> RankedSiteOrderSurface:
    ranked_by_kinase: dict[str, list[str]] = {}
    for kinase, group in surface.frame.groupby("kinase", sort=False):
        ordered = group.sort_values("rank", kind="mergesort").loc[:, "site_id"].tolist()
        ranked_by_kinase[str(kinase)] = [str(site_id) for site_id in ordered]
    return RankedSiteOrderSurface(
        surface_name="ranked_topk_export",
        ranked_by_kinase=ranked_by_kinase,
    )


def _collect_policy_divergence_metrics(
    *,
    stable_prediction_matrix: PredictionMatrixSurface,
    r_parity_prediction_matrix: PredictionMatrixSurface,
    stable_topk_export: RankedTopKExportSurface,
    r_parity_topk_export: RankedTopKExportSurface,
) -> PolicyDivergenceMetrics:
    stable_prediction_frame = _canonicalise_kinase_columns(
        stable_prediction_matrix.frame
    )
    r_parity_prediction_frame = _canonicalise_kinase_columns(
        r_parity_prediction_matrix.frame
    )
    stable_long = (
        _stack_frame(stable_prediction_frame)
        .rename("score_stable")
        .reset_index()
        .rename(
            columns={
                "level_1": "kinase",
                stable_prediction_frame.index.name or "level_0": "site_id",
            }
        )
    )
    r_parity_long = (
        _stack_frame(r_parity_prediction_frame)
        .rename("score_r_parity")
        .reset_index()
        .rename(
            columns={
                "level_1": "kinase",
                r_parity_prediction_frame.index.name or "level_0": "site_id",
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
    absolute_delta = (
        merged_non_null.loc[:, "score_stable"]
        - merged_non_null.loc[:, "score_r_parity"]
    ).abs()
    stable_pred_ranked = _rank_prediction_matrix_surface(stable_prediction_matrix)
    r_parity_pred_ranked = _rank_prediction_matrix_surface(r_parity_prediction_matrix)
    prediction_matrix_contract = RankingSurfaceContract(
        surface_name="prediction_matrix",
        observed_source="rewrite_prediction_matrix",
        expected_source="rewrite_prediction_matrix",
        observed_policy="r_parity",
        expected_policy="stable",
    )
    ranked_topk_contract = RankingSurfaceContract(
        surface_name="ranked_topk_export",
        observed_source="rewrite_ranked_topk_export",
        expected_source="rewrite_ranked_topk_export",
        observed_policy="r_parity",
        expected_policy="stable",
    )
    return PolicyDivergenceMetrics(
        prediction_matrix_score_corr=(
            _safe_corr(
                merged_non_null.loc[:, "score_stable"],
                merged_non_null.loc[:, "score_r_parity"],
                method="pearson",
            )
            or 0.0
        )
        if not merged_non_null.empty
        else 0.0,
        prediction_matrix_score_mae=float(absolute_delta.mean())
        if not absolute_delta.empty
        else 0.0,
        prediction_matrix_ranking=_collect_ranking_metrics(
            observed=r_parity_pred_ranked,
            expected=stable_pred_ranked,
        ),
        prediction_matrix_ranking_contract=prediction_matrix_contract,
        ranked_topk_export=_collect_ranking_metrics(
            observed=_rank_topk_export_surface(r_parity_topk_export),
            expected=_rank_topk_export_surface(stable_topk_export),
        ),
        ranked_topk_export_contract=ranked_topk_contract,
    )


@lru_cache(maxsize=1)
def _run_l6_workflow(adaptive_policy: str):
    dataset = build_rat_l6_dataset(n_sites=None)
    return KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=ReferencePreset.AUTO,
            scoring_config=KinaseScoringConfig(
                reliability_profile="custom",
                min_substrates=2,
                include_diagnostic_scoring_tables=True,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
            prediction_config=KinasePredictionConfig(
                top_k=30,
                deterministic_max_selected_kinases=10,
                adaptive_ensemble_runs=10,
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

    observed_profile = _canonicalise_kinase_columns(
        _with_dataset_display_index(
            stable_result.scoring_result.profile_scores,
            stable_result,
        )
    )
    expected_profile = _canonicalise_kinase_columns(
        load_l6_prediction_reference_profile_scores()
    )
    observed_combined = stable_result.scoring_result.rank_weighted_fusion_scores
    if observed_combined is None:
        raise AssertionError(
            "expected rank_weighted_fusion_scores to be present for L6 parity"
        )
    observed_combined = _canonicalise_kinase_columns(
        _with_dataset_display_index(observed_combined, stable_result)
    )
    expected_combined = _canonicalise_kinase_columns(
        load_l6_prediction_reference_rank_weighted_fusion_scores()
    )
    observed_weights = stable_result.scoring_result.score_fusion_weights
    if observed_weights is None:
        raise AssertionError("expected score_fusion_weights table for L6 parity")
    observed_weights = _canonicalise_kinase_index(observed_weights)
    expected_weights = _canonicalise_kinase_index(
        load_l6_prediction_reference_score_fusion_weights()
    )
    expected_candidates = _canonicalise_kinase_column(
        load_l6_prediction_reference_candidate_substrates(),
        column_name="kinase",
    )
    expected_pred_mat = _canonicalise_kinase_columns(
        load_l6_prediction_reference_predmat()
    )
    expected_topk_export_surface = RankedTopKExportSurface(
        frame=_normalize_prediction_topk_frame(load_l6_prediction_reference_top30())
    )

    stable_substrate_list = stable_result.prediction_result.substrate_list
    if stable_substrate_list is None:
        raise AssertionError("expected substrate_list to be present for L6 parity")
    stable_topk_export_surface = RankedTopKExportSurface(
        frame=_normalize_prediction_topk_frame(stable_substrate_list)
    )
    stable_candidates_frame = stable_topk_export_surface.frame.loc[
        :, ["kinase", "site_id"]
    ]
    r_parity_substrate_list = r_parity_result.prediction_result.substrate_list
    if r_parity_substrate_list is None:
        raise AssertionError("expected substrate_list to be present for L6 r_parity")
    r_parity_topk_export_surface = RankedTopKExportSurface(
        frame=_normalize_prediction_topk_frame(r_parity_substrate_list)
    )

    stable_prediction_matrix_surface = PredictionMatrixSurface(
        frame=_canonicalise_kinase_columns(
            _with_dataset_display_index(
                stable_result.prediction_result.pred_mat,
                stable_result,
            )
        )
    )
    r_parity_prediction_matrix_surface = PredictionMatrixSurface(
        frame=_canonicalise_kinase_columns(
            _with_dataset_display_index(
                r_parity_result.prediction_result.pred_mat,
                r_parity_result,
            )
        )
    )
    expected_prediction_matrix_surface = PredictionMatrixSurface(
        frame=expected_pred_mat
    )
    expected_pred_mat_ranked = _rank_prediction_matrix_surface(
        expected_prediction_matrix_surface
    )
    stable_pred_ranked = _rank_prediction_matrix_surface(
        stable_prediction_matrix_surface
    )
    expected_topk_ranked = _rank_topk_export_surface(expected_topk_export_surface)
    stable_topk_ranked = _rank_topk_export_surface(stable_topk_export_surface)
    prediction_matrix_contract = RankingSurfaceContract(
        surface_name="prediction_matrix",
        observed_source="rewrite_prediction_matrix",
        expected_source="promoted_reference_prediction_matrix",
        observed_policy="stable",
        expected_policy="stable",
    )
    ranked_topk_contract = RankingSurfaceContract(
        surface_name="ranked_topk_export",
        observed_source="rewrite_ranked_topk_export",
        expected_source="promoted_reference_ranked_topk_export",
        observed_policy="stable",
        expected_policy="stable",
    )

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
        prediction_matrix=_collect_table_parity_metrics(
            observed=stable_prediction_matrix_surface.frame,
            expected=expected_pred_mat,
        ),
        candidates=_collect_candidate_metrics(
            observed=stable_candidates_frame,
            expected=expected_candidates,
        ),
        prediction_matrix_ranking=_collect_ranking_metrics(
            observed=stable_pred_ranked,
            expected=expected_pred_mat_ranked,
        ),
        prediction_matrix_ranking_contract=prediction_matrix_contract,
        ranked_topk_export=_collect_ranking_metrics(
            observed=stable_topk_ranked,
            expected=expected_topk_ranked,
        ),
        ranked_topk_export_contract=ranked_topk_contract,
        policy_divergence=_collect_policy_divergence_metrics(
            stable_prediction_matrix=stable_prediction_matrix_surface,
            r_parity_prediction_matrix=r_parity_prediction_matrix_surface,
            stable_topk_export=stable_topk_export_surface,
            r_parity_topk_export=r_parity_topk_export_surface,
        ),
    )
