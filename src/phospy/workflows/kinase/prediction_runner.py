"""Prediction orchestration for kinase workflow execution."""

from __future__ import annotations

from collections.abc import Callable
from typing import NoReturn

import pandas as pd

from phospy.contracts.configs import (
    KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KinasePredictionConfig,
)
from phospy.errors.workflows import WorkflowBoundaryError, WorkflowStageError
from phospy.science.prediction.candidates import (
    CandidateShortfallDiagnostics,
    build_candidate_substrate_list,
    summarize_candidate_shortfall,
)
from phospy.science.prediction.execution import run_adaptive_ensemble_prediction
from phospy.science.prediction.models import KinasePredictionResult
from phospy.workflows.kinase.component_models import (
    CANDIDATE_MIN_INCLUSION,
    CANDIDATE_SCORE_THRESHOLD,
    KinaseScoringRunResult,
)
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.science import (
    build_prediction_outputs,
    rank_kinases_for_prediction,
)


class KinasePredictionRunner:
    """Run deterministic/adaptive prediction lanes from scoring outputs."""

    def __init__(
        self,
        *,
        build_candidates: Callable[..., dict[str, list[str]]] = (
            build_candidate_substrate_list
        ),
        summarize_candidate_shortfall_fn: Callable[
            ..., CandidateShortfallDiagnostics
        ] = (summarize_candidate_shortfall),
        run_adaptive_prediction: Callable[..., pd.DataFrame] = (
            run_adaptive_ensemble_prediction
        ),
        rank_kinases: Callable[..., pd.Series] = rank_kinases_for_prediction,
        build_outputs: Callable[..., tuple[pd.DataFrame, pd.DataFrame]] = (
            build_prediction_outputs
        ),
    ) -> None:
        self._build_candidates = build_candidates
        self._summarize_candidate_shortfall = summarize_candidate_shortfall_fn
        self._run_adaptive_prediction = run_adaptive_prediction
        self._rank_kinases = rank_kinases
        self._build_outputs = build_outputs

    def run(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        scoring_execution: KinaseScoringRunResult,
    ) -> KinasePredictionResult:
        if config.prediction_mode == KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING:
            return self._run_deterministic_prediction_lane(
                request=request,
                config=config,
                scoring_execution=scoring_execution,
            )
        if config.prediction_mode == KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE:
            return self._run_adaptive_prediction_lane(
                request=request,
                config=config,
                scoring_execution=scoring_execution,
            )
        raise WorkflowStageError(
            "kinase workflow internal invariant failed at seam="
            "kinase.executor.prediction_mode; "
            f"unsupported prediction mode: {config.prediction_mode}"
        )

    def _run_deterministic_prediction_lane(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        scoring_execution: KinaseScoringRunResult,
    ) -> KinasePredictionResult:
        dataset_sample_count = request.dataset._borrow_phospho_frame().shape[1]
        downstream_score_matrix = scoring_execution.downstream_score_matrix
        candidate_substrates = self._build_candidates(
            scores=downstream_score_matrix,
            top=config.prediction_top_k,
            score_threshold=CANDIDATE_SCORE_THRESHOLD,
            inclusion=CANDIDATE_MIN_INCLUSION,
        )
        kinase_ranking = self._rank_kinases(
            prediction_score_matrix=downstream_score_matrix,
            candidate_substrates=candidate_substrates,
        )
        selected_kinases = kinase_ranking.head(
            config.prediction_deterministic_max_selected_kinases
        ).index
        if selected_kinases.empty:
            candidate_shortfall = self._summarize_candidate_shortfall(
                scores=downstream_score_matrix,
                top=config.prediction_top_k,
                score_threshold=CANDIDATE_SCORE_THRESHOLD,
                inclusion=CANDIDATE_MIN_INCLUSION,
            )
            self._raise_boundary_error(
                seam="kinase.executor.prediction_ensemble",
                next_action=(
                    "provide dataset.phospho with at least two non-constant "
                    "sample columns or lower scoring_config.min_substrates "
                    "(scientific floor: min_substrates >= 2)"
                ),
                eligible_kinases=len(scoring_execution.quantified_substrates),
                ranked_kinases=int(kinase_ranking.size),
                prediction_config_deterministic_max_selected_kinases=(
                    config.prediction_deterministic_max_selected_kinases
                ),
                prediction_config_top_k=config.prediction_top_k,
                prediction_config_mode=config.prediction_mode,
                dataset_samples=dataset_sample_count,
                downstream_score_source=scoring_execution.downstream_score_source,
                candidate_qualifying_kinases=candidate_shortfall.qualifying_kinases,
                candidate_max_qualifying_sites=candidate_shortfall.max_qualifying_sites,
            )
        pred_mat, substrate_list = self._build_outputs(
            prediction_score_matrix=downstream_score_matrix,
            selected_kinases=selected_kinases,
            # Candidate keys and selected kinase labels are inherited from the
            # normalized downstream score matrix columns.
            candidate_substrates=candidate_substrates,
            top_k=config.prediction_top_k,
        )
        substrate_list = self._annotate_substrate_list(
            substrate_list=substrate_list,
            site_identity_map=_require_site_identity_map(request.site_identity_map),
        )
        return KinasePredictionResult._from_owned(
            pred_mat=pred_mat,
            substrate_list=substrate_list,
        )

    def _run_adaptive_prediction_lane(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        scoring_execution: KinaseScoringRunResult,
    ) -> KinasePredictionResult:
        dataset_sample_count = request.dataset._borrow_phospho_frame().shape[1]
        downstream_score_matrix = scoring_execution.downstream_score_matrix
        candidate_substrates = self._build_candidates(
            scores=downstream_score_matrix,
            top=config.prediction_top_k,
            score_threshold=CANDIDATE_SCORE_THRESHOLD,
            inclusion=CANDIDATE_MIN_INCLUSION,
        )
        if not candidate_substrates:
            candidate_shortfall = self._summarize_candidate_shortfall(
                scores=downstream_score_matrix,
                top=config.prediction_top_k,
                score_threshold=CANDIDATE_SCORE_THRESHOLD,
                inclusion=CANDIDATE_MIN_INCLUSION,
            )
            self._raise_boundary_error(
                seam="kinase.executor.prediction_adaptive_candidates",
                next_action=(
                    "provide dataset.phospho with at least two non-constant "
                    "sample columns or lower scoring_config.min_substrates "
                    "(scientific floor: min_substrates >= 2)"
                ),
                eligible_kinases=len(scoring_execution.quantified_substrates),
                candidate_kinases=0,
                prediction_config_mode=config.prediction_mode,
                prediction_config_top_k=config.prediction_top_k,
                prediction_config_adaptive_ensemble_runs=(
                    config.prediction_adaptive_ensemble_runs
                ),
                prediction_config_n_iterations=config.prediction_n_iterations,
                dataset_samples=dataset_sample_count,
                downstream_score_source=scoring_execution.downstream_score_source,
                candidate_qualifying_kinases=candidate_shortfall.qualifying_kinases,
                candidate_max_qualifying_sites=candidate_shortfall.max_qualifying_sites,
            )
        adaptive_random_state = self._require_adaptive_random_state(config)
        try:
            adaptive_scores = self._run_adaptive_prediction(
                prediction_score_matrix=downstream_score_matrix,
                candidate_substrates=candidate_substrates,
                prediction_config=self._as_prediction_config(
                    config,
                    random_state=adaptive_random_state,
                ),
                random_state=adaptive_random_state,
            )
        except ImportError as exc:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.prediction_adaptive_dependencies; "
                f"{exc}"
            ) from exc
        kinase_ranking = self._rank_kinases(
            prediction_score_matrix=adaptive_scores,
            candidate_substrates=candidate_substrates,
        )
        selected_kinases = kinase_ranking.index
        if selected_kinases.empty:
            self._raise_boundary_error(
                seam="kinase.executor.prediction_adaptive_ensemble",
                next_action=(
                    "lower prediction_config.top_k, increase dataset signal depth, or "
                    "review scoring-stage support for adaptive candidates"
                ),
                eligible_kinases=len(scoring_execution.quantified_substrates),
                candidate_kinases=len(candidate_substrates),
                ranked_kinases=0,
                prediction_config_mode=config.prediction_mode,
                prediction_config_top_k=config.prediction_top_k,
                prediction_config_adaptive_ensemble_runs=(
                    config.prediction_adaptive_ensemble_runs
                ),
                prediction_config_n_iterations=config.prediction_n_iterations,
            )
        pred_mat, substrate_list = self._build_outputs(
            prediction_score_matrix=adaptive_scores,
            selected_kinases=selected_kinases,
            # Candidate keys and selected kinase labels are inherited from the
            # normalized downstream score matrix columns.
            candidate_substrates=candidate_substrates,
            top_k=config.prediction_top_k,
            retain_full_scores=True,
        )
        substrate_list = self._annotate_substrate_list(
            substrate_list=substrate_list,
            site_identity_map=_require_site_identity_map(request.site_identity_map),
        )
        return KinasePredictionResult._from_owned(
            pred_mat=pred_mat,
            substrate_list=substrate_list,
        )

    @staticmethod
    def _annotate_substrate_list(
        *,
        substrate_list: pd.DataFrame,
        site_identity_map: pd.DataFrame,
    ) -> pd.DataFrame:
        if substrate_list.empty:
            annotated = substrate_list.copy(deep=True)
            annotated.loc[:, "site_key"] = pd.Series(dtype="object")
            annotated.loc[:, "display_id"] = pd.Series(dtype="object")
            return annotated
        display_lookup = {
            str(site_key): str(display_id)
            for site_key, display_id in site_identity_map.loc[
                :, ["site_key", "display_id"]
            ].itertuples(index=False)
        }
        annotated = substrate_list.copy(deep=True)
        site_keys = annotated.loc[:, "substrate_site"].astype(str)
        display_ids = site_keys.map(lambda value: display_lookup.get(value, value))
        annotated.loc[:, "site_key"] = site_keys
        annotated.loc[:, "display_id"] = display_ids
        annotated.loc[:, "substrate_site"] = display_ids
        return annotated

    @staticmethod
    def _as_prediction_config(
        config: ResolvedKinaseExecutionConfig,
        *,
        random_state: int,
    ) -> KinasePredictionConfig:
        return KinasePredictionConfig(
            top_k=config.prediction_top_k,
            deterministic_max_selected_kinases=(
                config.prediction_deterministic_max_selected_kinases
            ),
            adaptive_ensemble_runs=config.prediction_adaptive_ensemble_runs,
            mode=config.prediction_mode,
            adaptive_policy=config.prediction_adaptive_policy,
            n_iterations=config.prediction_n_iterations,
            random_state=random_state,
        )

    @staticmethod
    def _require_adaptive_random_state(config: ResolvedKinaseExecutionConfig) -> int:
        random_state = config.prediction_random_state
        if random_state is None:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.prediction_adaptive_random_state; "
                "prediction_config.random_state must be set when "
                "prediction_config.mode='adaptive_ensemble'"
            )
        return random_state

    @staticmethod
    def _raise_boundary_error(
        *,
        seam: str,
        next_action: str,
        **details: object,
    ) -> NoReturn:
        raise WorkflowBoundaryError(
            seam=seam,
            next_action=next_action,
            details=details,
            message_prefix="kinase workflow boundary validation failed",
        )


def _require_site_identity_map(site_identity_map: pd.DataFrame | None) -> pd.DataFrame:
    if site_identity_map is None:
        raise WorkflowBoundaryError(
            seam="kinase.prediction.site_identity_map",
            next_action="ensure kinase workflow interpretation resolves site identity mapping",
            message_prefix="kinase workflow boundary validation failed",
        )
    return site_identity_map


__all__ = ["KinasePredictionRunner"]
