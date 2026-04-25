"""Internal executor for the kinase workflow."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.activities.models import KinaseActivityResult
from phospy.activities.scoring import compute_activity_from_inputs
from phospy.api.configs import (
    KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE,
    KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING,
    KinasePredictionConfig,
)
from phospy.api.results import KinaseWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError, WorkflowStageError
from phospy.prediction.candidates import (
    build_candidate_substrate_list,
    summarize_candidate_shortfall,
)
from phospy.prediction.execution import run_adaptive_ensemble_prediction
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.prediction.motif_scoring import (
    DEFAULT_MOTIF_FLANK_SIZE,
    MotifScoringResult,
    build_motif_library,
    score_phosphosite_motifs,
)
from phospy.prediction.policies import resolve_prediction_sampling_policy
from phospy.prediction.scoring import (
    combine_profile_and_motif_scores,
    select_downstream_score_matrix,
)
from phospy.provenance.environment import collect_environment_provenance
from phospy.provenance.hashing import fingerprint_optional_table
from phospy.provenance.models import (
    PreprocessingStageProvenance,
    RunProvenance,
    TableFingerprint,
)
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    build_prediction_outputs,
    rank_kinases_for_prediction,
    score_profile_correlations,
)


@dataclass(frozen=True, slots=True)
class _ScoringExecution:
    scoring_result: KinaseScoringResult
    downstream_score_matrix: pd.DataFrame
    downstream_score_source: str
    quantified_substrates: dict[str, list[str]]


class KinaseWorkflowExecutor:
    """Run stage logic and assemble `KinaseWorkflowResult`."""

    def __init__(
        self,
        *,
        activity_input_validator: KinaseActivityInputValidator | None = None,
    ) -> None:
        self._activity_input_validator = (
            activity_input_validator or KinaseActivityInputValidator()
        )

    def run(self, request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
        config = request.execution_config
        scoring_execution = self._run_scoring_stage(
            request=request,
            config=config,
        )
        prediction_result = self._run_prediction_stage(
            request=request,
            config=config,
            scoring_execution=scoring_execution,
        )
        activity_result = self._run_activity_stage(
            request=request,
            config=config,
            prediction_result=prediction_result,
        )
        provenance = _build_kinase_run_provenance(
            request=request,
            config=config,
            scoring_result=scoring_execution.scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
        )
        return KinaseWorkflowResult(
            dataset=request.dataset,
            references=request.references,
            scoring_result=scoring_execution.scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
            provenance=provenance,
        )

    def _run_scoring_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
    ) -> _ScoringExecution:
        # Authoritative downstream route:
        # - profile correlations from quantified kinase substrates
        # - profile+motif combined scores (with profile fallback)
        #
        # Optional diagnostic tables:
        # - motif_scores
        # - weights
        include_diagnostic_tables = config.include_diagnostic_scoring_tables
        scoring_min_substrates = int(config.scoring_min_substrates)
        scoring_phospho = request.dataset.phospho.loc[request.scoring_site_index, :]
        profile_build = build_kinase_profiles(
            phospho=scoring_phospho,
            kinase_substrate_map=request.kinase_substrate_map,
            min_substrates=scoring_min_substrates,
            allow_single_substrate_profiles=False,
            profile_missing_value_strategy=config.profile_missing_value_strategy,
        )
        if profile_build.profile_matrix.empty:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.scoring_profiles; interpreter should reject "
                "requests with zero eligible kinases before scoring"
            )
        profile_scores = score_profile_correlations(
            phospho=scoring_phospho,
            profile_matrix=profile_build.profile_matrix,
        )
        eligible_kinases = set(profile_scores.columns.astype(str))
        motif_kinase_substrate_map = request.kinase_substrate_map.loc[
            request.kinase_substrate_map.loc[:, "kinase"]
            .astype(str)
            .isin(eligible_kinases)
        ]
        sequence_series = request.site_sequences.loc[:, "site_sequence"]
        motif_result = self._resolve_motif_scores(
            request=request,
            scoring_min_substrates=scoring_min_substrates,
            scoring_phospho=scoring_phospho,
            motif_kinase_substrate_map=motif_kinase_substrate_map,
            sequence_series=sequence_series,
        )
        try:
            combined_scores, weights = combine_profile_and_motif_scores(
                motif_scores=motif_result.motif_scores,
                profile_scores=profile_scores,
                motif_sizes=motif_result.motif_sizes,
                profile_sizes=profile_build.substrate_counts.astype(float),
                allow_profile_only_fallback=True,
                emit_weights=include_diagnostic_tables,
            )
        except ValueError as exc:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.combined_scoring; "
                f"{exc}"
            ) from exc
        scoring_result = KinaseScoringResult._from_owned(
            profile_scores=profile_scores,
            motif_scores=(
                motif_result.motif_scores if include_diagnostic_tables else None
            ),
            combined_scores=combined_scores,
            weights=weights,
        )
        downstream_score_matrix, downstream_score_source = (
            select_downstream_score_matrix(
                profile_scores=profile_scores,
                combined_scores=combined_scores,
            )
        )
        return _ScoringExecution(
            scoring_result=scoring_result,
            downstream_score_matrix=downstream_score_matrix,
            downstream_score_source=downstream_score_source,
            quantified_substrates=profile_build.quantified_substrates,
        )

    @staticmethod
    def _resolve_motif_scores(
        *,
        request: ResolvedKinaseWorkflowRequest,
        scoring_min_substrates: int,
        scoring_phospho: pd.DataFrame,
        motif_kinase_substrate_map: pd.DataFrame,
        sequence_series: pd.Series,
    ) -> MotifScoringResult:
        motif_frequency_matrices, motif_sizes = build_motif_library(
            kinase_substrate_map=motif_kinase_substrate_map,
            site_sequences=sequence_series,
            flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        )
        return score_phosphosite_motifs(
            site_sequences=sequence_series.loc[scoring_phospho.index],
            motif_frequency_matrices=motif_frequency_matrices,
            motif_sizes=motif_sizes,
            site_index=scoring_phospho.index,
            min_motif_size=scoring_min_substrates,
            flank_size=DEFAULT_MOTIF_FLANK_SIZE,
        )

    def _run_prediction_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        scoring_execution: _ScoringExecution,
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
        scoring_execution: _ScoringExecution,
    ) -> KinasePredictionResult:
        downstream_score_matrix = scoring_execution.downstream_score_matrix
        candidate_substrates = build_candidate_substrate_list(
            scores=downstream_score_matrix,
            top=config.prediction_top_k,
            score_threshold=0.0,
            inclusion=1,
        )
        kinase_ranking = rank_kinases_for_prediction(
            prediction_score_matrix=downstream_score_matrix,
            candidate_substrates=candidate_substrates,
        )
        selected_kinases = kinase_ranking.head(config.prediction_ensemble_size).index
        if selected_kinases.empty:
            candidate_shortfall = summarize_candidate_shortfall(
                scores=downstream_score_matrix,
                top=config.prediction_top_k,
                score_threshold=0.0,
                inclusion=1,
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
                prediction_config_ensemble_size=config.prediction_ensemble_size,
                prediction_config_top_k=config.prediction_top_k,
                prediction_config_mode=config.prediction_mode,
                dataset_samples=request.dataset.phospho.shape[1],
                downstream_score_source=scoring_execution.downstream_score_source,
                candidate_qualifying_kinases=candidate_shortfall.qualifying_kinases,
                candidate_max_qualifying_sites=candidate_shortfall.max_qualifying_sites,
            )
        pred_mat, substrate_list = build_prediction_outputs(
            prediction_score_matrix=downstream_score_matrix,
            selected_kinases=selected_kinases,
            candidate_substrates=candidate_substrates,
            top_k=config.prediction_top_k,
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
        scoring_execution: _ScoringExecution,
    ) -> KinasePredictionResult:
        downstream_score_matrix = scoring_execution.downstream_score_matrix
        candidate_substrates = build_candidate_substrate_list(
            scores=downstream_score_matrix,
            top=config.prediction_top_k,
            score_threshold=0.0,
            inclusion=1,
        )
        if not candidate_substrates:
            candidate_shortfall = summarize_candidate_shortfall(
                scores=downstream_score_matrix,
                top=config.prediction_top_k,
                score_threshold=0.0,
                inclusion=1,
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
                prediction_config_ensemble_size=config.prediction_ensemble_size,
                prediction_config_n_iterations=config.prediction_n_iterations,
                dataset_samples=request.dataset.phospho.shape[1],
                downstream_score_source=scoring_execution.downstream_score_source,
                candidate_qualifying_kinases=candidate_shortfall.qualifying_kinases,
                candidate_max_qualifying_sites=candidate_shortfall.max_qualifying_sites,
            )
        try:
            adaptive_scores = run_adaptive_ensemble_prediction(
                prediction_score_matrix=downstream_score_matrix,
                candidate_substrates=candidate_substrates,
                prediction_config=self._as_prediction_config(config),
            )
        except ImportError as exc:
            raise WorkflowStageError(
                "kinase workflow internal invariant failed at seam="
                "kinase.executor.prediction_adaptive_dependencies; "
                f"{exc}"
            ) from exc
        kinase_ranking = rank_kinases_for_prediction(
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
                prediction_config_ensemble_size=config.prediction_ensemble_size,
                prediction_config_n_iterations=config.prediction_n_iterations,
            )
        pred_mat, substrate_list = build_prediction_outputs(
            prediction_score_matrix=adaptive_scores,
            selected_kinases=selected_kinases,
            candidate_substrates=candidate_substrates,
            top_k=config.prediction_top_k,
            retain_full_scores=True,
        )
        return KinasePredictionResult._from_owned(
            pred_mat=pred_mat,
            substrate_list=substrate_list,
        )

    def _run_activity_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        prediction_result: KinasePredictionResult,
    ) -> KinaseActivityResult | None:
        activity_config = config.activity
        if activity_config is None:
            return None
        validated_inputs = self._activity_input_validator.run(
            pred_mat=prediction_result.pred_mat,
            phospho_matrix=request.activity_phospho_matrix,
            threshold=activity_config.threshold,
            min_substrates=activity_config.min_substrates,
            top_n_substrates=activity_config.top_n_substrates,
        )
        return compute_activity_from_inputs(validated_inputs)

    @staticmethod
    def _as_prediction_config(
        config: ResolvedKinaseExecutionConfig,
    ) -> KinasePredictionConfig:
        return KinasePredictionConfig(
            top_k=config.prediction_top_k,
            ensemble_size=config.prediction_ensemble_size,
            mode=config.prediction_mode,
            adaptive_policy=config.prediction_adaptive_policy,
            n_iterations=config.prediction_n_iterations,
            random_state=config.prediction_random_state,
        )

    @staticmethod
    def _raise_boundary_error(
        *,
        seam: str,
        next_action: str,
        **details: object,
    ) -> None:
        raise WorkflowBoundaryError(
            seam=seam,
            next_action=next_action,
            details=details,
            message_prefix="kinase workflow boundary validation failed",
        )


def _build_kinase_run_provenance(
    *,
    request: ResolvedKinaseWorkflowRequest,
    config: ResolvedKinaseExecutionConfig,
    scoring_result: KinaseScoringResult,
    prediction_result: KinasePredictionResult,
    activity_result: KinaseActivityResult | None,
) -> RunProvenance:
    input_tables = _collect_fingerprints(
        (
            ("dataset.phospho", request.dataset.phospho),
            ("dataset.site_metadata", request.dataset.site_metadata),
            ("dataset.sample_metadata", request.dataset.sample_metadata),
            ("dataset.total", request.dataset.total),
            ("dataset.comparisons", request.dataset.comparisons),
            (
                "references.kinase_substrate_map",
                request.references.kinase_substrate_map,
            ),
            ("references.site_sequences", request.references.site_sequences),
        )
    )
    output_tables = _collect_fingerprints(
        (
            ("outputs.scoring.profile_scores", scoring_result.profile_scores),
            ("outputs.scoring.motif_scores", scoring_result.motif_scores),
            ("outputs.scoring.combined_scores", scoring_result.combined_scores),
            ("outputs.scoring.weights", scoring_result.weights),
            ("outputs.prediction.pred_mat", prediction_result.pred_mat),
            ("outputs.prediction.substrate_list", prediction_result.substrate_list),
            (
                "outputs.activity.weighted_activity",
                None if activity_result is None else activity_result.weighted_activity,
            ),
            (
                "outputs.activity.ksea_scores",
                None if activity_result is None else activity_result.ksea_scores,
            ),
            (
                "outputs.activity.ksea_counts",
                (
                    None
                    if activity_result is None
                    else activity_result.ksea_counts.to_frame(name="n_substrates")
                ),
            ),
            (
                "outputs.activity.target_counts",
                (
                    None
                    if activity_result is None
                    else activity_result.target_counts.to_frame(name="n_targets")
                ),
            ),
            (
                "outputs.activity.target_table",
                None if activity_result is None else activity_result.target_table,
            ),
        )
    )
    return RunProvenance(
        environment=collect_environment_provenance(),
        input_tables=input_tables,
        preprocessing_stages=_dataset_preprocessing_stages(request),
        reference=request.references.provenance,
        workflow_name="kinase_workflow",
        workflow_parameters={
            "scoring_config": {
                "min_substrates": int(config.scoring_min_substrates),
                "include_diagnostic_scoring_tables": bool(
                    config.include_diagnostic_scoring_tables
                ),
                "profile_missing_value_strategy": str(
                    config.profile_missing_value_strategy
                ),
            },
            "prediction_config": {
                "top_k": int(config.prediction_top_k),
                "ensemble_size": int(config.prediction_ensemble_size),
                "mode": str(config.prediction_mode),
                "adaptive_policy": str(config.prediction_adaptive_policy),
                "n_iterations": int(config.prediction_n_iterations),
                "random_state": config.prediction_random_state,
            },
            "activity_config": (
                None
                if config.activity is None
                else {
                    "threshold": float(config.activity.threshold),
                    "min_substrates": int(config.activity.min_substrates),
                    "top_n_substrates": int(config.activity.top_n_substrates),
                }
            ),
        },
        random_state=config.prediction_random_state,
        random_seed_policy=_resolve_seed_policy(config),
        output_tables=output_tables,
    )


def _resolve_seed_policy(config: ResolvedKinaseExecutionConfig) -> str | None:
    if config.prediction_mode == KINASE_PREDICTION_MODE_DETERMINISTIC_RANKING:
        return None
    if config.prediction_mode != KINASE_PREDICTION_MODE_ADAPTIVE_ENSEMBLE:
        return None
    return resolve_prediction_sampling_policy(
        config.prediction_adaptive_policy
    ).seed_strategy


def _dataset_preprocessing_stages(
    request: ResolvedKinaseWorkflowRequest,
) -> tuple[PreprocessingStageProvenance, ...]:
    provenance = request.dataset.provenance
    if provenance is None:
        return ()
    return tuple(provenance.preprocessing_stages)


def _collect_fingerprints(
    entries: tuple[tuple[str, pd.DataFrame | None], ...],
) -> tuple[TableFingerprint, ...]:
    fingerprints: list[TableFingerprint] = []
    for name, table in entries:
        fingerprint = fingerprint_optional_table(table, name=name)
        if fingerprint is None:
            continue
        fingerprints.append(fingerprint)
    return tuple(fingerprints)
