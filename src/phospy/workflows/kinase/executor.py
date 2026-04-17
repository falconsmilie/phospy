"""Internal executor for the simple kinase workflow."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from phospy.activities.models import KinaseActivityResult
from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.errors.workflows import WorkflowStageError
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    build_prediction_outputs,
    rank_kinases_for_prediction,
    score_profile_correlations,
)


@dataclass(frozen=True, slots=True)
class _ScoringExecution:
    scoring_result: KinaseScoringResult
    combined_scores: pd.DataFrame
    quantified_substrates: dict[str, list[str]]


class SimpleKinaseWorkflowExecutor:
    """Run stage logic and assemble `SimpleKinaseWorkflowResult`."""

    _KINASE_COLUMN = "kinase"
    _SUBSTRATE_COLUMN = "substrate_site"

    def run(self, request: ResolvedKinaseWorkflowRequest) -> SimpleKinaseWorkflowResult:
        scoring_execution = self._run_scoring_stage(request)
        prediction_result = self._run_prediction_stage(
            request=request,
            scoring_execution=scoring_execution,
        )
        activity_result = self._run_activity_stage(
            request=request,
            prediction_result=prediction_result,
        )
        return SimpleKinaseWorkflowResult(
            dataset=request.dataset,
            references=request.references,
            scoring_result=scoring_execution.scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
        )

    def _run_scoring_stage(
        self, request: ResolvedKinaseWorkflowRequest
    ) -> _ScoringExecution:
        mapping = self._validated_mapping(request.references.kinase_substrate_map)
        profile_build = build_kinase_profiles(
            phospho=request.dataset.phospho,
            kinase_substrate_map=mapping,
            min_substrates=request.scoring_config.min_substrates,
        )
        if profile_build.profile_matrix.empty:
            raise WorkflowStageError(
                "no kinases have enough quantified substrates for scoring with the "
                "current dataset and scoring_config.min_substrates"
            )
        profile_scores = score_profile_correlations(
            phospho=request.dataset.phospho,
            profile_matrix=profile_build.profile_matrix,
        )
        combined_scores = profile_scores.copy(deep=True)
        weights = pd.DataFrame(
            {
                "motif_weight": 0.0,
                "profile_weight": 1.0,
            },
            index=combined_scores.columns.copy(),
        )
        weights.index.name = self._KINASE_COLUMN
        scoring_result = KinaseScoringResult(
            profile_scores=profile_scores,
            motif_scores=None,
            combined_scores=combined_scores,
            weights=weights,
        )
        return _ScoringExecution(
            scoring_result=scoring_result,
            combined_scores=combined_scores,
            quantified_substrates=profile_build.quantified_substrates,
        )

    def _run_prediction_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        scoring_execution: _ScoringExecution,
    ) -> KinasePredictionResult:
        kinase_ranking = rank_kinases_for_prediction(
            score_matrix=scoring_execution.combined_scores,
            quantified_substrates=scoring_execution.quantified_substrates,
        )
        selected_kinases = kinase_ranking.head(
            request.prediction_config.ensemble_size
        ).index
        pred_mat, substrate_list = build_prediction_outputs(
            score_matrix=scoring_execution.combined_scores,
            selected_kinases=selected_kinases,
            quantified_substrates=scoring_execution.quantified_substrates,
            top_k=request.prediction_config.top_k,
        )
        return KinasePredictionResult(
            pred_mat=pred_mat,
            substrate_list=substrate_list,
        )

    @staticmethod
    def _run_activity_stage(
        *,
        request: ResolvedKinaseWorkflowRequest,
        prediction_result: KinasePredictionResult,
    ) -> KinaseActivityResult | None:
        activity_config = request.activity_config
        if activity_config is None or not activity_config.enabled:
            return None
        pred_mat = prediction_result.pred_mat.astype(float)
        if pred_mat.empty:
            activity_table = pd.DataFrame(
                columns=[
                    "activity_score",
                    "weighted_signal",
                    "n_predicted_sites",
                    "is_active",
                ]
            )
            activity_table.index.name = "kinase"
            return KinaseActivityResult(activity_scores=activity_table)

        activity_scores = pred_mat.mean(axis=0).rename("activity_score")
        site_signal = request.dataset.phospho.astype(float).mean(axis=1)
        weighted_signal: dict[str, float] = {}
        predicted_counts: dict[str, int] = {}
        for kinase in pred_mat.columns:
            kinase_scores = pred_mat.loc[:, kinase].clip(lower=0.0)
            predicted_counts[str(kinase)] = int((kinase_scores > 0.0).sum())
            score_sum = float(kinase_scores.sum())
            if score_sum <= 0.0:
                weighted_signal[str(kinase)] = float("nan")
                continue
            weighted_signal[str(kinase)] = float(
                (kinase_scores * site_signal).sum() / score_sum
            )
        activity_table = activity_scores.to_frame()
        activity_table["weighted_signal"] = pd.Series(weighted_signal)
        activity_table["n_predicted_sites"] = pd.Series(predicted_counts).astype(
            "int64"
        )
        activity_table.index.name = "kinase"
        activity_table["is_active"] = (
            activity_table["activity_score"] >= activity_config.threshold
        )
        return KinaseActivityResult(activity_scores=activity_table)

    def _validated_mapping(self, mapping: pd.DataFrame) -> pd.DataFrame:
        required = {self._KINASE_COLUMN, self._SUBSTRATE_COLUMN}
        missing = required.difference(mapping.columns)
        if missing:
            missing_str = ", ".join(sorted(missing))
            raise WorkflowStageError(
                f"references.kinase_substrate_map is missing required columns: {missing_str}"
            )
        cleaned = mapping[[self._KINASE_COLUMN, self._SUBSTRATE_COLUMN]].copy(deep=True)
        if cleaned.isna().any(axis=None):
            raise WorkflowStageError(
                "references.kinase_substrate_map must not contain missing values"
            )
        cleaned.loc[:, self._KINASE_COLUMN] = (
            cleaned.loc[:, self._KINASE_COLUMN].astype(str).str.strip()
        )
        cleaned.loc[:, self._SUBSTRATE_COLUMN] = (
            cleaned.loc[:, self._SUBSTRATE_COLUMN].astype(str).str.strip()
        )
        if (cleaned.loc[:, self._KINASE_COLUMN] == "").any() or (
            cleaned.loc[:, self._SUBSTRATE_COLUMN] == ""
        ).any():
            raise WorkflowStageError(
                "references.kinase_substrate_map entries must be non-empty strings"
            )
        return cleaned.drop_duplicates(ignore_index=True)
