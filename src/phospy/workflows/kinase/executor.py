"""Internal executor for the simple kinase workflow."""

from __future__ import annotations

import pandas as pd

from phospy.activities.models import KinaseActivityResult
from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.errors.workflows import WorkflowStageError
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest


class SimpleKinaseWorkflowExecutor:
    """Run stage logic and assemble `SimpleKinaseWorkflowResult`."""

    _KINASE_COLUMN = "kinase"
    _SUBSTRATE_COLUMN = "substrate_site"

    def run(self, request: ResolvedKinaseWorkflowRequest) -> SimpleKinaseWorkflowResult:
        scoring_result, matched = self._run_scoring_stage(request)
        prediction_result = self._run_prediction_stage(
            request=request,
            matched_substrates=matched,
            scoring_result=scoring_result,
        )
        activity_result = self._run_activity_stage(
            request=request,
            prediction_result=prediction_result,
        )
        return SimpleKinaseWorkflowResult(
            dataset=request.dataset,
            references=request.references,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
        )

    def _run_scoring_stage(
        self, request: ResolvedKinaseWorkflowRequest
    ) -> tuple[KinaseScoringResult, pd.DataFrame]:
        mapping = self._validated_mapping(request.references.kinase_substrate_map)
        matched = mapping[
            mapping[self._SUBSTRATE_COLUMN].isin(request.dataset.phospho.index)
        ]
        support = (
            matched.groupby(self._KINASE_COLUMN)[self._SUBSTRATE_COLUMN]
            .nunique()
            .sort_values(ascending=False)
        )
        support = support[support >= request.scoring_config.min_substrates]

        profile_scores = support.rename("profile_score").to_frame()
        profile_scores.index.name = self._KINASE_COLUMN
        combined_scores = support.rename("combined_score").to_frame()
        combined_scores.index.name = self._KINASE_COLUMN

        if support.empty:
            weights = pd.DataFrame(columns=["weight"], index=profile_scores.index)
        else:
            weights = (support / support.sum()).rename("weight").to_frame()
        weights.index.name = self._KINASE_COLUMN

        return (
            KinaseScoringResult(
                profile_scores=profile_scores,
                motif_scores=None,
                combined_scores=combined_scores,
                weights=weights,
            ),
            matched,
        )

    def _run_prediction_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        matched_substrates: pd.DataFrame,
        scoring_result: KinaseScoringResult,
    ) -> KinasePredictionResult:
        selected_kinases = scoring_result.profile_scores.head(
            request.prediction_config.ensemble_size
        ).index
        pred_mat = pd.DataFrame(
            0.0,
            index=request.dataset.phospho.index.copy(),
            columns=selected_kinases.copy(),
        )
        pred_mat.index.name = request.dataset.phospho.index.name
        pred_mat.columns.name = self._KINASE_COLUMN

        if not pred_mat.empty:
            matched_pairs = (
                matched_substrates[
                    matched_substrates[self._KINASE_COLUMN].isin(selected_kinases)
                ][[self._KINASE_COLUMN, self._SUBSTRATE_COLUMN]]
                .drop_duplicates()
                .itertuples(index=False)
            )
            for kinase, substrate_site in matched_pairs:
                if substrate_site in pred_mat.index:
                    pred_mat.at[substrate_site, kinase] = 1.0

        substrate_list = (
            matched_substrates[
                matched_substrates[self._KINASE_COLUMN].isin(selected_kinases)
            ][[self._KINASE_COLUMN, self._SUBSTRATE_COLUMN]]
            .drop_duplicates()
            .head(request.prediction_config.top_k)
            .reset_index(drop=True)
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
        activity_scores = prediction_result.pred_mat.mean(axis=0).rename(
            "activity_score"
        )
        activity_table = activity_scores.to_frame()
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
        return mapping[[self._KINASE_COLUMN, self._SUBSTRATE_COLUMN]].copy(deep=True)
