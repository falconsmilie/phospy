"""Internal executor for the kinase workflow."""

from __future__ import annotations

from phospy.contracts.results import KinaseWorkflowResult
from phospy.science.activities.models import KinaseActivityResult
from phospy.science.prediction.candidates import (
    build_candidate_substrate_list,
    summarize_candidate_shortfall,
)
from phospy.science.prediction.execution import run_adaptive_ensemble_prediction
from phospy.science.prediction.models import KinasePredictionResult
from phospy.science.prediction.motif_scoring import (
    build_motif_library,
    get_motif_library_validation,
    score_phosphosite_motifs,
)
from phospy.science.prediction.scoring import (
    fuse_profile_and_motif_scores_by_rank_weight,
    select_downstream_score_matrix,
)
from phospy.validation.workflows.activity import KinaseActivityInputValidator
from phospy.workflows.kinase.activity_runner import KinaseActivityRunner
from phospy.workflows.kinase.component_models import KinaseScoringRunResult
from phospy.workflows.kinase.contracts import (
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.eligibility import KinaseEligibilityReportComposer
from phospy.workflows.kinase.prediction_runner import KinasePredictionRunner
from phospy.workflows.kinase.provenance import KinaseProvenanceBuilder
from phospy.workflows.kinase.result_assembly import KinaseResultAssembler
from phospy.workflows.kinase.science import (
    build_kinase_profiles,
    build_prediction_outputs,
    rank_kinases_for_prediction,
    score_profile_correlations,
)
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from phospy.workflows.kinase.site_attrition import KinaseSiteAttritionSummaryComposer


class KinaseWorkflowExecutor:
    """Run stage orchestration and assemble `KinaseWorkflowResult`."""

    def __init__(
        self,
        *,
        activity_input_validator: KinaseActivityInputValidator | None = None,
        scoring_runner: KinaseScoringRunner | None = None,
        prediction_runner: KinasePredictionRunner | None = None,
        activity_runner: KinaseActivityRunner | None = None,
        provenance_builder: KinaseProvenanceBuilder | None = None,
        eligibility_report_composer: KinaseEligibilityReportComposer | None = None,
        site_attrition_summary_composer: KinaseSiteAttritionSummaryComposer
        | None = None,
        result_assembler: KinaseResultAssembler | None = None,
    ) -> None:
        # Keep dependency wiring local so tests monkeypatching this module's symbols
        # can still intercept default runner behavior.
        self._scoring_runner = scoring_runner or KinaseScoringRunner(
            build_profiles=build_kinase_profiles,
            score_profiles=score_profile_correlations,
            build_motif_library_fn=build_motif_library,
            get_motif_library_validation_fn=get_motif_library_validation,
            score_motifs=score_phosphosite_motifs,
            fuse_scores=fuse_profile_and_motif_scores_by_rank_weight,
            select_downstream=select_downstream_score_matrix,
        )
        self._prediction_runner = prediction_runner or KinasePredictionRunner(
            build_candidates=build_candidate_substrate_list,
            summarize_candidate_shortfall_fn=summarize_candidate_shortfall,
            run_adaptive_prediction=run_adaptive_ensemble_prediction,
            rank_kinases=rank_kinases_for_prediction,
            build_outputs=build_prediction_outputs,
        )
        self._activity_runner = activity_runner or KinaseActivityRunner(
            activity_input_validator=activity_input_validator,
        )
        self._provenance_builder = provenance_builder or KinaseProvenanceBuilder()
        self._eligibility_report_composer = (
            eligibility_report_composer or KinaseEligibilityReportComposer()
        )
        self._site_attrition_summary_composer = (
            site_attrition_summary_composer or KinaseSiteAttritionSummaryComposer()
        )
        self._result_assembler = result_assembler or KinaseResultAssembler()

    def run(self, request: ResolvedKinaseWorkflowRequest) -> KinaseWorkflowResult:
        config = request.execution_config
        eligibility_report = self._eligibility_report_composer.run(
            request=request,
            config=config,
        )
        scoring_execution = self._scoring_runner.run(
            request=request,
            config=config,
            collect_substrate_contributions=bool(
                config.include_substrate_contributions
            ),
        )
        prediction_result = self._prediction_runner.run(
            request=request,
            config=config,
            scoring_execution=scoring_execution,
        )
        activity_result = self._activity_runner.run(
            request=request,
            config=config,
            prediction_result=prediction_result,
            scoring_execution=scoring_execution,
        )
        site_attrition_summary = self._site_attrition_summary_composer.run(
            request=request,
            scoring_execution=scoring_execution,
            prediction_result=prediction_result,
            activity_enabled=activity_result is not None,
        )
        provenance = self._provenance_builder.run(
            request=request,
            config=config,
            scoring_result=scoring_execution.scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
            substrate_contributions=(
                scoring_execution.substrate_contributions
                if config.include_substrate_contributions
                else None
            ),
        )
        return self._result_assembler.run(
            request=request,
            scoring_result=scoring_execution.scoring_result,
            prediction_result=prediction_result,
            eligibility_report=eligibility_report,
            site_attrition_summary=site_attrition_summary,
            activity_result=activity_result,
            provenance=provenance,
            substrate_contributions=(
                scoring_execution.substrate_contributions
                if config.include_substrate_contributions
                else None
            ),
        )

    # Compatibility hooks for existing tests that validate stage contracts directly.
    def _run_scoring_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
    ) -> KinaseScoringRunResult:
        return self._scoring_runner.run(
            request=request,
            config=config,
            collect_substrate_contributions=bool(
                config.include_substrate_contributions
            ),
        )

    def _run_prediction_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        scoring_execution: KinaseScoringRunResult,
    ) -> KinasePredictionResult:
        return self._prediction_runner.run(
            request=request,
            config=config,
            scoring_execution=scoring_execution,
        )

    def _run_activity_stage(
        self,
        *,
        request: ResolvedKinaseWorkflowRequest,
        config: ResolvedKinaseExecutionConfig,
        prediction_result: KinasePredictionResult,
    ) -> KinaseActivityResult | None:
        return self._activity_runner.run(
            request=request,
            config=config,
            prediction_result=prediction_result,
        )


__all__ = ["KinaseWorkflowExecutor"]
