"""Internal executor for the kinase workflow."""

from __future__ import annotations

from phospy.kinase.activity.runner import KinaseActivityRunner
from phospy.kinase.activity.validation import KinaseActivityInputValidator
from phospy.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.kinase.eligibility import KinaseEligibilityReportComposer
from phospy.kinase.prediction.candidates import build_candidate_substrate_list
from phospy.kinase.prediction.diagnostics import summarize_candidate_shortfall
from phospy.kinase.prediction.execution import run_adaptive_ensemble_prediction
from phospy.kinase.prediction.outputs import (
    build_prediction_outputs,
    rank_kinases_for_prediction,
)
from phospy.kinase.prediction.runner import KinasePredictionRunner
from phospy.kinase.provenance import KinaseProvenanceBuilder
from phospy.kinase.result_assembly import KinaseResultAssembler
from phospy.kinase.results import KinaseWorkflowResult
from phospy.kinase.scoring.execution import KinaseScoringRunner
from phospy.kinase.scoring.motif import (
    build_motif_library,
    get_motif_library_validation,
    score_phosphosite_motifs,
)
from phospy.kinase.scoring.profile import (
    build_kinase_profiles,
    score_profile_correlations,
)
from phospy.kinase.scoring.selection import select_downstream_score_matrix
from phospy.kinase.site_attrition import KinaseSiteAttritionSummaryComposer


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


__all__ = ["KinaseWorkflowExecutor"]
