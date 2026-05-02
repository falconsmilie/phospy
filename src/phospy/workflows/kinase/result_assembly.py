"""Public result assembly for kinase workflow execution."""

from __future__ import annotations

from phospy.activities.models import KinaseActivityResult
from phospy.api.results import KinaseWorkflowResult, KinaseWorkflowSiteAttritionSummary
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.provenance.models import RunProvenance
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest


class KinaseResultAssembler:
    """Assemble the public kinase workflow result."""

    @staticmethod
    def run(
        *,
        request: ResolvedKinaseWorkflowRequest,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult,
        site_attrition_summary: KinaseWorkflowSiteAttritionSummary,
        activity_result: KinaseActivityResult | None,
        provenance: RunProvenance,
    ) -> KinaseWorkflowResult:
        return KinaseWorkflowResult(
            dataset=request.dataset,
            references=request.references,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            site_attrition_summary=site_attrition_summary,
            activity_result=activity_result,
            provenance=provenance,
        )


__all__ = ["KinaseResultAssembler"]
