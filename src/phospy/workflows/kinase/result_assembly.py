"""Public result assembly for kinase workflow execution."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.results import (
    KinaseEligibilityReport,
    KinaseWorkflowResult,
    KinaseWorkflowSiteAttritionSummary,
)
from phospy.provenance.models import RunProvenance
from phospy.science.activities.models import KinaseActivityResult
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest


class KinaseResultAssembler:
    """Assemble the public kinase workflow result."""

    @staticmethod
    def run(
        *,
        request: ResolvedKinaseWorkflowRequest,
        scoring_result: KinaseScoringResult,
        prediction_result: KinasePredictionResult,
        eligibility_report: KinaseEligibilityReport | None,
        site_attrition_summary: KinaseWorkflowSiteAttritionSummary,
        activity_result: KinaseActivityResult | None,
        provenance: RunProvenance,
        substrate_contributions: pd.DataFrame | None = None,
    ) -> KinaseWorkflowResult:
        return KinaseWorkflowResult(
            dataset=request.dataset,
            references=request.references,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            eligibility_report=eligibility_report,
            site_attrition_summary=site_attrition_summary,
            activity_result=activity_result,
            provenance=provenance,
            substrate_contributions=substrate_contributions,
            _assume_owned=True,
        )


__all__ = ["KinaseResultAssembler"]
