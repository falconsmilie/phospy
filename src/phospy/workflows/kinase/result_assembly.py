"""Public result assembly for kinase workflow execution."""

from __future__ import annotations

from collections.abc import Mapping
from typing import cast

import pandas as pd

from phospy.contracts.results import (
    KinaseEligibilityReport,
    KinaseWorkflowAttritionProvenance,
    KinaseWorkflowCaveat,
    KinaseWorkflowResult,
    KinaseWorkflowSiteAttritionSummary,
)
from phospy.provenance.models import RunProvenance
from phospy.science.activities.models import KinaseActivityResult
from phospy.science.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.workflows.kinase.attrition_metrics import (
    build_kinase_attrition_provenance_payload,
)
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
        attrition_metrics = request.attrition_metrics
        if attrition_metrics is None:
            raise RuntimeError(
                "kinase result assembly requires resolved attrition metrics"
            )
        attrition_payload = build_kinase_attrition_provenance_payload(
            metrics=attrition_metrics,
            policy=request.execution_config.attrition_policy,
            violations=request.attrition_policy_violations,
        )
        policy_violations = cast(
            list[Mapping[str, object]],
            attrition_payload["policy_violations"],
        )
        warning_messages = cast(
            list[str],
            attrition_payload["warning_messages"],
        )
        return KinaseWorkflowResult(
            dataset=request.dataset,
            references=request.references,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            eligibility_report=eligibility_report,
            site_attrition_summary=site_attrition_summary,
            attrition_provenance=KinaseWorkflowAttritionProvenance(
                metrics=cast(Mapping[str, object], attrition_payload["metrics"]),
                policy=cast(Mapping[str, object], attrition_payload["policy"]),
                policy_outcome=str(attrition_payload["policy_outcome"]),
                policy_violations=tuple(policy_violations),
                warning_messages=tuple(warning_messages),
            ),
            activity_result=activity_result,
            provenance=provenance,
            substrate_contributions=substrate_contributions,
            caveats=tuple(
                KinaseWorkflowCaveat(
                    code=str(violation.to_payload()["code"]),
                    severity="warning",
                    message=violation.message,
                    details=violation.to_payload(),
                )
                for violation in request.attrition_policy_violations
            ),
            _assume_owned=True,
        )


__all__ = ["KinaseResultAssembler"]
