"""Public kinase workflow shell."""

from __future__ import annotations

import pandas as pd

from phospy.activities.models import KinaseActivityResult
from phospy.api.requests import SimpleKinaseWorkflowRequest
from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.references.models import ReferenceBundle
from phospy.references.resolution import BundledReferenceProvider, ReferenceResolver


class SimpleKinaseWorkflow:
    """Public entrypoint for the kinase workflow."""

    def __init__(self) -> None:
        self._reference_resolver = ReferenceResolver(
            provider=BundledReferenceProvider()
        )

    def run(self, request: SimpleKinaseWorkflowRequest) -> SimpleKinaseWorkflowResult:
        references = self._resolve_references(request)
        scoring_result = KinaseScoringResult(profile_scores=pd.DataFrame())
        prediction_result = KinasePredictionResult(
            pred_mat=pd.DataFrame(),
            substrate_list=pd.DataFrame(),
        )
        activity_result: KinaseActivityResult | None = None
        if request.activity_config is not None and request.activity_config.enabled:
            activity_result = KinaseActivityResult(activity_scores=pd.DataFrame())
        return SimpleKinaseWorkflowResult(
            dataset=request.dataset,
            references=references,
            scoring_result=scoring_result,
            prediction_result=prediction_result,
            activity_result=activity_result,
        )

    def _resolve_references(
        self, request: SimpleKinaseWorkflowRequest
    ) -> ReferenceBundle:
        return self._reference_resolver.run(
            request.references,
            dataset_organism=request.dataset.organism,
        )
