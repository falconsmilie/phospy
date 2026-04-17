"""Public kinase workflow shell."""

from __future__ import annotations

import pandas as pd

from phospy.activities.models import KinaseActivityResult
from phospy.api.requests import SimpleKinaseWorkflowRequest
from phospy.api.results import SimpleKinaseWorkflowResult
from phospy.errors.references import PhosPyReferenceError
from phospy.errors.workflows import PhosPyWorkflowError
from phospy.prediction.models import KinasePredictionResult, KinaseScoringResult
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset


class SimpleKinaseWorkflow:
    """Public entrypoint for the kinase workflow."""

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

    @staticmethod
    def _resolve_references(request: SimpleKinaseWorkflowRequest) -> ReferenceBundle:
        reference_input = request.references
        if isinstance(reference_input, ReferenceBundle):
            return reference_input
        if not isinstance(reference_input, ReferencePreset):
            raise PhosPyWorkflowError("unsupported reference input type")
        organism = SimpleKinaseWorkflow._resolve_preset(reference_input, request)
        return ReferenceBundle(
            organism=organism,
            kinase_substrate_map=pd.DataFrame(),
            site_sequences=pd.DataFrame(),
        )

    @staticmethod
    def _resolve_preset(
        preset: ReferencePreset,
        request: SimpleKinaseWorkflowRequest,
    ) -> Organism:
        if preset is ReferencePreset.AUTO:
            if request.dataset.organism is None:
                raise PhosPyReferenceError(
                    "ReferencePreset.AUTO requires dataset.organism"
                )
            return request.dataset.organism
        mapping = {
            ReferencePreset.HUMAN: Organism.HUMAN,
            ReferencePreset.MOUSE: Organism.MOUSE,
            ReferencePreset.RAT: Organism.RAT,
        }
        organism = mapping[preset]
        if (
            request.dataset.organism is not None
            and request.dataset.organism != organism
        ):
            raise PhosPyReferenceError(
                "dataset.organism and requested reference preset must match"
            )
        return organism
