"""Internal interpreter for simple kinase workflow requests."""

from __future__ import annotations

from phospy.api.requests import SimpleKinaseWorkflowRequest
from phospy.references.resolution import (
    BundledReferenceProvider,
    ReferenceResolver,
    ReferenceResolverContract,
)
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest


class SimpleKinaseWorkflowInterpreter:
    """Resolve workflow request defaults and references for execution."""

    def __init__(
        self, *, reference_resolver: ReferenceResolverContract | None = None
    ) -> None:
        self._reference_resolver = reference_resolver or ReferenceResolver(
            provider=BundledReferenceProvider()
        )

    def run(
        self, request: SimpleKinaseWorkflowRequest
    ) -> ResolvedKinaseWorkflowRequest:
        references = self._reference_resolver.run(
            request.references,
            dataset_organism=request.dataset.organism,
        )
        return ResolvedKinaseWorkflowRequest(
            dataset=request.dataset,
            references=references,
            scoring_config=request.scoring_config,
            prediction_config=request.prediction_config,
            activity_config=request.activity_config,
        )
