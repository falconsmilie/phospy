"""Public signalome workflow shell."""

from __future__ import annotations

import pandas as pd

from phospy.api.requests import SignalomeWorkflowRequest
from phospy.api.results import SignalomeWorkflowResult
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)


class SignalomeWorkflow:
    """Public entrypoint for the signalome workflow."""

    def run(self, request: SignalomeWorkflowRequest) -> SignalomeWorkflowResult:
        dataset = request.kinase_result.dataset
        return SignalomeWorkflowResult(
            dataset=dataset,
            kinase_result=request.kinase_result,
            module_assignments=SignalomeAssignments(table=pd.DataFrame()),
            signalome_modules=SignalomeModules(table=pd.DataFrame()),
            kinase_network=KinaseNetwork(edges=pd.DataFrame()),
            expanded_signalome=None,
        )
