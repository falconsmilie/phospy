"""Internal executor coordinator for enrichment workflow requests."""

from __future__ import annotations

from typing import cast

from phospy.contracts.results import EnrichmentWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.workflows.enrichment.models import (
    InterpretedEnrichmentWorkflowRequest,
    OraEngineContract,
)
from phospy.workflows.enrichment.ora_execution import EnrichmentOraRunner
from phospy.workflows.enrichment.provenance import EnrichmentRunProvenanceAssembler
from phospy.workflows.enrichment.result_assembly import EnrichmentResultAssembler
from phospy.workflows.enrichment.set_filtering import EnrichmentSetExecutionPreparer


class EnrichmentWorkflowExecutor:
    """Coordinate prepared ORA execution and public result assembly."""

    def __init__(
        self,
        *,
        ora_engine: OraEngineContract | None = None,
        set_preparer: EnrichmentSetExecutionPreparer | None = None,
        ora_runner: EnrichmentOraRunner | None = None,
        result_assembler: EnrichmentResultAssembler | None = None,
        provenance_assembler: EnrichmentRunProvenanceAssembler | None = None,
    ) -> None:
        self._set_preparer = set_preparer or EnrichmentSetExecutionPreparer()
        self._ora_runner = ora_runner or EnrichmentOraRunner(ora_engine=ora_engine)
        self._result_assembler = result_assembler or EnrichmentResultAssembler()
        self._provenance_assembler = (
            provenance_assembler or EnrichmentRunProvenanceAssembler()
        )

    def run(
        self,
        request: InterpretedEnrichmentWorkflowRequest,
    ) -> EnrichmentWorkflowResult:
        if not isinstance(cast(object, request), InterpretedEnrichmentWorkflowRequest):
            raise WorkflowBoundaryError(
                seam="enrichment.executor.interpreted_request_type",
                next_action=(
                    "pass interpreter output into EnrichmentWorkflowExecutor.run"
                ),
                message_prefix="enrichment workflow boundary validation failed",
            )

        set_size_filter_result = self._set_preparer.run(request)
        ora_result = self._ora_runner.run(
            request=request,
            set_size_filter_result=set_size_filter_result,
        )
        preliminary = self._result_assembler.run(
            request=request,
            ora_result=ora_result,
            set_size_filter_result=set_size_filter_result,
        )
        provenance = self._provenance_assembler.run(
            request=request,
            ora_result=ora_result,
            result_table=preliminary.table,
            set_size_filter_result=set_size_filter_result,
        )
        return self._result_assembler.run(
            request=request,
            ora_result=ora_result,
            set_size_filter_result=set_size_filter_result,
            provenance=provenance,
        )


__all__ = ["EnrichmentWorkflowExecutor"]
