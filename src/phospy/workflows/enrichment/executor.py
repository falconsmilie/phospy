"""Internal executor for enrichment workflow requests."""

from __future__ import annotations

from dataclasses import replace
from typing import cast

import pandas as pd

from phospy.contracts.results import EnrichmentWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance import (
    RunProvenance,
    collect_environment_provenance,
    fingerprint_table,
)
from phospy.science.enrichment.models import (
    MULTIPLE_TESTING_CORRECTION_NONE,
    EnrichmentResultRecord,
)
from phospy.science.enrichment.ora import OraEngine, OraResult, OraResultRecord
from phospy.science.statistics.multiple_testing import (
    MultipleTestingCorrection,
)
from phospy.science.statistics.multiple_testing import (
    run as run_multiple_testing_correction,
)
from phospy.workflows.enrichment.models import (
    InterpretedEnrichmentWorkflowRequest,
    MultipleTestingCorrectionRunner,
    OraEngineContract,
)

ENRICHMENT_OFFLINE_NO_ONLINE_RESOURCE_POLICY = "offline_user_supplied_collections_only"
ENRICHMENT_LIMITATIONS: tuple[str, ...] = (
    "offline over-representation analysis only",
    "gene-set or PTM-set collections must be supplied by the caller",
    "background universe is explicit and required; it is not inferred",
    "GO, KEGG, Reactome, and PTM-SEA resources are not bundled or fetched",
    "Enrichr, gseapy, and clusterProfiler online calls are not executed",
    "GSEA, ssGSEA, and PTM-SEA are not implemented by this workflow",
    "gene-level and site-level enrichment require explicit identifier semantics",
)


class EnrichmentWorkflowExecutor:
    """Run ORA and assemble a typed enrichment workflow result."""

    def __init__(
        self,
        *,
        ora_engine: OraEngineContract | None = None,
        correction_runner: MultipleTestingCorrectionRunner | None = None,
    ) -> None:
        self._ora_engine = ora_engine or OraEngine()
        self._correction_runner = correction_runner or run_multiple_testing_correction

    def run(
        self, request: InterpretedEnrichmentWorkflowRequest
    ) -> EnrichmentWorkflowResult:
        if not isinstance(cast(object, request), InterpretedEnrichmentWorkflowRequest):
            raise WorkflowBoundaryError(
                seam="enrichment.executor.interpreted_request_type",
                next_action=(
                    "pass interpreter output into EnrichmentWorkflowExecutor.run"
                ),
                message_prefix="enrichment workflow boundary validation failed",
            )

        raw_config = replace(
            request.method_config,
            multiple_testing_correction=MULTIPLE_TESTING_CORRECTION_NONE,
        )
        ora_result = self._ora_engine.run(
            selected_identifiers=request.selected_identifiers,
            background_universe=request.background_universe,
            enrichment_sets=request.set_collection,
            config=raw_config,
        )
        adjusted_p_values = self._correction_runner(
            tuple(record.p_value for record in ora_result.records),
            method=request.method_config.multiple_testing_correction,
        )
        if len(adjusted_p_values) != len(ora_result.records):
            raise WorkflowBoundaryError(
                seam="enrichment.executor.correction_result_length",
                next_action=(
                    "ensure the multiple-testing correction helper returns one "
                    "adjusted p-value per ORA record"
                ),
                details={
                    "record_count": len(ora_result.records),
                    "adjusted_count": len(adjusted_p_values),
                },
                message_prefix="enrichment workflow boundary validation failed",
            )

        records = _build_result_records(
            ora_records=ora_result.records,
            adjusted_p_values=adjusted_p_values,
            correction_method=request.method_config.multiple_testing_correction,
        )
        unmatched_identifiers = _unmatched_selected_identifiers(ora_result)
        background_summary = _execution_background_summary(
            request=request,
            ora_result=ora_result,
        )
        set_collection_summary = _execution_set_collection_summary(
            request=request,
            ora_result=ora_result,
        )
        diagnostics = _execution_diagnostics(
            request=request,
            ora_result=ora_result,
            unmatched_identifiers=unmatched_identifiers,
        )
        preliminary = EnrichmentWorkflowResult(
            identifier_kind=request.identifier_semantics.identifier_kind,
            set_collection=request.set_collection,
            config=request.config,
            records=records,
            unmatched_identifiers=unmatched_identifiers,
            diagnostics=diagnostics,
            method_metadata=request.method_metadata,
            background_summary=background_summary,
            set_collection_summary=set_collection_summary,
        )
        provenance = _build_run_provenance(
            request=request,
            result_table=preliminary.table,
        )
        return EnrichmentWorkflowResult(
            identifier_kind=request.identifier_semantics.identifier_kind,
            set_collection=request.set_collection,
            config=request.config,
            records=records,
            unmatched_identifiers=unmatched_identifiers,
            diagnostics=diagnostics,
            method_metadata=request.method_metadata,
            background_summary=background_summary,
            set_collection_summary=set_collection_summary,
            provenance=provenance,
        )


def _build_result_records(
    *,
    ora_records: tuple[OraResultRecord, ...],
    adjusted_p_values: tuple[float | None, ...],
    correction_method: MultipleTestingCorrection,
) -> tuple[EnrichmentResultRecord, ...]:
    records: list[EnrichmentResultRecord] = []
    for ora_record, adjusted_p_value in zip(
        ora_records,
        adjusted_p_values,
        strict=True,
    ):
        records.append(
            EnrichmentResultRecord(
                term_id=ora_record.set_id,
                term_name=ora_record.name,
                collection_kind=ora_record.collection_kind,
                identifier_kind=ora_record.identifier_kind,
                input_overlap_count=ora_record.overlap_size,
                background_overlap_count=ora_record.set_size,
                set_size=ora_record.raw_set_size,
                overlap_identifiers=ora_record.overlap_identifiers,
                p_value=ora_record.p_value,
                adjusted_p_value=adjusted_p_value,
                correction_method=correction_method,
                enrichment_ratio=ora_record.enrichment_ratio,
            )
        )
    return tuple(records)


def _unmatched_selected_identifiers(result: OraResult) -> tuple[str, ...]:
    matched = frozenset(
        identifier
        for record in result.records
        for identifier in record.overlap_identifiers
    )
    return tuple(
        identifier
        for identifier in result.selected_identifiers
        if identifier not in matched
    )


def _execution_background_summary(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
) -> dict[str, object]:
    summary = dict(request.background_summary)
    summary.update(
        {
            "universe_size": ora_result.background_size,
            "selected_in_background_count": ora_result.selected_size,
            "dropped_selected_count": len(ora_result.dropped_selected_identifiers),
            "dropped_selected_identifiers": ora_result.dropped_selected_identifiers,
        }
    )
    return summary


def _execution_set_collection_summary(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
) -> dict[str, object]:
    empty_after_background = tuple(
        record.set_id for record in ora_result.records if record.set_size == 0
    )
    summary = dict(request.set_collection_summary)
    summary.update(
        {
            "identifiers_outside_background_count": sum(
                record.set_identifiers_outside_background_count
                for record in ora_result.records
            ),
            "empty_after_background_count": len(empty_after_background),
            "empty_after_background_set_ids": empty_after_background,
        }
    )
    return summary


def _execution_diagnostics(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
    unmatched_identifiers: tuple[str, ...],
) -> dict[str, object]:
    diagnostics = dict(request.diagnostics)
    no_hit_set_ids = tuple(
        record.set_id for record in ora_result.records if record.overlap_size == 0
    )
    diagnostics.update(
        {
            "ora": {
                "method": ora_result.method,
                "record_count": len(ora_result.records),
                "selected_size": ora_result.selected_size,
                "background_size": ora_result.background_size,
                "no_hit_set_count": len(no_hit_set_ids),
                "no_hit_set_ids": no_hit_set_ids,
            },
            "multiple_testing_correction": {
                "method": request.method_config.multiple_testing_correction,
                "applied": True,
                "tested_record_count": len(ora_result.records),
            },
            "unmatched_selected_identifier_count": len(unmatched_identifiers),
        }
    )
    return diagnostics


def _build_run_provenance(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    result_table: pd.DataFrame,
) -> RunProvenance:
    workflow_parameters = {
        "method": request.config.method,
        "identifier_kind": request.identifier_semantics.identifier_kind,
        "identifier_column": request.identifier_semantics.identifier_column,
        "collection_kind": request.identifier_semantics.collection_kind,
        "analysis_level": request.identifier_semantics.analysis_level,
        "background_universe_source": "explicit",
        "background_universe_size": len(request.background_universe),
        "selected_identifier_count": len(request.selected_identifiers),
        "selected_identifier_source": request.selected_identifier_source,
        "multiple_testing_correction": (
            request.method_config.multiple_testing_correction
        ),
        "set_collection": _set_collection_provenance(request),
        "offline_no_online_resource_policy": (
            ENRICHMENT_OFFLINE_NO_ONLINE_RESOURCE_POLICY
        ),
        "online_resources_used": False,
        "limitations": ENRICHMENT_LIMITATIONS,
        "method_metadata": request.method_metadata,
        "background_summary": request.background_summary,
        "set_collection_summary": request.set_collection_summary,
    }
    return RunProvenance(
        environment=collect_environment_provenance(),
        input_tables=(
            fingerprint_table(
                _identifier_table(request.selected_identifiers),
                name="enrichment.selected_identifiers",
            ),
            fingerprint_table(
                _identifier_table(request.background_universe),
                name="enrichment.background_universe",
            ),
            fingerprint_table(
                _set_collection_table(request),
                name="enrichment.set_collection",
            ),
        ),
        preprocessing_stages=(),
        reference=None,
        workflow_name="enrichment",
        workflow_parameters=workflow_parameters,
        random_state=None,
        random_seed_policy=None,
        output_tables=(
            fingerprint_table(result_table, name="enrichment.result_table"),
        ),
    )


def _set_collection_provenance(
    request: InterpretedEnrichmentWorkflowRequest,
) -> dict[str, object]:
    return {
        "collection_kind": request.set_collection.collection_kind,
        "identifier_kind": request.set_collection.identifier_kind,
        "source_name": request.set_collection.source_name,
        "source_version": request.set_collection.source_version,
        "set_count": len(request.set_collection.enrichment_sets),
        "sets": tuple(
            {
                "set_id": enrichment_set.set_id,
                "name": enrichment_set.name,
                "source_name": enrichment_set.source_name,
                "source_version": enrichment_set.source_version,
            }
            for enrichment_set in request.set_collection.enrichment_sets
        ),
    }


def _identifier_table(identifiers: tuple[str, ...]) -> pd.DataFrame:
    return pd.DataFrame({"identifier": list(identifiers)})


def _set_collection_table(
    request: InterpretedEnrichmentWorkflowRequest,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for enrichment_set in request.set_collection.enrichment_sets:
        for identifier in enrichment_set.identifiers:
            rows.append(
                {
                    "set_id": enrichment_set.set_id,
                    "name": enrichment_set.name,
                    "identifier": identifier,
                    "identifier_kind": enrichment_set.identifier_kind,
                    "collection_kind": request.set_collection.collection_kind,
                    "source_name": enrichment_set.source_name,
                    "source_version": enrichment_set.source_version,
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "set_id",
            "name",
            "identifier",
            "identifier_kind",
            "collection_kind",
            "source_name",
            "source_version",
        ],
    )


__all__ = ["EnrichmentWorkflowExecutor"]
