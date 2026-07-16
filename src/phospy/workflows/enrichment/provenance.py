"""Enrichment workflow run provenance assembly."""

from __future__ import annotations

import pandas as pd

from phospy.contracts.enrichment_identifier_sets import (
    EnrichmentIdentifierSetProvenance,
)
from phospy.provenance import (
    RowAttritionRecord,
    RowAttritionReport,
    RunProvenance,
    collect_environment_provenance,
    fingerprint_table,
)
from phospy.science.enrichment.ora import OraResult
from phospy.workflows.enrichment.models import InterpretedEnrichmentWorkflowRequest
from phospy.workflows.enrichment.result_assembly import (
    retained_foreground_fraction,
    unmatched_set_identifier_summary_for_request,
)
from phospy.workflows.enrichment.set_filtering import (
    SET_SIZE_FILTER_REASON_ABOVE_MAX,
    SET_SIZE_FILTER_REASON_BELOW_MIN,
    SetSizeFilterResult,
    dropped_set_reason_counts,
    set_size_filters_configured,
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


class EnrichmentRunProvenanceAssembler:
    """Assemble deterministic run provenance for enrichment workflow results."""

    def run(
        self,
        *,
        request: InterpretedEnrichmentWorkflowRequest,
        ora_result: OraResult,
        result_table: pd.DataFrame,
        set_size_filter_result: SetSizeFilterResult,
    ) -> RunProvenance:
        multiple_testing_method = ora_result.config.multiple_testing_correction
        number_of_tests = len(ora_result.records)
        row_attrition_metrics = _row_attrition_metrics(
            request=request,
            ora_result=ora_result,
            set_size_filter_result=set_size_filter_result,
        )
        row_attrition = _row_attrition_reports(
            request=request,
            ora_result=ora_result,
            set_size_filter_result=set_size_filter_result,
        )
        workflow_parameters = {
            "method": request.config.method,
            "identifier_kind": request.identifier_semantics.identifier_kind,
            "identifier_column": request.identifier_semantics.identifier_column,
            "collection_kind": request.identifier_semantics.collection_kind,
            "analysis_level": request.identifier_semantics.analysis_level,
            "background_universe_source": "explicit",
            "background_identifiers_provided": request.background_identifier_input_count,
            "background_universe_size": len(request.background_universe),
            "background_identifiers_retained_in_universe": len(
                request.background_universe
            ),
            "selected_identifiers_provided": request.selected_identifier_input_count,
            "selected_identifier_count": len(request.selected_identifiers),
            "selected_identifiers_retained_in_universe": ora_result.selected_size,
            "retained_foreground_fraction": retained_foreground_fraction(
                retained_count=ora_result.selected_size,
                selected_count=len(request.selected_identifiers),
            ),
            "selected_identifier_source": request.selected_identifier_source,
            "selected_outside_background_policy": (
                request.method_config.selected_outside_background_policy
            ),
            "set_member_outside_background_policy": (
                request.config.set_member_outside_background_policy
            ),
            "minimum_retained_foreground_fraction": (
                request.config.minimum_retained_foreground_fraction
            ),
            "multiple_testing_correction": (
                request.method_config.multiple_testing_correction
            ),
            "multiple_testing_method": multiple_testing_method,
            "number_of_tests": number_of_tests,
            "correction_owner": "ora_engine",
            "set_collection": _set_collection_provenance(request),
            "selected_identifier_provenance": _identifier_set_provenance_payload(
                request.selected_identifier_provenance
            ),
            "background_identifier_provenance": _identifier_set_provenance_payload(
                request.background_identifier_provenance
            ),
            "offline_no_online_resource_policy": (
                ENRICHMENT_OFFLINE_NO_ONLINE_RESOURCE_POLICY
            ),
            "online_resources_used": False,
            "limitations": ENRICHMENT_LIMITATIONS,
            "method_metadata": request.method_metadata,
            "background_summary": request.background_summary,
            "set_collection_summary": request.set_collection_summary,
            "universe_policy": _universe_policy_provenance(
                request=request,
                ora_result=ora_result,
            ),
            "row_attrition_metrics": row_attrition_metrics,
        }
        if row_attrition:
            workflow_parameters["row_attrition"] = row_attrition
        if set_size_filters_configured(request):
            workflow_parameters["set_size_filter"] = {
                "min_set_size": request.config.min_set_size,
                "max_set_size": request.config.max_set_size,
                "applied_after_background_intersection": True,
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


def _row_attrition_metrics(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
    set_size_filter_result: SetSizeFilterResult,
) -> dict[str, object]:
    dropped_reason_counts = dropped_set_reason_counts(
        set_size_filter_result.dropped_sets
    )
    return {
        "selected_identifiers_provided": int(request.selected_identifier_input_count),
        "selected_identifiers_prepared": int(len(request.selected_identifiers)),
        "selected_identifiers_retained_in_universe": int(ora_result.selected_size),
        "retained_foreground_fraction": retained_foreground_fraction(
            retained_count=ora_result.selected_size,
            selected_count=len(request.selected_identifiers),
        ),
        "selected_identifiers_dropped_before_universe_intersection": int(
            max(
                request.selected_identifier_input_count
                - len(request.selected_identifiers),
                0,
            )
        ),
        "selected_identifiers_dropped_outside_universe": int(
            len(ora_result.dropped_selected_identifiers)
        ),
        "background_identifiers_provided": int(
            request.background_identifier_input_count
        ),
        "background_identifiers_retained_in_universe": int(
            len(request.background_universe)
        ),
        "background_identifiers_dropped_before_universe_use": int(
            max(
                request.background_identifier_input_count
                - len(request.background_universe),
                0,
            )
        ),
        "sets_provided": int(set_size_filter_result.input_set_count),
        "sets_tested": int(set_size_filter_result.tested_set_count),
        "sets_skipped_due_to_min_max_size": int(
            len(set_size_filter_result.dropped_sets)
        ),
        "sets_skipped_due_to_min_size": int(
            dropped_reason_counts[SET_SIZE_FILTER_REASON_BELOW_MIN]
        ),
        "sets_skipped_due_to_max_size": int(
            dropped_reason_counts[SET_SIZE_FILTER_REASON_ABOVE_MAX]
        ),
    }


def _row_attrition_reports(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
    set_size_filter_result: SetSizeFilterResult,
) -> dict[str, object]:
    reports: dict[str, object] = {}
    selected_report = _selected_identifier_attrition_report(
        request=request,
        ora_result=ora_result,
    )
    if selected_report is not None:
        reports["selected_identifiers"] = selected_report.to_payload()
    background_report = _background_identifier_attrition_report(request)
    if background_report is not None:
        reports["background_identifiers"] = background_report.to_payload()
    set_report = _set_attrition_report(set_size_filter_result)
    if set_report is not None:
        reports["sets"] = set_report.to_payload()
    return reports


def _selected_identifier_attrition_report(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
) -> RowAttritionReport | None:
    records: list[RowAttritionRecord] = []
    current_rows = int(request.selected_identifier_input_count)
    prepared_count = int(len(request.selected_identifiers))
    normalized_drop_count = max(current_rows - prepared_count, 0)
    if normalized_drop_count:
        current_rows = _append_row_attrition_record(
            records,
            stage="enrichment_selected_identifier_validation",
            input_rows=current_rows,
            removed_rows=normalized_drop_count,
            reason="selected_identifiers_dropped_before_universe_intersection",
            examples=(),
        )
    outside_background_count = int(len(ora_result.dropped_selected_identifiers))
    if outside_background_count:
        current_rows = _append_row_attrition_record(
            records,
            stage="enrichment_selected_identifier_universe_intersection",
            input_rows=current_rows,
            removed_rows=outside_background_count,
            reason="selected_identifiers_not_retained_in_universe",
            examples=_examples(ora_result.dropped_selected_identifiers),
        )
    if not records:
        return None
    _ = current_rows
    return RowAttritionReport.from_records(records)


def _background_identifier_attrition_report(
    request: InterpretedEnrichmentWorkflowRequest,
) -> RowAttritionReport | None:
    input_rows = int(request.background_identifier_input_count)
    retained_rows = int(len(request.background_universe))
    removed_rows = max(input_rows - retained_rows, 0)
    if removed_rows <= 0:
        return None
    return RowAttritionReport.from_records(
        (
            RowAttritionRecord(
                stage="enrichment_background_identifier_validation",
                input_rows=input_rows,
                output_rows=retained_rows,
                removed_rows=removed_rows,
                reason="background_identifiers_dropped_before_universe_use",
                examples=(),
            ),
        )
    )


def _set_attrition_report(
    set_size_filter_result: SetSizeFilterResult,
) -> RowAttritionReport | None:
    removed_rows = int(len(set_size_filter_result.dropped_sets))
    if removed_rows <= 0:
        return None
    return RowAttritionReport.from_records(
        (
            RowAttritionRecord(
                stage="enrichment_set_size_filter",
                input_rows=int(set_size_filter_result.input_set_count),
                output_rows=int(set_size_filter_result.tested_set_count),
                removed_rows=removed_rows,
                reason="sets_skipped_due_to_min_max_size",
                examples=_examples(
                    tuple(
                        dropped_set.set_id
                        for dropped_set in set_size_filter_result.dropped_sets
                    )
                ),
            ),
        )
    )


def _append_row_attrition_record(
    records: list[RowAttritionRecord],
    *,
    stage: str,
    input_rows: int,
    removed_rows: int,
    reason: str,
    examples: tuple[str, ...],
) -> int:
    output_rows = int(input_rows) - int(removed_rows)
    records.append(
        RowAttritionRecord(
            stage=stage,
            input_rows=int(input_rows),
            output_rows=output_rows,
            removed_rows=int(removed_rows),
            reason=reason,
            examples=examples,
        )
    )
    return output_rows


def _examples(values: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(str(value) for value in values[:5])


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


def _universe_policy_provenance(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
) -> dict[str, object]:
    set_identifier_summary = unmatched_set_identifier_summary_for_request(
        request=request
    )
    return {
        "selected_outside_background_policy": (
            request.method_config.selected_outside_background_policy
        ),
        "set_member_outside_background_policy": (
            request.config.set_member_outside_background_policy
        ),
        "minimum_retained_foreground_fraction": (
            request.config.minimum_retained_foreground_fraction
        ),
        "selected_identifier_count": len(request.selected_identifiers),
        "selected_identifiers_retained_in_background_count": ora_result.selected_size,
        "selected_identifiers_outside_background_count": len(
            ora_result.dropped_selected_identifiers
        ),
        "selected_identifiers_outside_background": (
            ora_result.dropped_selected_identifiers
        ),
        "retained_foreground_fraction": retained_foreground_fraction(
            retained_count=ora_result.selected_size,
            selected_count=len(request.selected_identifiers),
        ),
        "set_identifiers_outside_background_count": (
            set_identifier_summary["set_identifiers_missing_from_background_count"]
        ),
        "set_identifiers_outside_background": (
            set_identifier_summary["set_identifiers_missing_from_background"]
        ),
        "set_identifiers_outside_background_truncated": (
            set_identifier_summary["set_identifiers_missing_from_background_truncated"]
        ),
        "set_identifiers_outside_background_preview_limit": (
            set_identifier_summary[
                "set_identifiers_missing_from_background_preview_limit"
            ]
        ),
    }


def _identifier_set_provenance_payload(
    provenance: EnrichmentIdentifierSetProvenance | None,
) -> dict[str, object] | None:
    if provenance is None:
        return None
    return provenance.to_payload()


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


__all__ = [
    "ENRICHMENT_LIMITATIONS",
    "ENRICHMENT_OFFLINE_NO_ONLINE_RESOURCE_POLICY",
    "EnrichmentRunProvenanceAssembler",
]
