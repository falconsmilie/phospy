"""Enrichment workflow public result assembly."""

from __future__ import annotations

from phospy.contracts.results import EnrichmentWorkflowResult
from phospy.provenance import RunProvenance
from phospy.science.enrichment.models import (
    MULTIPLE_TESTING_CORRECTION_NONE,
    EnrichmentResultRecord,
)
from phospy.science.enrichment.ora import OraResult, OraResultRecord
from phospy.science.statistics.multiple_testing import MultipleTestingCorrection
from phospy.workflows.enrichment.caveats import build_enrichment_result_caveats
from phospy.workflows.enrichment.models import InterpretedEnrichmentWorkflowRequest
from phospy.workflows.enrichment.set_filtering import (
    SetSizeFilterResult,
    dropped_set_reason_counts,
)

MAX_UNMATCHED_SET_IDENTIFIER_DIAGNOSTIC_COUNT = 50


class EnrichmentResultAssembler:
    """Assemble public enrichment workflow results from ORA outputs."""

    def run(
        self,
        *,
        request: InterpretedEnrichmentWorkflowRequest,
        ora_result: OraResult,
        set_size_filter_result: SetSizeFilterResult,
        provenance: RunProvenance | None = None,
    ) -> EnrichmentWorkflowResult:
        records = _build_result_records(ora_records=ora_result.records)
        unmatched_identifiers = _unmatched_selected_identifiers(ora_result)
        background_summary = _execution_background_summary(
            request=request,
            ora_result=ora_result,
        )
        set_collection_summary = _execution_set_collection_summary(
            request=request,
            ora_result=ora_result,
            set_size_filter_result=set_size_filter_result,
        )
        diagnostics = _execution_diagnostics(
            request=request,
            ora_result=ora_result,
            unmatched_identifiers=unmatched_identifiers,
            set_size_filter_result=set_size_filter_result,
        )
        caveats = build_enrichment_result_caveats(
            request=request,
            background_summary=background_summary,
            set_collection_summary=set_collection_summary,
        )
        return EnrichmentWorkflowResult(
            identifier_kind=request.identifier_semantics.identifier_kind,
            set_collection=request.set_collection,
            config=request.config,
            records=records,
            unmatched_identifiers=unmatched_identifiers,
            caveats=caveats,
            diagnostics=diagnostics,
            method_metadata=request.method_metadata,
            background_summary=background_summary,
            set_collection_summary=set_collection_summary,
            selected_identifier_provenance=request.selected_identifier_provenance,
            background_identifier_provenance=request.background_identifier_provenance,
            provenance=provenance,
        )


def _build_result_records(
    *,
    ora_records: tuple[OraResultRecord, ...],
) -> tuple[EnrichmentResultRecord, ...]:
    records: list[EnrichmentResultRecord] = []
    for ora_record in ora_records:
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
                adjusted_p_value=ora_record.adjusted_p_value,
                correction_method=ora_record.correction_method,
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
    _, missing_foreground_identifiers = foreground_background_intersection(
        selected_identifiers=request.selected_identifiers,
        background_universe=request.background_universe,
    )
    summary.update(
        {
            "universe_size": ora_result.background_size,
            "selected_in_background_count": ora_result.selected_size,
            "dropped_selected_count": len(ora_result.dropped_selected_identifiers),
            "dropped_selected_identifiers": ora_result.dropped_selected_identifiers,
            "foreground_size_before_intersection": len(request.selected_identifiers),
            "usable_foreground_size_after_background_intersection": (
                ora_result.selected_size
            ),
            "foreground_identifiers_missing_from_background_count": len(
                missing_foreground_identifiers
            ),
            "foreground_identifiers_missing_from_background": (
                missing_foreground_identifiers
            ),
            "retained_foreground_fraction": retained_foreground_fraction(
                retained_count=ora_result.selected_size,
                selected_count=len(request.selected_identifiers),
            ),
            "selected_outside_background_policy": (
                request.method_config.selected_outside_background_policy
            ),
            "set_member_outside_background_policy": (
                request.config.set_member_outside_background_policy
            ),
            "minimum_retained_foreground_fraction": (
                request.config.minimum_retained_foreground_fraction
            ),
        }
    )
    return summary


def _execution_set_collection_summary(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
    set_size_filter_result: SetSizeFilterResult,
) -> dict[str, object]:
    empty_after_background = tuple(
        record.set_id for record in ora_result.records if record.set_size == 0
    )
    summary = dict(request.set_collection_summary)
    identifiers_outside_background_count = sum(
        record.set_identifiers_outside_background_count for record in ora_result.records
    )
    if set_size_filter_result.filtering_configured:
        identifiers_outside_background_count = (
            set_size_filter_result.identifiers_outside_background_count
        )
    summary.update(
        {
            "identifiers_outside_background_count": identifiers_outside_background_count,
            "empty_after_background_count": len(empty_after_background),
            "empty_after_background_set_ids": empty_after_background,
        }
    )
    if set_size_filter_result.filtering_configured:
        summary.update(
            {
                "tested_set_count": set_size_filter_result.tested_set_count,
                "dropped_set_count": len(set_size_filter_result.dropped_sets),
                "dropped_set_ids": tuple(
                    dropped_set.set_id
                    for dropped_set in set_size_filter_result.dropped_sets
                ),
                "dropped_set_reason_counts": dropped_set_reason_counts(
                    set_size_filter_result.dropped_sets
                ),
            }
        )
    return summary


def _execution_diagnostics(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
    unmatched_identifiers: tuple[str, ...],
    set_size_filter_result: SetSizeFilterResult,
) -> dict[str, object]:
    diagnostics = dict(request.diagnostics)
    no_hit_set_ids = tuple(
        record.set_id for record in ora_result.records if record.overlap_size == 0
    )
    diagnostics.update(
        {
            "foreground_background": foreground_background_diagnostics(
                request=request,
                set_size_filter_result=set_size_filter_result,
            ),
            "ora": {
                "method": ora_result.method,
                "record_count": len(ora_result.records),
                "selected_size": ora_result.selected_size,
                "background_size": ora_result.background_size,
                "no_hit_set_count": len(no_hit_set_ids),
                "no_hit_set_ids": no_hit_set_ids,
            },
            "multiple_testing_correction": (
                _multiple_testing_correction_diagnostics(
                    method=request.method_config.multiple_testing_correction,
                    tested_record_count=len(ora_result.records),
                )
            ),
            "unmatched_selected_identifier_count": len(unmatched_identifiers),
        }
    )
    if set_size_filter_result.filtering_configured:
        diagnostics["set_size_filter"] = _set_size_filter_diagnostics(
            request=request,
            set_size_filter_result=set_size_filter_result,
        )
    return diagnostics


def _multiple_testing_correction_diagnostics(
    *,
    method: MultipleTestingCorrection,
    tested_record_count: int,
) -> dict[str, object]:
    return {
        "method": method,
        "applied": (
            method != MULTIPLE_TESTING_CORRECTION_NONE and tested_record_count > 0
        ),
        "tested_record_count": tested_record_count,
    }


def foreground_background_intersection(
    *,
    selected_identifiers: tuple[str, ...],
    background_universe: tuple[str, ...],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    background = frozenset(background_universe)
    usable_foreground = tuple(
        identifier for identifier in selected_identifiers if identifier in background
    )
    missing_foreground = tuple(
        identifier
        for identifier in selected_identifiers
        if identifier not in background
    )
    return usable_foreground, missing_foreground


def retained_foreground_fraction(
    *,
    retained_count: int,
    selected_count: int,
) -> float | None:
    if selected_count <= 0:
        return None
    return float(retained_count) / float(selected_count)


def foreground_background_diagnostics(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    set_size_filter_result: SetSizeFilterResult,
) -> dict[str, object]:
    usable_foreground, missing_foreground = foreground_background_intersection(
        selected_identifiers=request.selected_identifiers,
        background_universe=request.background_universe,
    )
    unmatched_set_identifier_summary = unmatched_set_identifier_summary_for_request(
        request=request
    )
    return {
        "identifier_kind": request.identifier_semantics.identifier_kind,
        "foreground_size_before_intersection": len(request.selected_identifiers),
        "background_size": len(request.background_universe),
        "usable_foreground_size_after_background_intersection": len(usable_foreground),
        "retained_foreground_fraction": retained_foreground_fraction(
            retained_count=len(usable_foreground),
            selected_count=len(request.selected_identifiers),
        ),
        "foreground_identifiers_missing_from_background_count": len(missing_foreground),
        "foreground_identifiers_missing_from_background": missing_foreground,
        "selected_outside_background_policy": (
            request.method_config.selected_outside_background_policy
        ),
        "set_member_outside_background_policy": (
            request.config.set_member_outside_background_policy
        ),
        "minimum_retained_foreground_fraction": (
            request.config.minimum_retained_foreground_fraction
        ),
        "tested_set_count": set_size_filter_result.tested_set_count,
        "dropped_set_count": len(set_size_filter_result.dropped_sets),
        **unmatched_set_identifier_summary,
    }


def unmatched_set_identifier_summary_for_request(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
) -> dict[str, object]:
    background = frozenset(request.background_universe)
    missing_identifiers: set[str] = set()
    for enrichment_set in request.set_collection.enrichment_sets:
        for identifier in enrichment_set.identifiers:
            if identifier in background:
                continue
            missing_identifiers.add(identifier)

    ordered_missing_identifiers = tuple(sorted(missing_identifiers))
    preview = ordered_missing_identifiers[
        :MAX_UNMATCHED_SET_IDENTIFIER_DIAGNOSTIC_COUNT
    ]
    return {
        "set_identifiers_missing_from_background_count": len(
            ordered_missing_identifiers
        ),
        "set_identifiers_missing_from_background": preview,
        "set_identifiers_missing_from_background_truncated": (
            len(ordered_missing_identifiers)
            > MAX_UNMATCHED_SET_IDENTIFIER_DIAGNOSTIC_COUNT
        ),
        "set_identifiers_missing_from_background_preview_limit": (
            MAX_UNMATCHED_SET_IDENTIFIER_DIAGNOSTIC_COUNT
        ),
    }


def _set_size_filter_diagnostics(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    set_size_filter_result: SetSizeFilterResult,
) -> dict[str, object]:
    return {
        "applied_after_background_intersection": True,
        "min_set_size": request.config.min_set_size,
        "max_set_size": request.config.max_set_size,
        "input_set_count": set_size_filter_result.input_set_count,
        "tested_set_count": set_size_filter_result.tested_set_count,
        "dropped_set_count": len(set_size_filter_result.dropped_sets),
        "dropped_set_reason_counts": dropped_set_reason_counts(
            set_size_filter_result.dropped_sets
        ),
        "dropped_sets": tuple(
            {
                "set_id": dropped_set.set_id,
                "term_name": dropped_set.name,
                "reason": dropped_set.reason,
                "raw_set_size": dropped_set.raw_set_size,
                "background_overlap_count": dropped_set.background_overlap_count,
                "identifiers_outside_background_count": (
                    dropped_set.identifiers_outside_background_count
                ),
            }
            for dropped_set in set_size_filter_result.dropped_sets
        ),
    }


__all__ = [
    "EnrichmentResultAssembler",
    "foreground_background_diagnostics",
    "foreground_background_intersection",
    "retained_foreground_fraction",
    "unmatched_set_identifier_summary_for_request",
]
