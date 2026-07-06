"""Internal executor for enrichment workflow requests."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal, cast

import pandas as pd

from phospy.contracts.results import EnrichmentWorkflowResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.provenance import (
    RunProvenance,
    collect_environment_provenance,
    fingerprint_table,
)
from phospy.science.enrichment.models import (
    ENRICHMENT_METHOD_OVER_REPRESENTATION,
    MULTIPLE_TESTING_CORRECTION_NONE,
    EnrichmentResultRecord,
    EnrichmentSet,
    EnrichmentSetCollection,
)
from phospy.science.enrichment.ora import (
    OraConfig,
    OraEngine,
    OraResult,
    OraResultRecord,
)
from phospy.science.statistics.multiple_testing import (
    MultipleTestingCorrection,
)
from phospy.science.statistics.multiple_testing import (
    run as run_multiple_testing_correction,
)
from phospy.workflows.enrichment.caveats import build_enrichment_result_caveats
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
SetSizeFilterDropReason = Literal["below_min_set_size", "above_max_set_size"]
SET_SIZE_FILTER_REASON_BELOW_MIN: SetSizeFilterDropReason = "below_min_set_size"
SET_SIZE_FILTER_REASON_ABOVE_MAX: SetSizeFilterDropReason = "above_max_set_size"
MAX_UNMATCHED_SET_IDENTIFIER_DIAGNOSTIC_COUNT = 50


@dataclass(frozen=True, slots=True)
class _DroppedEnrichmentSet:
    set_id: str
    name: str
    reason: SetSizeFilterDropReason
    raw_set_size: int
    background_overlap_count: int
    identifiers_outside_background_count: int


@dataclass(frozen=True, slots=True)
class _SetSizeFilterResult:
    filtering_configured: bool
    tested_set_collection: EnrichmentSetCollection | None
    input_set_count: int
    tested_set_count: int
    dropped_sets: tuple[_DroppedEnrichmentSet, ...]
    identifiers_outside_background_count: int


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
        set_size_filter_result = _prepare_set_collection_for_execution(request)
        if set_size_filter_result.tested_set_collection is None:
            ora_result = _empty_ora_result(request=request, config=raw_config)
        else:
            ora_result = self._ora_engine.run(
                selected_identifiers=request.selected_identifiers,
                background_universe=request.background_universe,
                enrichment_sets=set_size_filter_result.tested_set_collection,
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
        preliminary = EnrichmentWorkflowResult(
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
            caveats=caveats,
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


def _prepare_set_collection_for_execution(
    request: InterpretedEnrichmentWorkflowRequest,
) -> _SetSizeFilterResult:
    if not _set_size_filters_configured(request):
        return _SetSizeFilterResult(
            filtering_configured=False,
            tested_set_collection=request.set_collection,
            input_set_count=len(request.set_collection.enrichment_sets),
            tested_set_count=len(request.set_collection.enrichment_sets),
            dropped_sets=(),
            identifiers_outside_background_count=0,
        )

    background = frozenset(request.background_universe)
    tested_sets: list[EnrichmentSet] = []
    dropped_sets: list[_DroppedEnrichmentSet] = []
    identifiers_outside_background_count = 0
    for enrichment_set in request.set_collection.enrichment_sets:
        raw_set_size = len(enrichment_set.identifiers)
        background_overlap_count = sum(
            1 for identifier in enrichment_set.identifiers if identifier in background
        )
        set_outside_background_count = raw_set_size - background_overlap_count
        identifiers_outside_background_count += set_outside_background_count
        drop_reason = _set_size_filter_drop_reason(
            background_overlap_count=background_overlap_count,
            min_set_size=request.config.min_set_size,
            max_set_size=request.config.max_set_size,
        )
        if drop_reason is None:
            tested_sets.append(enrichment_set)
            continue
        dropped_sets.append(
            _DroppedEnrichmentSet(
                set_id=enrichment_set.set_id,
                name=enrichment_set.name,
                reason=drop_reason,
                raw_set_size=raw_set_size,
                background_overlap_count=background_overlap_count,
                identifiers_outside_background_count=set_outside_background_count,
            )
        )

    tested_set_collection = (
        None
        if not tested_sets
        else EnrichmentSetCollection(
            sets=tuple(tested_sets),
            identifier_kind=request.set_collection.identifier_kind,
            collection_kind=request.set_collection.collection_kind,
            source_name=request.set_collection.source_name,
            source_version=request.set_collection.source_version,
        )
    )
    return _SetSizeFilterResult(
        filtering_configured=True,
        tested_set_collection=tested_set_collection,
        input_set_count=len(request.set_collection.enrichment_sets),
        tested_set_count=len(tested_sets),
        dropped_sets=tuple(dropped_sets),
        identifiers_outside_background_count=identifiers_outside_background_count,
    )


def _set_size_filters_configured(
    request: InterpretedEnrichmentWorkflowRequest,
) -> bool:
    return (
        request.config.min_set_size is not None
        or request.config.max_set_size is not None
    )


def _set_size_filter_drop_reason(
    *,
    background_overlap_count: int,
    min_set_size: int | None,
    max_set_size: int | None,
) -> SetSizeFilterDropReason | None:
    if min_set_size is not None and background_overlap_count < min_set_size:
        return SET_SIZE_FILTER_REASON_BELOW_MIN
    if max_set_size is not None and background_overlap_count > max_set_size:
        return SET_SIZE_FILTER_REASON_ABOVE_MAX
    return None


def _empty_ora_result(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    config: OraConfig,
) -> OraResult:
    usable_selected_identifiers, missing_selected_identifiers = (
        _foreground_background_intersection(
            selected_identifiers=request.selected_identifiers,
            background_universe=request.background_universe,
        )
    )
    return OraResult(
        method=ENRICHMENT_METHOD_OVER_REPRESENTATION,
        config=config,
        background_size=len(request.background_universe),
        selected_size=len(usable_selected_identifiers),
        selected_identifiers=usable_selected_identifiers,
        dropped_selected_identifiers=tuple(sorted(missing_selected_identifiers)),
        records=(),
    )


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
    _, missing_foreground_identifiers = _foreground_background_intersection(
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
        }
    )
    return summary


def _execution_set_collection_summary(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    ora_result: OraResult,
    set_size_filter_result: _SetSizeFilterResult,
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
                "dropped_set_reason_counts": _dropped_set_reason_counts(
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
    set_size_filter_result: _SetSizeFilterResult,
) -> dict[str, object]:
    diagnostics = dict(request.diagnostics)
    no_hit_set_ids = tuple(
        record.set_id for record in ora_result.records if record.overlap_size == 0
    )
    diagnostics.update(
        {
            "foreground_background": _foreground_background_diagnostics(
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


def _foreground_background_intersection(
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


def _foreground_background_diagnostics(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    set_size_filter_result: _SetSizeFilterResult,
) -> dict[str, object]:
    usable_foreground, missing_foreground = _foreground_background_intersection(
        selected_identifiers=request.selected_identifiers,
        background_universe=request.background_universe,
    )
    unmatched_set_identifier_summary = _unmatched_set_identifier_summary(
        request=request
    )
    return {
        "identifier_kind": request.identifier_semantics.identifier_kind,
        "foreground_size_before_intersection": len(request.selected_identifiers),
        "background_size": len(request.background_universe),
        "usable_foreground_size_after_background_intersection": len(usable_foreground),
        "foreground_identifiers_missing_from_background_count": len(missing_foreground),
        "foreground_identifiers_missing_from_background": missing_foreground,
        "tested_set_count": set_size_filter_result.tested_set_count,
        "dropped_set_count": len(set_size_filter_result.dropped_sets),
        **unmatched_set_identifier_summary,
    }


def _unmatched_set_identifier_summary(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
) -> dict[str, object]:
    background = frozenset(request.background_universe)
    missing_identifiers: list[str] = []
    seen: set[str] = set()
    for enrichment_set in request.set_collection.enrichment_sets:
        for identifier in enrichment_set.identifiers:
            if identifier in background or identifier in seen:
                continue
            seen.add(identifier)
            missing_identifiers.append(identifier)

    preview = tuple(missing_identifiers[:MAX_UNMATCHED_SET_IDENTIFIER_DIAGNOSTIC_COUNT])
    return {
        "set_identifiers_missing_from_background_count": len(missing_identifiers),
        "set_identifiers_missing_from_background": preview,
        "set_identifiers_missing_from_background_truncated": (
            len(missing_identifiers) > MAX_UNMATCHED_SET_IDENTIFIER_DIAGNOSTIC_COUNT
        ),
        "set_identifiers_missing_from_background_preview_limit": (
            MAX_UNMATCHED_SET_IDENTIFIER_DIAGNOSTIC_COUNT
        ),
    }


def _set_size_filter_diagnostics(
    *,
    request: InterpretedEnrichmentWorkflowRequest,
    set_size_filter_result: _SetSizeFilterResult,
) -> dict[str, object]:
    return {
        "applied_after_background_intersection": True,
        "min_set_size": request.config.min_set_size,
        "max_set_size": request.config.max_set_size,
        "input_set_count": set_size_filter_result.input_set_count,
        "tested_set_count": set_size_filter_result.tested_set_count,
        "dropped_set_count": len(set_size_filter_result.dropped_sets),
        "dropped_set_reason_counts": _dropped_set_reason_counts(
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


def _dropped_set_reason_counts(
    dropped_sets: tuple[_DroppedEnrichmentSet, ...],
) -> dict[str, int]:
    return {
        SET_SIZE_FILTER_REASON_BELOW_MIN: sum(
            1
            for dropped_set in dropped_sets
            if dropped_set.reason == SET_SIZE_FILTER_REASON_BELOW_MIN
        ),
        SET_SIZE_FILTER_REASON_ABOVE_MAX: sum(
            1
            for dropped_set in dropped_sets
            if dropped_set.reason == SET_SIZE_FILTER_REASON_ABOVE_MAX
        ),
    }


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
    if _set_size_filters_configured(request):
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
