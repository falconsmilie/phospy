"""Enrichment execution set filtering."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from phospy.science.enrichment.models import EnrichmentSet, EnrichmentSetCollection
from phospy.workflows.enrichment.models import InterpretedEnrichmentWorkflowRequest

SetSizeFilterDropReason = Literal["below_min_set_size", "above_max_set_size"]
SET_SIZE_FILTER_REASON_BELOW_MIN: SetSizeFilterDropReason = "below_min_set_size"
SET_SIZE_FILTER_REASON_ABOVE_MAX: SetSizeFilterDropReason = "above_max_set_size"


@dataclass(frozen=True, slots=True)
class DroppedEnrichmentSet:
    """One set removed by enrichment min/max set-size filtering."""

    set_id: str
    name: str
    reason: SetSizeFilterDropReason
    raw_set_size: int
    background_overlap_count: int
    identifiers_outside_background_count: int


@dataclass(frozen=True, slots=True)
class SetSizeFilterResult:
    """Execution set collection after optional set-size filtering."""

    filtering_configured: bool
    tested_set_collection: EnrichmentSetCollection | None
    input_set_count: int
    tested_set_count: int
    dropped_sets: tuple[DroppedEnrichmentSet, ...]
    identifiers_outside_background_count: int


class EnrichmentSetExecutionPreparer:
    """Prepare enrichment sets for ORA execution."""

    def run(
        self,
        request: InterpretedEnrichmentWorkflowRequest,
    ) -> SetSizeFilterResult:
        if not set_size_filters_configured(request):
            return SetSizeFilterResult(
                filtering_configured=False,
                tested_set_collection=request.set_collection,
                input_set_count=len(request.set_collection.enrichment_sets),
                tested_set_count=len(request.set_collection.enrichment_sets),
                dropped_sets=(),
                identifiers_outside_background_count=0,
            )

        background = frozenset(request.background_universe)
        tested_sets: list[EnrichmentSet] = []
        dropped_sets: list[DroppedEnrichmentSet] = []
        identifiers_outside_background_count = 0
        for enrichment_set in request.set_collection.enrichment_sets:
            raw_set_size = len(enrichment_set.identifiers)
            background_overlap_count = sum(
                1
                for identifier in enrichment_set.identifiers
                if identifier in background
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
                DroppedEnrichmentSet(
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
        return SetSizeFilterResult(
            filtering_configured=True,
            tested_set_collection=tested_set_collection,
            input_set_count=len(request.set_collection.enrichment_sets),
            tested_set_count=len(tested_sets),
            dropped_sets=tuple(dropped_sets),
            identifiers_outside_background_count=identifiers_outside_background_count,
        )


def set_size_filters_configured(
    request: InterpretedEnrichmentWorkflowRequest,
) -> bool:
    return (
        request.config.min_set_size is not None
        or request.config.max_set_size is not None
    )


def dropped_set_reason_counts(
    dropped_sets: tuple[DroppedEnrichmentSet, ...],
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


__all__ = [
    "DroppedEnrichmentSet",
    "EnrichmentSetExecutionPreparer",
    "SET_SIZE_FILTER_REASON_ABOVE_MAX",
    "SET_SIZE_FILTER_REASON_BELOW_MIN",
    "SetSizeFilterDropReason",
    "SetSizeFilterResult",
    "dropped_set_reason_counts",
    "set_size_filters_configured",
]
