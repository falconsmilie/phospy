"""Enrichment ORA execution coordination."""

from __future__ import annotations

from phospy.science.enrichment.models import ENRICHMENT_METHOD_OVER_REPRESENTATION
from phospy.science.enrichment.ora import OraEngine, OraResult
from phospy.workflows.enrichment.models import (
    InterpretedEnrichmentWorkflowRequest,
    OraEngineContract,
)
from phospy.workflows.enrichment.set_filtering import SetSizeFilterResult


class EnrichmentOraRunner:
    """Run the pure ORA kernel for prepared enrichment inputs."""

    def __init__(
        self,
        *,
        ora_engine: OraEngineContract | None = None,
    ) -> None:
        self._ora_engine = ora_engine or OraEngine()

    def run(
        self,
        *,
        request: InterpretedEnrichmentWorkflowRequest,
        set_size_filter_result: SetSizeFilterResult,
    ) -> OraResult:
        if set_size_filter_result.tested_set_collection is None:
            return _empty_ora_result(request)
        return self._ora_engine.run(
            selected_identifiers=request.selected_identifiers,
            background_universe=request.background_universe,
            enrichment_sets=set_size_filter_result.tested_set_collection,
            config=request.method_config,
        )


def _empty_ora_result(
    request: InterpretedEnrichmentWorkflowRequest,
) -> OraResult:
    usable_selected_identifiers, missing_selected_identifiers = (
        _foreground_background_intersection(
            selected_identifiers=request.selected_identifiers,
            background_universe=request.background_universe,
        )
    )
    return OraResult(
        method=ENRICHMENT_METHOD_OVER_REPRESENTATION,
        config=request.method_config,
        background_size=len(request.background_universe),
        selected_size=len(usable_selected_identifiers),
        selected_identifiers=usable_selected_identifiers,
        dropped_selected_identifiers=tuple(sorted(missing_selected_identifiers)),
        records=(),
    )


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


__all__ = ["EnrichmentOraRunner"]
