"""Immutable Enrichment plans and collaborator contracts."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal, Protocol

from phospy.enrichment.contracts.collections import EnrichmentSetCollection
from phospy.enrichment.contracts.config import EnrichmentConfig
from phospy.enrichment.contracts.identifier_sets import (
    EnrichmentIdentifierSetProvenance,
)
from phospy.enrichment.contracts.requests import EnrichmentWorkflowRequest
from phospy.enrichment.contracts.results import EnrichmentWorkflowResult
from phospy.enrichment.contracts.vocabulary import (
    EnrichmentCollectionKind,
    EnrichmentIdentifierKind,
)
from phospy.enrichment.ora import OraConfig, OraResult

EnrichmentAnalysisLevel = Literal["gene", "ptm"]
EnrichmentSelectedIdentifierSource = Literal["selected_identifiers", "input_table"]


@dataclass(frozen=True, slots=True)
class ValidatedEnrichmentWorkflowRequest:
    """Validated enrichment request passed to interpretation."""

    request: EnrichmentWorkflowRequest
    identifier_column: str
    identifier_kind: EnrichmentIdentifierKind
    set_collection: EnrichmentSetCollection
    background_universe: tuple[str, ...]
    selected_identifiers: tuple[str, ...]
    config: EnrichmentConfig
    selected_identifier_source: EnrichmentSelectedIdentifierSource
    selected_identifier_input_count: int
    background_identifier_input_count: int
    selected_identifier_provenance: EnrichmentIdentifierSetProvenance | None
    background_identifier_provenance: EnrichmentIdentifierSetProvenance | None


@dataclass(frozen=True, slots=True)
class EnrichmentIdentifierSemantics:
    """Execution-time identifier meaning resolved by interpretation."""

    identifier_column: str
    identifier_kind: EnrichmentIdentifierKind
    collection_kind: EnrichmentCollectionKind
    analysis_level: EnrichmentAnalysisLevel


@dataclass(frozen=True, slots=True)
class InterpretedEnrichmentWorkflowRequest:
    """Execution-ready enrichment request produced by interpretation."""

    selected_identifiers: tuple[str, ...]
    background_universe: tuple[str, ...]
    set_collection: EnrichmentSetCollection
    method_config: OraConfig
    identifier_semantics: EnrichmentIdentifierSemantics
    config: EnrichmentConfig
    selected_identifier_source: EnrichmentSelectedIdentifierSource
    method_metadata: dict[str, object]
    background_summary: dict[str, object]
    set_collection_summary: dict[str, object]
    diagnostics: dict[str, object]
    selected_identifier_input_count: int
    background_identifier_input_count: int
    selected_identifier_provenance: EnrichmentIdentifierSetProvenance | None
    background_identifier_provenance: EnrichmentIdentifierSetProvenance | None


class EnrichmentWorkflowValidatorContract(Protocol):
    """Internal contract for enrichment workflow validation."""

    def run(self, request: object) -> ValidatedEnrichmentWorkflowRequest: ...


class EnrichmentWorkflowInterpreterContract(Protocol):
    """Internal contract for enrichment workflow interpretation."""

    def run(
        self, request: ValidatedEnrichmentWorkflowRequest
    ) -> InterpretedEnrichmentWorkflowRequest: ...


class OraEngineContract(Protocol):
    """Internal contract for ORA execution."""

    def run(
        self,
        *,
        selected_identifiers: Sequence[str],
        background_universe: Sequence[str],
        enrichment_sets: EnrichmentSetCollection,
        config: OraConfig | None = None,
    ) -> OraResult: ...


class EnrichmentWorkflowExecutorContract(Protocol):
    """Internal contract for enrichment workflow execution."""

    def run(
        self, request: InterpretedEnrichmentWorkflowRequest
    ) -> EnrichmentWorkflowResult: ...


__all__ = [
    "EnrichmentAnalysisLevel",
    "EnrichmentIdentifierSemantics",
    "EnrichmentSelectedIdentifierSource",
    "EnrichmentWorkflowExecutorContract",
    "EnrichmentWorkflowInterpreterContract",
    "EnrichmentWorkflowValidatorContract",
    "InterpretedEnrichmentWorkflowRequest",
    "OraEngineContract",
    "ValidatedEnrichmentWorkflowRequest",
]
