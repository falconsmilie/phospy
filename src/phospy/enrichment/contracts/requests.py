"""Passive Enrichment workflow request contract."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path

import pandas as pd

from phospy.enrichment.contracts.collections import EnrichmentSetCollection
from phospy.enrichment.contracts.config import EnrichmentConfig
from phospy.enrichment.contracts.identifier_sets import (
    EnrichmentIdentifierSetProvenance,
)
from phospy.enrichment.contracts.vocabulary import EnrichmentIdentifierKind


@dataclass(frozen=True, slots=True, eq=False)
class EnrichmentWorkflowRequest:
    """Request for native enrichment.

    Construction stores the payload only. Identifier semantics, collection and
    background compatibility, input-source selection, input emptiness, and
    config type are validated by ``EnrichmentWorkflow.run(...)`` before
    interpretation or execution.

    Python equality and hashing are identity-based. The workflow owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

    identifier_column: str
    identifier_kind: EnrichmentIdentifierKind
    set_collection: EnrichmentSetCollection
    background_universe: Sequence[str]
    config: EnrichmentConfig = field(default_factory=EnrichmentConfig)
    input_table: pd.DataFrame | str | Path | PathLike[str] | None = None
    selected_identifiers: Sequence[str] | None = None
    selected_identifier_provenance: EnrichmentIdentifierSetProvenance | None = None
    background_identifier_provenance: EnrichmentIdentifierSetProvenance | None = None


__all__ = ["EnrichmentWorkflowRequest"]
