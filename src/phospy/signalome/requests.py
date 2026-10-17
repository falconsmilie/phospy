"""Public Signalome workflow request."""

from __future__ import annotations

from dataclasses import dataclass, field

from phospy.kinase.results import KinaseWorkflowResult
from phospy.signalome.config import SignalomeConfig


@dataclass(frozen=True, slots=True, eq=False)
class SignalomeWorkflowRequest:
    """Request for the public signalome workflow.

    Construction stores the payload only. The upstream kinase result, signalome
    config, score/prediction matrix alignment, site identity, and protein
    grouping metadata are validated by ``SignalomeWorkflow.run(...)`` before
    interpretation or execution.

    Signalome execution requires explicit protein grouping metadata per
    interpreted site via ``dataset.site_metadata.protein_group_id``. The legacy
    ``dataset.site_metadata.protein_id`` column is accepted only as a migration
    alias.

    Python equality and hashing are identity-based. The workflow owns semantic
    request interpretation; this payload does not define content equality.
    """

    __hash__ = object.__hash__

    kinase_result: KinaseWorkflowResult
    config: SignalomeConfig = field(default_factory=SignalomeConfig)


__all__ = ["SignalomeWorkflowRequest"]
