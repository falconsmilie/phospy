"""Internal protocol for signalome clustering backends."""

from __future__ import annotations

from typing import Protocol

from phospy.signalomes.clustering.models import (
    SignalomeClusteringBackendRequest,
    SignalomeClusteringBackendResult,
)


class SignalomeClusteringBackend(Protocol):
    """Contract implemented by numerical clustering backends."""

    name: str
    version: str

    def run(
        self,
        request: SignalomeClusteringBackendRequest,
    ) -> SignalomeClusteringBackendResult:
        """Cluster sites and derive module assignments."""


__all__ = ["SignalomeClusteringBackend"]
