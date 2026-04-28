"""Internal protocol for signalome clustering backends."""

from __future__ import annotations

from typing import Protocol

from phospy.signalomes.clustering.models import (
    SignalomeClusteringBackendRequest,
    SignalomeClusteringBackendResult,
)


class SignalomeClusteringBackend(Protocol):
    """Contract implemented by numerical clustering backends."""

    @property
    def name(self) -> str: ...

    @property
    def version(self) -> str: ...

    def run(
        self,
        request: SignalomeClusteringBackendRequest,
    ) -> SignalomeClusteringBackendResult: ...


__all__ = ["SignalomeClusteringBackend"]
