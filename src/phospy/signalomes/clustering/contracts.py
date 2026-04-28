"""Shared contracts for signalome clustering orchestration and engines."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Protocol

import numpy as np

from phospy.signalomes.clustering.models import (
    SignalomeClusteringBackendRequest,
    SignalomeClusteringBackendResult,
)


class SignalomeClusteringBackend(Protocol):
    """Contract implemented by top-level clustering backends."""

    @property
    def name(self) -> str: ...

    @property
    def version(self) -> str: ...

    def run(
        self,
        request: SignalomeClusteringBackendRequest,
    ) -> SignalomeClusteringBackendResult: ...


class ClusterTreeEngine(Protocol):
    """Low-level tree engine contract used by shared clustering orchestration."""

    @property
    def name(self) -> str: ...

    @property
    def version(self) -> str: ...

    def build_tree(self, values: np.ndarray) -> object: ...

    def labels_for_counts(
        self,
        *,
        tree: object,
        counts: Iterable[int],
    ) -> dict[int, np.ndarray]: ...


__all__ = ["ClusterTreeEngine", "SignalomeClusteringBackend"]
