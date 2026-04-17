"""Signalome domain models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True, slots=True)
class SignalomeAssignments:
    """Signalome module assignment table."""

    table: pd.DataFrame


@dataclass(frozen=True, slots=True)
class SignalomeModules:
    """Signalome module table."""

    table: pd.DataFrame


@dataclass(frozen=True, slots=True)
class KinaseNetwork:
    """Kinase network tables derived from signalome analysis."""

    edges: pd.DataFrame
    nodes: pd.DataFrame | None = None
