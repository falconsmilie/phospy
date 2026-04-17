"""Reference domain models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import pandas as pd


class Organism(str, Enum):
    """Supported organism identifiers."""

    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"


class ReferencePreset(str, Enum):
    """Built-in reference selection presets."""

    AUTO = "auto"
    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"


@dataclass(frozen=True, slots=True)
class ReferenceBundle:
    """Resolved workflow reference resources."""

    organism: Organism
    kinase_substrate_map: pd.DataFrame
    site_sequences: pd.DataFrame
