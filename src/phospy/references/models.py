"""Reference domain models."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import pandas as pd


class Organism(str, Enum):
    """Supported organism identifiers for dataset/reference contracts."""

    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"


class ReferencePreset(str, Enum):
    """Built-in organism presets for bundled-reference resolution.

    Enum values define the public organism lanes. A release may bundle only a
    subset of those lanes.
    """

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

    def __post_init__(self) -> None:
        kinase_substrate_map = _copy_frame(self.kinase_substrate_map)
        site_sequences = _copy_frame(self.site_sequences)

        from phospy.validation.references.bundle import ReferenceBundleValidator

        ReferenceBundleValidator().run(
            organism=self.organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
        )
        object.__setattr__(self, "kinase_substrate_map", kinase_substrate_map)
        object.__setattr__(self, "site_sequences", site_sequences)


def _copy_frame(value: object) -> object:
    if not isinstance(value, pd.DataFrame):
        return value
    return value.copy(deep=True)
