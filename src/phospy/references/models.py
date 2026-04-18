"""Reference domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass
from enum import Enum

import pandas as pd

from phospy._frame_ownership import own_dataframe
from phospy.errors.validation import ReferenceValidationError
from phospy.site_ids import canonicalize_site_index, canonicalize_site_series


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
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        kinase_substrate_map = own_dataframe(
            self.kinase_substrate_map,
            field_name="references.kinase_substrate_map",
            error_type=ReferenceValidationError,
            assume_owned=_assume_owned,
        )
        site_sequences = own_dataframe(
            self.site_sequences,
            field_name="references.site_sequences",
            error_type=ReferenceValidationError,
            assume_owned=_assume_owned,
        )
        if "kinase" in kinase_substrate_map.columns:
            kinase_substrate_map = kinase_substrate_map.assign(
                kinase=kinase_substrate_map.loc[:, "kinase"].astype(str).str.strip()
            )
        if "substrate_site" in kinase_substrate_map.columns:
            kinase_substrate_map = kinase_substrate_map.assign(
                substrate_site=canonicalize_site_series(
                    kinase_substrate_map.loc[:, "substrate_site"],
                    field_name="references.kinase_substrate_map.substrate_site",
                    error_type=ReferenceValidationError,
                )
            )
        if {"kinase", "substrate_site"}.issubset(kinase_substrate_map.columns):
            kinase_substrate_map = kinase_substrate_map.drop_duplicates(
                subset=["kinase", "substrate_site"],
                ignore_index=True,
            )
        site_sequences.index = canonicalize_site_index(
            site_sequences.index,
            field_name="references.site_sequences.index",
            error_type=ReferenceValidationError,
            index_name="site_id",
        )
        if "site_sequence" in site_sequences.columns:
            site_sequences = site_sequences.assign(
                site_sequence=site_sequences.loc[:, "site_sequence"]
                .astype(str)
                .str.strip()
            )

        from phospy.validation.references.bundle import ReferenceBundleValidator

        ReferenceBundleValidator().run(
            organism=self.organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
        )
        object.__setattr__(self, "kinase_substrate_map", kinase_substrate_map)
        object.__setattr__(self, "site_sequences", site_sequences)

    @classmethod
    def _from_owned(
        cls,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
    ) -> ReferenceBundle:
        return cls(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            _assume_owned=True,
        )
