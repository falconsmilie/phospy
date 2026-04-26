"""Reference domain models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass
from enum import Enum

import pandas as pd

from phospy._frame_ownership import own_dataframe
from phospy.errors.validation import ReferenceValidationError
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import ReferenceProvenance
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference


class Organism(str, Enum):
    """Public organism identifiers used in dataset/reference contracts.

    Enum membership defines contract syntax. Bundled runtime scientific support
    may be narrower in a given release.
    """

    HUMAN = "human"
    MOUSE = "mouse"
    RAT = "rat"


class ReferencePreset(str, Enum):
    """Built-in organism presets for bundled-reference resolution.

    Enum values define public organism lanes accepted by request contracts.
    Bundled runtime references may cover only a subset in a given release.
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
    provenance: ReferenceProvenance | None = None
    _assume_owned: InitVar[bool] = False

    def __post_init__(self, _assume_owned: bool) -> None:
        if not isinstance(self.organism, Organism):
            raise ReferenceValidationError(
                "references.organism must be an Organism enum value"
            )
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
        kinase_substrate_reference = KinaseSubstrateReference(
            frame=kinase_substrate_map,
            _assume_owned=True,
        )
        site_sequence_reference = SiteSequenceReference(
            frame=site_sequences,
            _assume_owned=True,
        )
        substrate_sites = set(
            kinase_substrate_reference.frame.loc[:, "substrate_site"].tolist()
        )
        known_sites = set(site_sequence_reference.frame.index.tolist())
        missing_sequences = sorted(substrate_sites.difference(known_sites))
        if missing_sequences:
            missing_sample = ", ".join(missing_sequences[:10])
            raise ReferenceValidationError(
                "references.site_sequences is missing sequence entries for "
                f"substrate sites in references.kinase_substrate_map: {missing_sample}"
            )
        provenance = self.provenance
        if provenance is None:
            provenance = ReferenceProvenance(
                source_type="explicit",
                organism=self.organism.value,
                bundle_id=None,
                table_fingerprints=(
                    fingerprint_table(
                        kinase_substrate_reference.frame,
                        name="references.kinase_substrate_map",
                    ),
                    fingerprint_table(
                        site_sequence_reference.frame,
                        name="references.site_sequences",
                    ),
                ),
            )
        elif not isinstance(provenance, ReferenceProvenance):
            raise ReferenceValidationError(
                "references.provenance must be ReferenceProvenance or None"
            )
        object.__setattr__(
            self,
            "kinase_substrate_map",
            kinase_substrate_reference.frame,
        )
        object.__setattr__(
            self,
            "site_sequences",
            site_sequence_reference.frame,
        )
        object.__setattr__(self, "provenance", provenance)

    @classmethod
    def _from_owned(
        cls,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
        provenance: ReferenceProvenance | None = None,
    ) -> ReferenceBundle:
        return cls(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            _assume_owned=True,
        )
