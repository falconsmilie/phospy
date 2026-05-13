"""Reference bundle validator."""

from __future__ import annotations

import pandas as pd

from phospy.errors.validation import ReferenceValidationError
from phospy.science.references.models import Organism
from phospy.tables.references import KinaseSubstrateReference, SiteSequenceReference


class ReferenceBundleValidator:
    """Validate the stable `ReferenceBundle` contract."""

    def run(
        self,
        *,
        organism: Organism,
        kinase_substrate_map: pd.DataFrame,
        site_sequences: pd.DataFrame,
    ) -> None:
        if not isinstance(organism, Organism):
            raise ReferenceValidationError(
                "references.organism must be an Organism enum value"
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
