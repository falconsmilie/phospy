"""Reference bundle validation service."""

from __future__ import annotations

from phospy.references.manifest import ReferenceManifest
from phospy.references.models import (
    ReferenceBundleValidationResult,
    run_reference_bundle_validation,
)
from phospy.references.provenance.models import ReferenceProvenance
from phospy.taxonomy.organisms import Organism


class ReferenceBundleValidator:
    """Validate the stable `ReferenceBundle` contract."""

    def run(
        self,
        *,
        organism: Organism,
        kinase_substrate_map: object,
        site_sequences: object,
        provenance: ReferenceProvenance | None = None,
        manifest: ReferenceManifest | None = None,
    ) -> ReferenceBundleValidationResult:
        return run_reference_bundle_validation(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
        )


__all__ = ["ReferenceBundleValidationResult", "ReferenceBundleValidator"]
