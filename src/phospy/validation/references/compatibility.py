"""Shared reference compatibility validation."""

from __future__ import annotations

from phospy.errors.references import (
    ReferenceCompatibilityError,
    ReferenceResolutionError,
)
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset


class ReferenceCompatibilityValidator:
    """Validate dataset/reference compatibility for preset or explicit bundle input."""

    def run(
        self,
        reference_input: object,
        *,
        dataset_organism: Organism | None,
    ) -> None:
        if not isinstance(reference_input, (ReferencePreset, ReferenceBundle)):
            raise ReferenceResolutionError(
                "reference input must be a ReferencePreset or ReferenceBundle"
            )
        if isinstance(reference_input, ReferenceBundle):
            self.run_bundle_organism(
                reference_organism=reference_input.organism,
                dataset_organism=dataset_organism,
            )
            return
        if reference_input is ReferencePreset.AUTO:
            if dataset_organism is None:
                raise ReferenceResolutionError(
                    "ReferencePreset.AUTO requires dataset.organism"
                )
            return
        target_organism = {
            ReferencePreset.HUMAN: Organism.HUMAN,
            ReferencePreset.MOUSE: Organism.MOUSE,
            ReferencePreset.RAT: Organism.RAT,
        }[reference_input]
        if dataset_organism is not None and dataset_organism is not target_organism:
            raise ReferenceCompatibilityError(
                "dataset.organism and requested reference preset must match"
            )

    @staticmethod
    def run_bundle_organism(
        *,
        reference_organism: Organism,
        dataset_organism: Organism | None,
        error_type: type[Exception] = ReferenceCompatibilityError,
    ) -> None:
        if dataset_organism is not None and dataset_organism is not reference_organism:
            raise error_type(
                "references.organism must match dataset.organism when both are present"
            )
