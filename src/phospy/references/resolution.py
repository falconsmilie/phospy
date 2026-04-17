"""Reference resolution contracts and bundled reference resolution."""

from __future__ import annotations

from typing import Protocol

from phospy.errors.references import (
    ReferenceCompatibilityError,
    ReferenceResolutionError,
    UnsupportedOrganismError,
)
from phospy.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.references.resources import (
    load_bundled_kinase_substrate_map,
    load_bundled_site_sequences,
)
from phospy.validation.references.bundle import ReferenceBundleValidator


class ReferenceProvider(Protocol):
    """Internal contract for loading a concrete bundle for an organism."""

    def run(self, organism: Organism) -> ReferenceBundle:
        """Return a concrete `ReferenceBundle` for the requested organism."""


class ReferenceResolverContract(Protocol):
    """Internal contract for resolving preset/bundle workflow inputs."""

    def run(
        self,
        reference_input: ReferencePreset | ReferenceBundle,
        *,
        dataset_organism: Organism | None,
    ) -> ReferenceBundle:
        """Resolve workflow reference input into a validated bundle."""


class BundledReferenceProvider:
    """Load packaged bundled reference data for supported organisms."""

    def run(self, organism: Organism) -> ReferenceBundle:
        if not isinstance(organism, Organism):
            raise UnsupportedOrganismError(
                "bundled reference provider requires a supported Organism"
            )
        kinase_substrate_map = load_bundled_kinase_substrate_map(organism)
        site_sequences = load_bundled_site_sequences(organism)
        return ReferenceBundle(
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
        )


class ReferenceResolver:
    """Resolve public reference inputs into concrete validated bundles."""

    def __init__(
        self,
        *,
        provider: ReferenceProvider | None = None,
        validator: ReferenceBundleValidator | None = None,
    ) -> None:
        self._provider = provider or BundledReferenceProvider()
        self._validator = validator or ReferenceBundleValidator()

    def run(
        self,
        reference_input: ReferencePreset | ReferenceBundle,
        *,
        dataset_organism: Organism | None,
    ) -> ReferenceBundle:
        if isinstance(reference_input, ReferenceBundle):
            self._validator.run(
                organism=reference_input.organism,
                kinase_substrate_map=reference_input.kinase_substrate_map,
                site_sequences=reference_input.site_sequences,
                dataset_organism=dataset_organism,
            )
            return reference_input
        if not isinstance(reference_input, ReferencePreset):
            raise ReferenceResolutionError(
                "reference input must be a ReferencePreset or ReferenceBundle"
            )
        organism = self._resolve_target_organism(
            preset=reference_input,
            dataset_organism=dataset_organism,
        )
        bundle = self._provider.run(organism)
        self._validator.run(
            organism=bundle.organism,
            kinase_substrate_map=bundle.kinase_substrate_map,
            site_sequences=bundle.site_sequences,
            dataset_organism=dataset_organism,
        )
        return bundle

    @staticmethod
    def _resolve_target_organism(
        *,
        preset: ReferencePreset,
        dataset_organism: Organism | None,
    ) -> Organism:
        if preset is ReferencePreset.AUTO:
            if dataset_organism is None:
                raise ReferenceResolutionError(
                    "ReferencePreset.AUTO requires dataset.organism"
                )
            return dataset_organism
        mapping = {
            ReferencePreset.HUMAN: Organism.HUMAN,
            ReferencePreset.MOUSE: Organism.MOUSE,
            ReferencePreset.RAT: Organism.RAT,
        }
        target_organism = mapping[preset]
        if dataset_organism is not None and dataset_organism is not target_organism:
            raise ReferenceCompatibilityError(
                "dataset.organism and requested reference preset must match"
            )
        return target_organism
