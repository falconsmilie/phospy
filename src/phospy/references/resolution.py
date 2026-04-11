from __future__ import annotations

from dataclasses import dataclass

from ..errors import InputCompatibilityError
from ..internal.constants import (
    BUNDLED_REFERENCE_ALIASES,
    BUNDLED_REFERENCE_AUTO,
    BUNDLED_REFERENCE_DEFAULTS,
    BUNDLED_REFERENCE_PROVIDER_NAME,
    BUNDLED_REFERENCE_SOURCE,
    BUNDLED_REFERENCE_SPECIES_ALIASES,
    BUNDLED_REFERENCE_VERSION,
)
from .models import (
    ReferenceBundle,
    ReferenceBundleProvenance,
    ReferenceBundleSourceMetadata,
)
from .resources import (
    build_reference_motif_sequences,
    load_bundled_site_sequences,
    load_bundled_substrate_map,
)


@dataclass(frozen=True, slots=True)
class BundledReferenceProvider:
    """Resolve packaged kinase priors for the supported bundled species lanes."""

    source: str = BUNDLED_REFERENCE_SOURCE
    version: str = BUNDLED_REFERENCE_VERSION

    def resolve(
        self,
        *,
        species: str,
        reference: str = BUNDLED_REFERENCE_AUTO,
    ) -> ReferenceBundle:
        resolved_species = normalize_bundled_species(species)
        resolved_reference = normalize_bundled_reference(
            species=resolved_species,
            reference=reference,
        )
        substrate_map = load_bundled_substrate_map(
            species=resolved_species,
            reference=resolved_reference,
        )
        site_sequences = load_bundled_site_sequences(
            species=resolved_species,
            reference=resolved_reference,
        )
        motif_sequences = build_reference_motif_sequences(
            substrate_map=substrate_map,
            site_sequences=site_sequences,
            species=resolved_species,
            reference=resolved_reference,
        )
        return ReferenceBundle(
            substrate_map=substrate_map,
            motif_sequences=motif_sequences,
            species=resolved_species,
            source_metadata=ReferenceBundleSourceMetadata(
                source=self.source,
                reference=resolved_reference,
                version=self.version,
            ),
            provenance=ReferenceBundleProvenance(
                provider=BUNDLED_REFERENCE_PROVIDER_NAME,
                notes=(
                    f"resolved species={resolved_species}",
                    f"resolved reference={resolved_reference}",
                ),
            ),
        )

    @classmethod
    def supported_species(cls) -> tuple[str, ...]:
        return tuple(BUNDLED_REFERENCE_DEFAULTS)

    @classmethod
    def supported_references_for_species(cls, species: str) -> tuple[str, ...]:
        resolved_species = normalize_bundled_species(species)
        canonical_references = {
            resolved_reference
            for resolved_reference in BUNDLED_REFERENCE_ALIASES[
                resolved_species
            ].values()
            if resolved_reference != BUNDLED_REFERENCE_AUTO
        }
        return tuple(sorted(canonical_references))


def normalize_bundled_species(species: str) -> str:
    """Normalize a bundled reference species selector to its canonical key."""

    normalized = str(species).strip().lower()
    resolved_species = BUNDLED_REFERENCE_SPECIES_ALIASES.get(normalized)
    if resolved_species is None:
        supported = ", ".join(sorted(BUNDLED_REFERENCE_DEFAULTS))
        msg = (
            f"Unsupported bundled reference species '{species}'. "
            f"Supported species: {supported}"
        )
        raise InputCompatibilityError(msg)
    return resolved_species


def normalize_bundled_reference(*, species: str, reference: str) -> str:
    """Normalize a bundled reference selector for one canonical species lane."""

    resolved_reference = BUNDLED_REFERENCE_ALIASES[species].get(
        str(reference).strip().lower()
    )
    if resolved_reference is None:
        supported = ", ".join(
            BundledReferenceProvider.supported_references_for_species(species)
        )
        msg = (
            f"Unsupported bundled reference '{reference}' for species '{species}'. "
            f"Supported references: {supported}"
        )
        raise InputCompatibilityError(msg)
    return resolved_reference
