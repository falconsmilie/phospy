"""Reference resolution contracts and bundled reference resolution."""

from __future__ import annotations

from typing import Protocol

from phospy.errors.references import (
    ReferenceCompatibilityError,
    ReferenceResolutionError,
)
from phospy.provenance.hashing import fingerprint_table
from phospy.provenance.models import ReferenceProvenance
from phospy.science.references.models import (
    Organism,
    ReferenceBundle,
    ReferencePreset,
    reference_context_from_manifest_if_complete,
)
from phospy.science.references.resources import (
    bundled_reference_name_for_organism,
    load_bundled_kinase_substrate_map,
    load_bundled_reference_manifest,
    load_bundled_site_sequences,
)
from phospy.science.sites.organisms import normalize_optional_organism


class ReferenceProvider(Protocol):
    """Internal contract for loading a concrete bundle for an organism."""

    def run(self, organism: Organism) -> ReferenceBundle:
        """Return a concrete `ReferenceBundle` for the requested organism."""
        ...


class ReferenceResolverContract(Protocol):
    """Internal contract for resolving preset/bundle workflow inputs."""

    def run(
        self,
        reference_input: ReferencePreset | ReferenceBundle,
        *,
        dataset_organism: Organism | None,
    ) -> ReferenceBundle:
        """Resolve workflow reference input into a validated bundle."""
        ...


class ReferenceCompatibilityValidatorProtocol(Protocol):
    """Reference compatibility checks consumed by reference resolution."""

    def run(
        self,
        reference_input: object,
        *,
        dataset_organism: object | None,
    ) -> None: ...

    def resolve_preset_organism(
        self,
        *,
        preset: ReferencePreset,
        dataset_organism: object | None,
    ) -> Organism: ...


class ReferenceCompatibilityValidator:
    """Validate dataset/reference compatibility for preset or explicit bundle input."""

    _PRESET_TO_ORGANISM = {
        ReferencePreset.HUMAN: Organism.HUMAN,
        ReferencePreset.MOUSE: Organism.MOUSE,
        ReferencePreset.RAT: Organism.RAT,
    }

    def run(
        self,
        reference_input: object,
        *,
        dataset_organism: object | None,
    ) -> None:
        resolved_dataset_organism = _canonical_dataset_organism(dataset_organism)
        if not isinstance(reference_input, (ReferencePreset, ReferenceBundle)):
            raise ReferenceResolutionError(
                "reference input must be a ReferencePreset or ReferenceBundle"
            )
        if isinstance(reference_input, ReferenceBundle):
            self.run_bundle_organism(
                reference_organism=reference_input.organism,
                dataset_organism=resolved_dataset_organism,
            )
            return
        if reference_input is ReferencePreset.AUTO:
            if resolved_dataset_organism is None:
                raise ReferenceResolutionError(
                    "ReferencePreset.AUTO requires dataset.organism"
                )
            return
        target_organism = self._PRESET_TO_ORGANISM[reference_input]
        if (
            resolved_dataset_organism is not None
            and resolved_dataset_organism is not target_organism
        ):
            raise ReferenceCompatibilityError(
                "dataset.organism and requested reference preset must match; "
                f"dataset.organism={resolved_dataset_organism.value!r}; "
                f"references={reference_input.value!r} "
                f"resolved_to={target_organism.value!r}"
            )

    def resolve_preset_organism(
        self,
        *,
        preset: ReferencePreset,
        dataset_organism: object | None,
    ) -> Organism:
        self.run(preset, dataset_organism=dataset_organism)
        resolved_dataset_organism = _canonical_dataset_organism(dataset_organism)
        if preset is ReferencePreset.AUTO:
            if resolved_dataset_organism is None:
                raise ReferenceResolutionError(
                    "ReferencePreset.AUTO requires dataset.organism"
                )
            return resolved_dataset_organism
        return self._PRESET_TO_ORGANISM[preset]

    @staticmethod
    def run_bundle_organism(
        *,
        reference_organism: Organism,
        dataset_organism: object | None,
        error_type: type[Exception] = ReferenceCompatibilityError,
    ) -> None:
        resolved_dataset_organism = normalize_optional_organism(
            dataset_organism,
            field_name="dataset.organism",
            error_type=error_type,
        )
        if (
            resolved_dataset_organism is not None
            and resolved_dataset_organism is not reference_organism
        ):
            raise error_type(
                "references.organism must match dataset.organism when both are "
                "present; "
                f"dataset.organism={resolved_dataset_organism.value!r}; "
                f"references.organism={reference_organism.value!r}"
            )


def _canonical_dataset_organism(dataset_organism: object | None) -> Organism | None:
    return normalize_optional_organism(
        dataset_organism,
        field_name="dataset.organism",
        error_type=ReferenceResolutionError,
    )


class BundledReferenceProvider:
    """Load packaged bundled reference data for supported organisms.

    In the current cutover release, bundled runtime coverage is rat-only.
    """

    def run(self, organism: Organism) -> ReferenceBundle:
        bundle_id = bundled_reference_name_for_organism(organism)
        manifest = load_bundled_reference_manifest(organism)
        if manifest.bundle_id != bundle_id:
            raise ReferenceResolutionError(
                "bundled reference manifest bundle_id does not match bundled lane name "
                f"for organism '{organism.value}': expected '{bundle_id}', "
                f"got '{manifest.bundle_id}'"
            )
        source_version = manifest.source_version
        if source_version is None:
            raise ReferenceResolutionError(
                "bundled reference manifest is missing source_version for "
                f"{organism.value}/{manifest.bundle_id}"
            )
        kinase_substrate_map = load_bundled_kinase_substrate_map(organism)
        site_sequences = load_bundled_site_sequences(organism)
        provenance = ReferenceProvenance(
            source_type="bundled",
            organism=organism.value,
            bundle_id=bundle_id,
            source_name=manifest.source_name,
            source_version=source_version,
            retrieved_at=manifest.retrieved_at.isoformat(),
            identifier_namespace=manifest.identifier_namespace,
            sequence_window=manifest.sequence_window.to_payload(),
            manifest=manifest.to_payload(),
            table_fingerprints=(
                fingerprint_table(
                    kinase_substrate_map,
                    name="references.kinase_substrate_map",
                ),
                fingerprint_table(
                    site_sequences,
                    name="references.site_sequences",
                ),
            ),
            reference_context=reference_context_from_manifest_if_complete(manifest),
        )
        return ReferenceBundle._from_owned(  # pyright: ignore[reportPrivateUsage] - trusted internal constructor avoids redundant frame copies
            organism=organism,
            kinase_substrate_map=kinase_substrate_map,
            site_sequences=site_sequences,
            provenance=provenance,
            manifest=manifest,
        )


class ReferenceResolver:
    """Resolve public reference inputs into concrete validated bundles."""

    def __init__(
        self,
        *,
        provider: ReferenceProvider | None = None,
        compatibility_validator: ReferenceCompatibilityValidatorProtocol | None = None,
    ) -> None:
        self._provider = provider or BundledReferenceProvider()
        self._compatibility_validator = (
            compatibility_validator or ReferenceCompatibilityValidator()
        )

    def run(
        self,
        reference_input: ReferencePreset | ReferenceBundle,
        *,
        dataset_organism: Organism | None,
    ) -> ReferenceBundle:
        if isinstance(reference_input, ReferenceBundle):
            # Structural bundle validation is owned by ReferenceBundle construction.
            self._compatibility_validator.run(
                reference_input,
                dataset_organism=dataset_organism,
            )
            return reference_input
        if not isinstance(reference_input, ReferencePreset):  # pyright: ignore[reportUnnecessaryIsInstance] - runtime guard for untyped/boundary inputs
            raise ReferenceResolutionError(
                "reference input must be a ReferencePreset or ReferenceBundle; "
                f"got {type(reference_input).__name__}"
            )
        organism = self._compatibility_validator.resolve_preset_organism(
            preset=reference_input,
            dataset_organism=dataset_organism,
        )
        return self._provider.run(organism)
