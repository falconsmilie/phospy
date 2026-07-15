"""Shared reference compatibility validation."""

from __future__ import annotations

from phospy.errors.references import (
    ReferenceCompatibilityError,
    ReferenceResolutionError,
)
from phospy.science.references.models import Organism, ReferenceBundle, ReferencePreset
from phospy.science.sites.organisms import normalize_optional_organism


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
                "dataset.organism and requested reference preset must match"
            )

    def resolve_preset_organism(
        self,
        *,
        preset: ReferencePreset,
        dataset_organism: object | None,
    ) -> Organism:
        """Resolve a compatible preset into the concrete organism to load."""

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
                "references.organism must match dataset.organism when both are present"
            )


def _canonical_dataset_organism(dataset_organism: object | None) -> Organism | None:
    return normalize_optional_organism(
        dataset_organism,
        field_name="dataset.organism",
        error_type=ReferenceResolutionError,
    )
