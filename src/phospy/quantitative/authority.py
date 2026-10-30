"""Owner-internal authority contract for quantitative state establishment.

This module is package-internal and defines the only approved lanes that can
mint an established intensity scale state.
"""

from __future__ import annotations

from dataclasses import dataclass

from phospy.core.errors.transformations import InvalidTransformationStateError


@dataclass(frozen=True, slots=True, eq=False)
class EstablishmentAuthority:
    """Opaque capability object for approved establishment lanes."""

    source: str


@dataclass(frozen=True, slots=True, eq=False)
class QuantitativeMeaningTransitionAuthority:
    """Opaque capability object for quantitative-meaning transitions."""

    source: str


_DATASET_RESOLVER_AUTHORITY = EstablishmentAuthority(
    source="phospy.science.datasets.builders.transformation_resolver"
)
_SPS_REFERENCE_ESTABLISHMENT_AUTHORITY = EstablishmentAuthority(
    # This historical stored authority value is part of predecessor evidence.
    source="phospy.science.datasets.builders.transformation_resolver"
)
_IDENTITY_TRANSFORMER_AUTHORITY = EstablishmentAuthority(
    source="phospy.science.transformations.transformers.identity"
)
_BUNDLE_RECONSTRUCTION_AUTHORITY = EstablishmentAuthority(
    source="phospy.io.bundles._shared.intensity_scale_state"
)
_DATASET_QUANTITATIVE_MEANING_AUTHORITY = QuantitativeMeaningTransitionAuthority(
    source="phospy.science.datasets.preprocessing.state_builder"
)
_SPS_REFERENCE_QUANTITATIVE_MEANING_AUTHORITY = QuantitativeMeaningTransitionAuthority(
    source="phospy.science.batch_correction.sps_discovery"
)
_BUNDLE_QUANTITATIVE_MEANING_RESTORATION_AUTHORITY = (
    QuantitativeMeaningTransitionAuthority(
        source="phospy.io.bundles._shared.intensity_scale_state"
    )
)


def dataset_resolver_establishment_authority() -> EstablishmentAuthority:
    """Return authority for the dataset intensity-scale resolver lane."""

    return _DATASET_RESOLVER_AUTHORITY


def identity_transformer_establishment_authority() -> EstablishmentAuthority:
    """Return authority for the internal identity-transformer lane."""

    return _IDENTITY_TRANSFORMER_AUTHORITY


def sps_reference_establishment_authority() -> EstablishmentAuthority:
    """Return the dedicated SPS external-assertion establishment capability."""

    return _SPS_REFERENCE_ESTABLISHMENT_AUTHORITY


def bundle_reconstruction_establishment_authority() -> EstablishmentAuthority:
    """Return authority for supported bundle reconstruction lanes."""

    return _BUNDLE_RECONSTRUCTION_AUTHORITY


def dataset_quantitative_meaning_transition_authority() -> (
    QuantitativeMeaningTransitionAuthority
):
    """Return authority for dataset-builder semantic meaning transitions."""

    return _DATASET_QUANTITATIVE_MEANING_AUTHORITY


def bundle_quantitative_meaning_restoration_authority() -> (
    QuantitativeMeaningTransitionAuthority
):
    """Return authority for restoring serialized quantitative meaning provenance."""

    return _BUNDLE_QUANTITATIVE_MEANING_RESTORATION_AUTHORITY


def sps_reference_quantitative_meaning_transition_authority() -> (
    QuantitativeMeaningTransitionAuthority
):
    """Return authority for the supported SPS reference assertion boundary."""

    return _SPS_REFERENCE_QUANTITATIVE_MEANING_AUTHORITY


def resolve_establishment_authority_source(authority: object | None) -> str:
    """Validate authority and return its owning source lane."""

    if authority is _DATASET_RESOLVER_AUTHORITY:
        return _DATASET_RESOLVER_AUTHORITY.source
    if authority is _SPS_REFERENCE_ESTABLISHMENT_AUTHORITY:
        return _SPS_REFERENCE_ESTABLISHMENT_AUTHORITY.source
    if authority is _IDENTITY_TRANSFORMER_AUTHORITY:
        return _IDENTITY_TRANSFORMER_AUTHORITY.source
    if authority is _BUNDLE_RECONSTRUCTION_AUTHORITY:
        return _BUNDLE_RECONSTRUCTION_AUTHORITY.source
    raise InvalidTransformationStateError(
        "intensity scale state can be established only through supported PhosPy "
        "builder/transformer or bundle reconstruction paths"
    )


def resolve_quantitative_meaning_transition_authority_source(
    authority: object | None,
) -> str:
    """Validate quantitative-meaning authority and return its owning lane."""

    if authority is _DATASET_QUANTITATIVE_MEANING_AUTHORITY:
        return _DATASET_QUANTITATIVE_MEANING_AUTHORITY.source
    if authority is _BUNDLE_QUANTITATIVE_MEANING_RESTORATION_AUTHORITY:
        return _BUNDLE_QUANTITATIVE_MEANING_RESTORATION_AUTHORITY.source
    if authority is _SPS_REFERENCE_QUANTITATIVE_MEANING_AUTHORITY:
        return _SPS_REFERENCE_QUANTITATIVE_MEANING_AUTHORITY.source
    raise InvalidTransformationStateError(
        "quantitative meaning can be established or transitioned only through "
        "supported PhosPy dataset-builder or bundle reconstruction paths"
    )


__all__ = [
    "EstablishmentAuthority",
    "QuantitativeMeaningTransitionAuthority",
    "bundle_quantitative_meaning_restoration_authority",
    "bundle_reconstruction_establishment_authority",
    "dataset_quantitative_meaning_transition_authority",
    "dataset_resolver_establishment_authority",
    "identity_transformer_establishment_authority",
    "resolve_establishment_authority_source",
    "resolve_quantitative_meaning_transition_authority_source",
    "sps_reference_establishment_authority",
    "sps_reference_quantitative_meaning_transition_authority",
]
