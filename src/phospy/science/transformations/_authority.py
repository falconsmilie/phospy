"""Internal authority capabilities for intensity-scale-state establishment.

This module is package-internal and defines the only approved lanes that can
mint an established intensity scale state.
"""

from __future__ import annotations

from dataclasses import dataclass

from phospy.errors.transformations import InvalidTransformationStateError


@dataclass(frozen=True, slots=True, eq=False)
class _EstablishmentAuthority:
    """Opaque capability object for approved establishment lanes."""

    source: str


@dataclass(frozen=True, slots=True, eq=False)
class _QuantitativeMeaningTransitionAuthority:
    """Opaque capability object for quantitative-meaning transitions."""

    source: str


_DATASET_RESOLVER_AUTHORITY = _EstablishmentAuthority(
    source="phospy.science.datasets.builders.transformation_resolver"
)
_IDENTITY_TRANSFORMER_AUTHORITY = _EstablishmentAuthority(
    source="phospy.science.transformations.transformers.identity"
)
_BUNDLE_RECONSTRUCTION_AUTHORITY = _EstablishmentAuthority(
    source="phospy.io.bundles._shared.intensity_scale_state"
)
_DATASET_QUANTITATIVE_MEANING_AUTHORITY = _QuantitativeMeaningTransitionAuthority(
    source="phospy.science.datasets.preprocessing.state_builder"
)
_BUNDLE_QUANTITATIVE_MEANING_RESTORATION_AUTHORITY = (
    _QuantitativeMeaningTransitionAuthority(
        source="phospy.io.bundles._shared.intensity_scale_state"
    )
)


def _dataset_resolver_establishment_authority() -> _EstablishmentAuthority:
    """Return authority for the dataset intensity-scale resolver lane."""

    return _DATASET_RESOLVER_AUTHORITY


def _identity_transformer_establishment_authority() -> _EstablishmentAuthority:
    """Return authority for the internal identity-transformer lane."""

    return _IDENTITY_TRANSFORMER_AUTHORITY


def _bundle_reconstruction_establishment_authority() -> _EstablishmentAuthority:
    """Return authority for supported bundle reconstruction lanes."""

    return _BUNDLE_RECONSTRUCTION_AUTHORITY


def _dataset_quantitative_meaning_transition_authority() -> (
    _QuantitativeMeaningTransitionAuthority
):
    """Return authority for dataset-builder semantic meaning transitions."""

    return _DATASET_QUANTITATIVE_MEANING_AUTHORITY


def _bundle_quantitative_meaning_restoration_authority() -> (
    _QuantitativeMeaningTransitionAuthority
):
    """Return authority for restoring serialized quantitative meaning provenance."""

    return _BUNDLE_QUANTITATIVE_MEANING_RESTORATION_AUTHORITY


def _resolve_establishment_authority_source(authority: object | None) -> str:
    """Validate authority and return its owning source lane."""

    if authority is _DATASET_RESOLVER_AUTHORITY:
        return _DATASET_RESOLVER_AUTHORITY.source
    if authority is _IDENTITY_TRANSFORMER_AUTHORITY:
        return _IDENTITY_TRANSFORMER_AUTHORITY.source
    if authority is _BUNDLE_RECONSTRUCTION_AUTHORITY:
        return _BUNDLE_RECONSTRUCTION_AUTHORITY.source
    raise InvalidTransformationStateError(
        "intensity scale state can be established only through supported PhosPy "
        "builder/transformer or bundle reconstruction paths"
    )


def _resolve_quantitative_meaning_transition_authority_source(
    authority: object | None,
) -> str:
    """Validate quantitative-meaning authority and return its owning lane."""

    if authority is _DATASET_QUANTITATIVE_MEANING_AUTHORITY:
        return _DATASET_QUANTITATIVE_MEANING_AUTHORITY.source
    if authority is _BUNDLE_QUANTITATIVE_MEANING_RESTORATION_AUTHORITY:
        return _BUNDLE_QUANTITATIVE_MEANING_RESTORATION_AUTHORITY.source
    raise InvalidTransformationStateError(
        "quantitative meaning can be established or transitioned only through "
        "supported PhosPy dataset-builder or bundle reconstruction paths"
    )


# Public aliases for cross-module internal typing/authority checks.
EstablishmentAuthority = _EstablishmentAuthority
QuantitativeMeaningTransitionAuthority = _QuantitativeMeaningTransitionAuthority
dataset_resolver_establishment_authority = _dataset_resolver_establishment_authority
identity_transformer_establishment_authority = (
    _identity_transformer_establishment_authority
)
bundle_reconstruction_establishment_authority = (
    _bundle_reconstruction_establishment_authority
)
dataset_quantitative_meaning_transition_authority = (
    _dataset_quantitative_meaning_transition_authority
)
bundle_quantitative_meaning_restoration_authority = (
    _bundle_quantitative_meaning_restoration_authority
)
resolve_establishment_authority_source = _resolve_establishment_authority_source
resolve_quantitative_meaning_transition_authority_source = (
    _resolve_quantitative_meaning_transition_authority_source
)
