from __future__ import annotations

from collections.abc import Mapping, Sequence

from ...references import ReferenceBundle
from ..errors import RequestValidationError


def resolve_reference_bundle_inputs(
    *,
    substrate_map: Mapping[str, Sequence[str]] | None,
    motif_sequences: Mapping[str, Sequence[str]] | None,
    reference_bundle: ReferenceBundle | None,
) -> tuple[Mapping[str, Sequence[str]], Mapping[str, Sequence[str]] | None]:
    """Resolve workflow reference inputs from either explicit maps or a bundle."""

    if reference_bundle is None:
        if substrate_map is None:
            msg = "substrate_map must be provided when reference_bundle is not used"
            raise RequestValidationError(msg)
        return substrate_map, motif_sequences

    if substrate_map is not None:
        msg = "Pass either reference_bundle or substrate_map, not both"
        raise RequestValidationError(msg)
    if motif_sequences is not None:
        msg = "Pass either reference_bundle or motif_sequences, not both"
        raise RequestValidationError(msg)
    return reference_bundle.substrate_map, reference_bundle.motif_sequences


__all__ = ["resolve_reference_bundle_inputs"]
