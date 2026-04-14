from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from ..errors import InputCompatibilityError
from ..validation.values.identifiers import parse_canonical_site_id

__all__ = [
    "parse_supported_signalome_site_ids",
    "protein_id_from_supported_signalome_site_id",
    "resolve_signalome_site_to_protein",
]


def protein_id_from_supported_signalome_site_id(site_id: object) -> str | None:
    """Extract the protein identifier from a supported signalome site ID.

    Supported IDs use canonical ``ENTITY;SITE;`` shape where the site token
    contains at least one letter followed by at least one digit, such as
    ``S1`` or ``T308``.
    """

    parsed = parse_canonical_site_id(site_id)
    if parsed is None:
        return None
    entity, _ = parsed
    return entity


def parse_supported_signalome_site_ids(
    site_ids: Sequence[str],
    *,
    invalid_site_id_context: str,
) -> pd.Series:
    """Parse supported phosphosite identifiers into protein IDs."""

    protein_ids: list[str] = []
    invalid_site_ids: list[str] = []
    for site_id in site_ids:
        protein_id = protein_id_from_supported_signalome_site_id(site_id)
        if protein_id is None:
            invalid_site_ids.append(str(site_id))
            continue
        protein_ids.append(protein_id)

    if invalid_site_ids:
        preview = ", ".join(invalid_site_ids[:3])
        msg = f"{invalid_site_id_context}: {preview}"
        if len(invalid_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    return _build_site_to_protein_series(site_ids=site_ids, protein_ids=protein_ids)


def resolve_signalome_site_to_protein(
    *,
    site_ids: Sequence[str],
    site_to_protein: Mapping[str, str] | pd.Series | None,
    missing_mapping_context: str,
    invalid_mapping_context: str,
    invalid_site_id_context: str,
) -> pd.Series:
    """Resolve aligned phosphosite IDs to validated protein IDs."""

    if site_to_protein is None:
        return parse_supported_signalome_site_ids(
            site_ids,
            invalid_site_id_context=invalid_site_id_context,
        )

    if isinstance(site_to_protein, pd.Series):
        mapping: Mapping[str, str] = {
            str(site_id): str(protein_id)
            for site_id, protein_id in site_to_protein.items()
        }
    else:
        mapping = site_to_protein

    missing_site_ids = [str(site_id) for site_id in site_ids if site_id not in mapping]
    if missing_site_ids:
        preview = ", ".join(missing_site_ids[:3])
        msg = f"{missing_mapping_context}: {preview}"
        if len(missing_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    protein_ids = [str(mapping[site_id]).strip() for site_id in site_ids]
    invalid_site_ids = [
        str(site_id)
        for site_id, protein_id in zip(site_ids, protein_ids, strict=True)
        if not protein_id
    ]
    if invalid_site_ids:
        preview = ", ".join(invalid_site_ids[:3])
        msg = f"{invalid_mapping_context}: {preview}"
        if len(invalid_site_ids) > 3:
            msg += ", ..."
        raise InputCompatibilityError(msg)

    return _build_site_to_protein_series(site_ids=site_ids, protein_ids=protein_ids)


def _build_site_to_protein_series(
    *,
    site_ids: Sequence[str],
    protein_ids: Sequence[str],
) -> pd.Series:
    series = pd.Series(
        [str(protein_id) for protein_id in protein_ids],
        index=pd.Index([str(site_id) for site_id in site_ids], dtype=object),
        dtype=object,
    )
    series.index.name = "site_id"
    series.name = "protein_id"
    return series
