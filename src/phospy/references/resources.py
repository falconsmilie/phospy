from __future__ import annotations

from collections.abc import Mapping, Sequence
from importlib import resources
from pathlib import Path

import pandas as pd

from ..validation.errors import InputCompatibilityError


def bundled_reference_resource_path(
    *,
    species: str,
    reference: str,
    filename: str,
) -> Path:
    """Resolve one packaged bundled-reference asset to a local path."""

    resource = resources.files("phospy").joinpath(
        "data",
        "reference_bundles",
        species,
        reference,
        filename,
    )
    if not resource.is_file():
        msg = (
            "BundledReferenceProvider could not find packaged reference data for "
            f"species '{species}' and reference '{reference}' ({filename})"
        )
        raise InputCompatibilityError(msg)

    with resources.as_file(resource) as resolved_path:
        return resolved_path


def load_grouped_mapping_file(
    path: Path,
    *,
    group_column: str,
    value_column: str,
) -> dict[str, tuple[str, ...]]:
    """Load a two-column grouped mapping file into a kinase-keyed mapping."""

    frame = pd.read_csv(path)
    grouped: dict[str, list[str]] = {}
    for group, value in frame.loc[:, [group_column, value_column]].itertuples(
        index=False
    ):
        grouped.setdefault(str(group).strip(), []).append(str(value).strip())
    return {key: tuple(values) for key, values in grouped.items()}


def load_string_mapping_file(
    path: Path,
    *,
    key_column: str,
    value_column: str,
) -> dict[str, str]:
    """Load a two-column string mapping file."""

    frame = pd.read_csv(path)
    return {
        str(key).strip(): str(value).strip()
        for key, value in frame.loc[:, [key_column, value_column]].itertuples(
            index=False
        )
    }


def load_bundled_substrate_map(
    *,
    species: str,
    reference: str,
) -> dict[str, tuple[str, ...]]:
    """Load the packaged substrate map for one bundled species/reference lane."""

    return load_grouped_mapping_file(
        bundled_reference_resource_path(
            species=species,
            reference=reference,
            filename="substrate_map.csv",
        ),
        group_column="kinase",
        value_column="site_id",
    )


def load_bundled_site_sequences(
    *,
    species: str,
    reference: str,
) -> dict[str, str]:
    """Load the packaged site-sequence table for one bundled species/reference lane."""

    return load_string_mapping_file(
        bundled_reference_resource_path(
            species=species,
            reference=reference,
            filename="site_sequences.csv",
        ),
        key_column="site_id",
        value_column="centralized_sequence",
    )


def build_reference_motif_sequences(
    *,
    substrate_map: Mapping[str, Sequence[str]],
    site_sequences: Mapping[str, str],
    species: str,
    reference: str,
) -> dict[str, tuple[str, ...]]:
    """Build bundled motif-sequence inputs from substrate and site resources."""

    motif_sequences: dict[str, tuple[str, ...]] = {}
    missing_sites: set[str] = set()
    for kinase, site_ids in substrate_map.items():
        sequences: list[str] = []
        for site_id in site_ids:
            sequence = site_sequences.get(str(site_id))
            if sequence is None:
                missing_sites.add(str(site_id))
                continue
            sequences.append(str(sequence))
        motif_sequences[str(kinase)] = tuple(sequences)

    if missing_sites:
        missing_text = ", ".join(sorted(missing_sites))
        msg = (
            "BundledReferenceProvider reference data is incomplete for "
            f"species '{species}' and reference '{reference}'; "
            f"missing site sequences for: {missing_text}"
        )
        raise InputCompatibilityError(msg)

    return motif_sequences
