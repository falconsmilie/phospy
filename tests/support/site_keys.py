from __future__ import annotations

import pandas as pd

from phospy.science.sites.site_keys import (
    ProteinScopedPhosphositeKey,
    decode_site_key,
    encode_site_key,
)


def protein_site_key(
    *,
    protein_identifier: str,
    site: str,
    organism: str = "rat",
    protein_namespace: str = "protein_id",
) -> str:
    return encode_site_key(
        ProteinScopedPhosphositeKey(
            organism=organism,
            protein_namespace=protein_namespace,
            protein_identifier=protein_identifier,
            residue=site[0],
            position=int(site[1:]),
        )
    )


def protein_site_key_index(
    *,
    protein_identifiers: list[str],
    sites: list[str],
    organism: str = "rat",
    protein_namespace: str = "protein_id",
) -> pd.Index:
    return pd.Index(
        [
            protein_site_key(
                protein_identifier=protein_identifier,
                site=site,
                organism=organism,
                protein_namespace=protein_namespace,
            )
            for protein_identifier, site in zip(protein_identifiers, sites, strict=True)
        ],
        name="site_key",
    )


def site_key_from_display_id(
    display_id: str,
    *,
    organism: str = "rat",
    protein_namespace: str = "protein_id",
) -> str:
    protein_identifier, site, *_ = display_id.split(";")
    return protein_site_key(
        protein_identifier=protein_identifier,
        site=site,
        organism=organism,
        protein_namespace=protein_namespace,
    )


def site_key_index_from_display_ids(
    display_ids: list[str],
    *,
    organism: str = "rat",
    protein_namespace: str = "protein_id",
) -> pd.Index:
    return pd.Index(
        [
            site_key_from_display_id(
                display_id,
                organism=organism,
                protein_namespace=protein_namespace,
            )
            for display_id in display_ids
        ],
        name="site_key",
    )


def site_key_context_columns(site_keys: list[str] | pd.Index) -> dict[str, list[str]]:
    decoded_keys = [
        decode_site_key(
            value,
            field_name="tests.support.site_keys.site_key",
            error_type=ValueError,
        )
        for value in pd.Index(site_keys).astype(str).tolist()
    ]
    return {
        "organism": [key.organism.value for key in decoded_keys],
        "protein_namespace": [key.protein_namespace for key in decoded_keys],
        "protein_identifier": [key.protein_identifier for key in decoded_keys],
    }
