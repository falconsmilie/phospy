from __future__ import annotations

import pandas as pd
import pytest

from phospy.errors import DatasetValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.science.references.models import Organism
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    decode_site_key,
    encode_site_key,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _site_key(
    *,
    organism: object,
    protein_identifier: str,
    site: str,
) -> str:
    key = build_protein_scoped_site_key(
        organism=organism,
        protein_namespace="protein_id",
        protein_identifier=protein_identifier,
        residue=site[0],
        position=int(site[1:]),
        field_name="test.site_key",
        error_type=ValueError,
    )
    return encode_site_key(key)


def _noncanonical_site_key_organism(
    site_key: str,
    *,
    encoded_organism: str,
) -> str:
    return site_key.replace("organism=human", f"organism={encoded_organism}")


def _payload(
    *,
    row_organisms: list[object],
    site_keys: list[str] | None = None,
    organism: Organism | None = None,
) -> dict[str, object]:
    protein_identifiers = ["MAPK14", "AKT1"][: len(row_organisms)]
    sites = ["Y182", "T308"][: len(row_organisms)]
    resolved_site_keys = site_keys or [
        _site_key(
            organism=row_organism,
            protein_identifier=protein_identifier,
            site=site,
        )
        for row_organism, protein_identifier, site in zip(
            row_organisms,
            protein_identifiers,
            sites,
            strict=True,
        )
    ]
    index = pd.Index(resolved_site_keys, name="site_key")
    display_ids = [
        f"{protein_identifier};{site};"
        for protein_identifier, site in zip(protein_identifiers, sites, strict=True)
    ]
    return {
        "phospho": pd.DataFrame(
            {
                "sample_a": [
                    float(row_position + 1)
                    for row_position in range(len(row_organisms))
                ],
                "sample_b": [
                    float(row_position + 2)
                    for row_position in range(len(row_organisms))
                ],
            },
            index=index.copy(),
        ),
        "site_metadata": pd.DataFrame(
            {
                "site_key": resolved_site_keys,
                "display_id": display_ids,
                "organism": row_organisms,
                "protein_namespace": ["protein_id"] * len(row_organisms),
                "protein_identifier": protein_identifiers,
                "gene_symbol": protein_identifiers,
                "site": sites,
                "site_sequence": [("A" * 15) + site[0] + ("A" * 15) for site in sites],
                "protein_id": protein_identifiers,
            },
            index=index.copy(),
        ),
        "organism": organism,
        "intensity_scale_state": supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        "processing_state": supported_linear_processing_state(has_total_matrix=False),
    }


def test_uniform_row_organism_infers_dataset_organism_and_normalizes_aliases() -> None:
    canonical_key = _site_key(
        organism=Organism.HUMAN,
        protein_identifier="MAPK14",
        site="Y182",
    )
    payload = _payload(
        row_organisms=["HUMAN"],
        site_keys=[
            _noncanonical_site_key_organism(
                canonical_key,
                encoded_organism="Homo%20sapiens",
            )
        ],
        organism=None,
    )

    dataset = AnalysisReadyPhosphoDataset(**payload)

    assert dataset.organism is Organism.HUMAN
    assert dataset.phospho.index.tolist() == [canonical_key]
    assert dataset.site_metadata.index.tolist() == [canonical_key]
    assert dataset.site_metadata.loc[:, "site_key"].tolist() == [canonical_key]
    assert dataset.site_metadata.loc[:, "organism"].tolist() == ["human"]


def test_explicit_dataset_organism_must_agree_with_every_row() -> None:
    payload = _payload(row_organisms=["human"], organism=Organism.RAT)

    with pytest.raises(
        DatasetValidationError,
        match="dataset\\.organism must match every .* row_examples",
    ):
        AnalysisReadyPhosphoDataset(**payload)


def test_decoded_site_key_organism_must_agree_with_row_metadata() -> None:
    human_site_key = _site_key(
        organism=Organism.HUMAN,
        protein_identifier="MAPK14",
        site="Y182",
    )
    payload = _payload(
        row_organisms=["rat"],
        site_keys=[human_site_key],
        organism=Organism.RAT,
    )

    with pytest.raises(DatasetValidationError, match="metadata-derived"):
        AnalysisReadyPhosphoDataset(**payload)


def test_mixed_row_organisms_fail_with_actionable_row_identifiers() -> None:
    payload = _payload(row_organisms=["human", "rat"], organism=None)

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(**payload)

    message = str(exc_info.value)
    assert (
        "mixed-organism AnalysisReadyPhosphoDataset rows are not supported" in message
    )
    assert "organism='human'" in message
    assert "organism='rat'" in message
    for site_key in payload["phospho"].index.astype(str).tolist():
        assert site_key in message


def test_trusted_construction_enforces_dataset_organism_coherence() -> None:
    payload = _payload(row_organisms=["human"], organism=Organism.RAT)

    with pytest.raises(DatasetValidationError, match="dataset\\.organism"):
        AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=payload["phospho"],
            site_metadata=payload["site_metadata"],
            organism=payload["organism"],
            intensity_scale_state=payload["intensity_scale_state"],
            processing_state=payload["processing_state"],
        )


def test_decode_site_key_returns_canonical_organism_enum() -> None:
    key = decode_site_key(
        _site_key(organism="homo_sapiens", protein_identifier="MAPK14", site="Y182"),
        field_name="test.site_key",
        error_type=ValueError,
    )

    assert key.organism is Organism.HUMAN
