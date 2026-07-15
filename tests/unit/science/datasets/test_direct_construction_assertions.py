from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import Organism
from phospy.errors import DatasetValidationError
from phospy.provenance import TrustedDatasetConstructionAssertions
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns


def _site_index() -> pd.Index:
    return protein_site_key_index(
        protein_identifiers=["MAPK14", "AKT1"],
        sites=["Y182", "T308"],
    )


def _phospho(index: pd.Index) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [1.5, 2.5],
        },
        index=index.copy(),
    )


def _site_metadata(index: pd.Index) -> pd.DataFrame:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    return pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14", "AKT1"],
            "protein_id": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            ],
            "localisation_confidence": [0.95, 0.9],
        },
        index=index.copy(),
    )


def _trusted_dataset(
    assertions: TrustedDatasetConstructionAssertions | None = None,
) -> AnalysisReadyPhosphoDataset:
    index = _site_index()
    return AnalysisReadyPhosphoDataset.from_trusted_tables(
        phospho=_phospho(index),
        site_metadata=_site_metadata(index),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
        trusted_construction_assertions=assertions,
    )


def test_from_trusted_tables_records_typed_construction_assertions() -> None:
    assertions = TrustedDatasetConstructionAssertions(
        sequence_user_asserted=True,
        identity_user_asserted=True,
        quantitative_meaning_user_asserted=True,
        reference_context_user_asserted=False,
        asserted_by="unit-test",
        assertion_source="curated analysis-ready export",
    )

    dataset = _trusted_dataset(assertions)

    assert dataset.trusted_construction_assertions == assertions
    assert dataset.provenance is not None
    construction = dataset.provenance.workflow_parameters["construction"]
    assert isinstance(construction, dict)
    payload = construction["trusted_construction_assertions"]
    assert isinstance(payload, dict)
    assert payload["assertion_metadata_provided"] is True
    assert payload["sequence_user_asserted"] is True
    assert payload["identity_user_asserted"] is True
    assert payload["quantitative_meaning_user_asserted"] is True
    assert payload["reference_context_user_asserted"] is False
    assert payload["asserted_by"] == "unit-test"
    assert payload["assertion_source"] == "curated analysis-ready export"
    assert construction["missing_trusted_assertions"] == [
        "reference_context_user_asserted"
    ]


def test_from_trusted_tables_records_missing_assertion_metadata_visibly() -> None:
    dataset = _trusted_dataset()

    assertions = dataset.trusted_construction_assertions
    assert assertions is not None
    assert assertions.assertion_metadata_provided is False
    assert assertions.missing_assertions == (
        "sequence_user_asserted",
        "identity_user_asserted",
        "quantitative_meaning_user_asserted",
        "reference_context_user_asserted",
    )

    assert dataset.provenance is not None
    construction = dataset.provenance.workflow_parameters["construction"]
    assert isinstance(construction, dict)
    payload = construction["trusted_construction_assertions"]
    assert isinstance(payload, dict)
    assert payload["assertion_metadata_provided"] is False
    assert payload["sequence_user_asserted"] is False
    assert payload["identity_user_asserted"] is False
    assert payload["quantitative_meaning_user_asserted"] is False
    assert payload["reference_context_user_asserted"] is False
    assert construction["trusted_assertion_metadata_provided"] is False
    assert "assertion_warning" in construction


def test_from_trusted_tables_rejects_untyped_assertion_mapping() -> None:
    index = _site_index()

    with pytest.raises(
        DatasetValidationError,
        match="TrustedDatasetConstructionAssertions",
    ):
        AnalysisReadyPhosphoDataset.from_trusted_tables(
            phospho=_phospho(index),
            site_metadata=_site_metadata(index),
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
            trusted_construction_assertions={
                "sequence_user_asserted": True,
                "identity_user_asserted": True,
                "quantitative_meaning_user_asserted": True,
                "reference_context_user_asserted": True,
            },
        )
