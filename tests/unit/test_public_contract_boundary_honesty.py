from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import Organism
from phospy.errors import DatasetValidationError
from phospy.science.sites.site_keys import (
    build_protein_scoped_site_key,
    encode_site_key,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _encoded_site_key(
    *,
    organism: str = "human",
    protein_namespace: str = "uniprot",
    protein_identifier: str = "P31749",
    residue: str = "T",
    position: int = 308,
) -> str:
    return encode_site_key(
        build_protein_scoped_site_key(
            organism=organism,
            protein_namespace=protein_namespace,
            protein_identifier=protein_identifier,
            residue=residue,
            position=position,
            isoform_id=None,
            field_name="test.site_key",
            error_type=ValueError,
        )
    )


def _coherent_site_identity_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    index = site_key_index_from_display_ids(display_ids)
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "T308"]
            ],
        },
        index=index,
    )
    return phospho, site_metadata


def _single_akt1_t308_inputs(
    *,
    encoded_organism: str = "human",
    encoded_protein_namespace: str = "uniprot",
    encoded_protein_identifier: str = "P31749",
    encoded_residue: str = "T",
    encoded_position: int = 308,
) -> tuple[pd.DataFrame, pd.DataFrame, str, str]:
    encoded_site_key = _encoded_site_key(
        organism=encoded_organism,
        protein_namespace=encoded_protein_namespace,
        protein_identifier=encoded_protein_identifier,
        residue=encoded_residue,
        position=encoded_position,
    )
    metadata_derived_site_key = _encoded_site_key()
    index = pd.Index([encoded_site_key], name="site_key")
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0],
            "sample_b": [2.0],
        },
        index=index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": [encoded_site_key],
            "display_id": ["AKT1;T308;"],
            "organism": ["human"],
            "protein_namespace": ["uniprot"],
            "protein_identifier": ["P31749"],
            "gene_symbol": ["AKT1"],
            "site": ["T308"],
            "site_sequence": [("A" * 15) + "T" + ("A" * 15)],
        },
        index=index.copy(),
    )
    return phospho, site_metadata, encoded_site_key, metadata_derived_site_key


def _construct_analysis_ready_dataset(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
    organism: Organism = Organism.RAT,
) -> AnalysisReadyPhosphoDataset:
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=organism,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def test_dataset_boundary_accepts_coherent_site_identity_rows() -> None:
    phospho, site_metadata = _coherent_site_identity_inputs()

    dataset = _construct_analysis_ready_dataset(
        phospho=phospho,
        site_metadata=site_metadata,
    )

    assert dataset.phospho.index.tolist() == phospho.index.tolist()
    assert dataset.site_metadata.index.tolist() == phospho.index.tolist()
    assert dataset.site_metadata["display_id"].tolist() == [
        "MAPK14;Y182;",
        "AKT1;T308;",
    ]


def test_dataset_boundary_rejects_encoded_t309_when_metadata_says_t308() -> None:
    phospho, site_metadata, encoded_site_key, metadata_derived_site_key = (
        _single_akt1_t308_inputs(encoded_position=309)
    )

    assert phospho.index.tolist() == [encoded_site_key]
    assert site_metadata.index.tolist() == [encoded_site_key]
    assert site_metadata.loc[:, "site_key"].tolist() == [encoded_site_key]

    with pytest.raises(DatasetValidationError) as exc_info:
        _construct_analysis_ready_dataset(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.HUMAN,
        )

    message = str(exc_info.value)
    assert "dataset.site_metadata.site_key" in message
    assert "metadata-derived ProteinScopedPhosphositeKey" in message
    assert "observed=" in message
    assert "expected=" in message
    assert encoded_site_key in message
    assert metadata_derived_site_key in message
    assert "position=309" in message
    assert "position=308" in message


@pytest.mark.parametrize(
    "encoded_overrides",
    [
        pytest.param({"encoded_position": 309}, id="position"),
        pytest.param({"encoded_residue": "S"}, id="residue"),
        pytest.param({"encoded_organism": "mouse"}, id="organism"),
        pytest.param({"encoded_protein_namespace": "refseq"}, id="protein_namespace"),
        pytest.param(
            {"encoded_protein_identifier": "Q9Y243"},
            id="protein_identifier",
        ),
    ],
)
def test_dataset_boundary_rejects_encoded_site_key_metadata_identity_mismatches(
    encoded_overrides: dict[str, object],
) -> None:
    phospho, site_metadata, encoded_site_key, metadata_derived_site_key = (
        _single_akt1_t308_inputs(**encoded_overrides)
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        _construct_analysis_ready_dataset(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.HUMAN,
        )

    message = str(exc_info.value)
    assert "site_key must match metadata-derived" in message
    assert encoded_site_key in message
    assert metadata_derived_site_key in message


@pytest.mark.parametrize(
    ("column_name", "position_value"),
    [
        pytest.param("position", 308, id="position-python-int"),
        pytest.param("position", np.int64(308), id="position-numpy-int"),
        pytest.param("site_position", 308, id="site-position-python-int"),
        pytest.param("site_position", np.int64(308), id="site-position-numpy-int"),
    ],
)
def test_dataset_boundary_accepts_strict_integer_position_metadata(
    column_name: str,
    position_value: object,
) -> None:
    phospho, site_metadata = _coherent_site_identity_inputs()
    site_metadata.loc[:, column_name] = pd.Series(
        [182, position_value],
        index=site_metadata.index,
        dtype="object",
    )

    dataset = _construct_analysis_ready_dataset(
        phospho=phospho,
        site_metadata=site_metadata,
    )

    assert dataset.site_metadata.loc[site_metadata.index[1], column_name] == 308


@pytest.mark.parametrize("column_name", ["position", "site_position"])
@pytest.mark.parametrize(
    "invalid_position",
    [
        pytest.param(True, id="boolean"),
        pytest.param(0, id="zero"),
        pytest.param(-1, id="negative"),
        pytest.param(308.0, id="float"),
        pytest.param("308", id="numeric-string"),
        pytest.param(None, id="missing"),
        pytest.param([308], id="list"),
        pytest.param((308,), id="tuple"),
        pytest.param(np.array([308]), id="array"),
    ],
)
def test_dataset_boundary_rejects_loose_explicit_position_metadata(
    column_name: str,
    invalid_position: object,
) -> None:
    phospho, site_metadata = _coherent_site_identity_inputs()
    site_metadata.loc[:, column_name] = pd.Series(
        [182, invalid_position],
        index=site_metadata.index,
        dtype="object",
    )

    with pytest.raises(DatasetValidationError, match=column_name):
        _construct_analysis_ready_dataset(
            phospho=phospho,
            site_metadata=site_metadata,
        )


def test_dataset_boundary_rejects_site_identity_semantic_disagreement_with_details() -> (
    None
):
    phospho, site_metadata = _coherent_site_identity_inputs()
    mapk14_site_key, akt1_site_key = site_metadata.index.tolist()
    site_metadata.loc[mapk14_site_key, "gene_symbol"] = "MAPK1"
    site_metadata.loc[akt1_site_key, "site"] = "S473"

    with pytest.raises(DatasetValidationError) as exc_info:
        _construct_analysis_ready_dataset(
            phospho=phospho,
            site_metadata=site_metadata,
        )

    message = str(exc_info.value)
    assert (
        "dataset.site_metadata phosphosite identity metadata validation failed"
        in message
    )
    assert (
        "site_sequence central residue must agree with site/residue metadata" in message
    )
    assert str(akt1_site_key) in message


def test_dataset_boundary_rejects_invalid_site_key_before_coherence_checks() -> None:
    index = pd.Index(["MAPK14-Y182"], name="site_key")
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": ["MAPK14-Y182"],
            "display_id": ["MAPK14;Y182;"],
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182"]
            ],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        _construct_analysis_ready_dataset(
            phospho=phospho,
            site_metadata=site_metadata,
        )

    message = str(exc_info.value)
    assert "must contain valid PhosPy site_key values" in message
    assert "must start with 'phospy:v1'" in message
