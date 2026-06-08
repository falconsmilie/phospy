from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import Organism
from phospy.errors import DatasetValidationError
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
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


def _construct_analysis_ready_dataset(
    *,
    phospho: pd.DataFrame,
    site_metadata: pd.DataFrame,
) -> AnalysisReadyPhosphoDataset:
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
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
