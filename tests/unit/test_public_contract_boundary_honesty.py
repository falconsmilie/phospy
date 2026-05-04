from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import Organism
from phospy.errors import DatasetValidationError
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _coherent_site_identity_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    index = pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id")
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0],
            "sample_b": [3.0, 4.0],
        },
        index=index,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": ["A" * 31, "B" * 31],
        },
        index=index,
    )
    return phospho, site_metadata


def test_dataset_boundary_accepts_coherent_site_identity_rows() -> None:
    phospho, site_metadata = _coherent_site_identity_inputs()

    dataset = AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )

    assert dataset.phospho.index.tolist() == ["MAPK14;Y182;", "AKT1;T308;"]
    assert dataset.site_metadata.index.tolist() == ["MAPK14;Y182;", "AKT1;T308;"]


def test_dataset_boundary_rejects_site_identity_semantic_disagreement_with_details() -> (
    None
):
    phospho, site_metadata = _coherent_site_identity_inputs()
    site_metadata.loc["MAPK14;Y182;", "gene_symbol"] = "MAPK1"
    site_metadata.loc["AKT1;T308;", "site"] = "S473"

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )

    message = str(exc_info.value)
    assert "dataset site-identity coherence failed" in message
    assert (
        "MAPK14;Y182; expected(gene_symbol='MAPK14', site='Y182') "
        "observed(gene_symbol='MAPK1', site='Y182')"
    ) in message
    assert (
        "AKT1;T308; expected(gene_symbol='AKT1', site='T308') "
        "observed(gene_symbol='AKT1', site='S473')"
    ) in message


def test_dataset_boundary_rejects_unparseable_site_ids_before_coherence_checks() -> (
    None
):
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14-Y182"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "site_sequence": ["A" * 31],
        },
        index=phospho.index.copy(),
    )

    with pytest.raises(DatasetValidationError) as exc_info:
        AnalysisReadyPhosphoDataset(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            intensity_scale_state=supported_linear_intensity_scale_state(
                has_total_matrix=False
            ),
            processing_state=supported_linear_processing_state(has_total_matrix=False),
        )

    message = str(exc_info.value)
    assert "site identifiers must use 'GENE;SITE;' format" in message
    assert "'MAPK14-Y182'" in message
