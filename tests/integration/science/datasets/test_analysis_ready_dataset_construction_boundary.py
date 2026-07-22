from __future__ import annotations

from collections.abc import Mapping

import pandas as pd
import pytest

from phospy import AnalysisReadyDatasetBuilder
from phospy.api import DatasetBuildRequest, Organism
from phospy.errors import DatasetValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import protein_site_key_index, site_key_context_columns

pytestmark = pytest.mark.integration

_CENTRED_Y_SEQUENCE = "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"


def _direct_constructor_payload() -> dict[str, object]:
    site_index = protein_site_key_index(
        protein_identifiers=["MAPK14"],
        sites=["Y182"],
    )
    return {
        "phospho": pd.DataFrame(
            {"sample_a": [1.0], "sample_b": [2.0]},
            index=site_index.copy(),
        ),
        "site_metadata": pd.DataFrame(
            {
                "site_key": site_index.tolist(),
                "display_id": ["MAPK14;Y182;"],
                **site_key_context_columns(site_index),
                "gene_symbol": ["MAPK14"],
                "site": ["Y182"],
                "protein_id": ["MAPK14"],
                "site_sequence": [_CENTRED_Y_SEQUENCE],
            },
            index=site_index.copy(),
        ),
        "sample_metadata": pd.DataFrame(
            {"condition": ["control", "treated"]},
            index=pd.Index(["sample_a", "sample_b"], name="sample_id"),
        ),
        "organism": Organism.RAT,
        "intensity_scale_state": supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        "processing_state": supported_linear_processing_state(has_total_matrix=False),
    }


def test_supported_construction_boundary_uses_builder_path_with_provenance() -> None:
    phospho = pd.DataFrame(
        {"sample_a": [1.0], "sample_b": [2.0]},
        index=pd.Index(["MAPK14;Y182;"], name="site_id"),
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14"],
            "site": ["Y182"],
            "protein_id": ["MAPK14"],
            "site_sequence": [_CENTRED_Y_SEQUENCE],
            "localisation_confidence": [0.95],
        },
        index=phospho.index.copy(),
    )

    dataset = AnalysisReadyDatasetBuilder().run(
        DatasetBuildRequest(
            phospho=phospho,
            site_metadata=site_metadata,
            organism=Organism.RAT,
            input_intensity_scale="linear",
        )
    )

    assert isinstance(dataset, AnalysisReadyPhosphoDataset)
    assert dataset.provenance is not None
    construction = dataset.provenance.workflow_parameters.get("construction")
    assert isinstance(construction, Mapping)
    assert construction["method"] == "AnalysisReadyDatasetBuilder.run"
    assert dataset.phospho.index.name == "site_key"
    assert dataset.site_metadata.index.equals(dataset.phospho.index)


def test_direct_constructor_fails_immediately_with_supported_paths() -> None:
    payload = _direct_constructor_payload()

    with pytest.raises(
        TypeError,
        match="AnalysisReadyDatasetBuilder.*from_trusted_tables",
    ):
        AnalysisReadyPhosphoDataset(**payload)


def test_trusted_factory_still_rejects_invalid_tables() -> None:
    payload = _direct_constructor_payload()
    site_metadata = payload["site_metadata"]
    assert isinstance(site_metadata, pd.DataFrame)
    payload["site_metadata"] = site_metadata.drop(columns=["site_sequence"])

    with pytest.raises(DatasetValidationError, match="site_sequence"):
        trusted_analysis_ready_dataset_from_tables(**payload)


def test_trusted_factory_rejects_invalid_nested_processing_state() -> None:
    payload = _direct_constructor_payload()
    processing_state = payload["processing_state"]
    object.__setattr__(processing_state, "ruv_readiness", object())

    with pytest.raises(DatasetValidationError, match="ruv_readiness"):
        trusted_analysis_ready_dataset_from_tables(**payload)
